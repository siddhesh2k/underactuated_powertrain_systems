%% setup_HGa_nonlinear.m
%

clear; clc; close all;

%% ========================================================================
%  1. Run the H-infinity design script to get K, Ga, Tpwr1
% ========================================================================
fprintf('\n=== Running Hinf4_null_allocator.m ===\n');
run('Hinf4_null_allocator.m');

% Verify required variables
needed = {'Tpwr1','Ktfpwr1','Gppwr2','Ga'};
for k = 1:numel(needed)
    assert(exist(needed{k},'var') == 1, ...
        '%s not found. Check Hinf4_null_allocator.m ran correctly.', needed{k});
end
fprintf('Tpwr1, Ktfpwr1, Gppwr2, Ga all present.\n');

%% ========================================================================
%  2. Run nonlinear plant parameter file
% ========================================================================
fprintf('\n=== Running parameter file ===\n');
% Save controller variables before the parameter script's clear all
Tpwr1_save   = Tpwr1;
Ktfpwr1_save = Ktfpwr1;
Gppwr2_save  = Gppwr2;
Ga_save      = Ga;

run('PTsystem_Openloop_simulink_parameter_v03.m');

% Restore controller variables
Tpwr1   = Tpwr1_save;
Ktfpwr1 = Ktfpwr1_save;
Gppwr2  = Gppwr2_save;
Ga      = Ga_save;
clear Tpwr1_save Ktfpwr1_save Gppwr2_save Ga_save;

%% ========================================================================
%  3. Load reference/demand data
% ========================================================================
fprintf('\n=== Loading data ===\n');

data_file = 'FCHEV_SDP_target_pdem_data_UDDS.mat';
D = load(data_file);

Ts   = 1;
SOC0 = 0.6;

SOC_ref = D.SOC_Target(2,:).';
H2_ref  = D.H2_cons_target(2,:).';
P_dem   = D.P_dem_simulink(2,:).';

N = min([numel(SOC_ref), numel(H2_ref), numel(P_dem)]);
SOC_ref = SOC_ref(1:N);
H2_ref  = H2_ref(1:N);
P_dem   = P_dem(1:N);

t = (0:N-1).' * Ts;
Tsim = t(end);

SOC_ref_delta = SOC_ref - SOC0;

% From Workspace signals [time, value]
SOC_ref_delta_fw = [t, SOC_ref_delta];
SOC_ref_fw       = [t, SOC_ref];
H2_ref_fw        = [t, H2_ref];
Pdem_fw          = [t, P_dem];

fprintf('  N = %d, Ts = %g s, Tsim = %g s\n', N, Ts, Tsim);

%% ========================================================================
%  4. Build H prefilter at omega0 = 0.04 rad/s
% ========================================================================
fprintf('\n=== Building H prefilter ===\n');

omega0 = 0.04;
Ty_w0  = evalfr(Tpwr1, 1j*omega0);
Ty_mag = abs(Ty_w0);
det_T_mag = Ty_mag(1,1)*Ty_mag(2,2) - Ty_mag(1,2)*Ty_mag(2,1);

H_mat = [ Ty_mag(1,1)*Ty_mag(2,2)/det_T_mag,  0;
         -Ty_mag(1,1)*Ty_mag(2,1)/det_T_mag,  1];

fprintf('  H_mat = [%.4e, 0; %.4f, 1]\n', H_mat(1,1), H_mat(2,1));
fprintf('  omega0 = %.4f rad/s\n', omega0);

%% ========================================================================
%  5. Discretize Ga(1,1) and K for state-space blocks
% ========================================================================
fprintf('\n=== Discretizing Ga(1,1) and K ===\n');

% Ga(1,1): the null-space allocator prefilter on SOC ref
Ga11_d = c2d(ss(Ga(1,1)), Ts, 'tustin');
[Aga, Bga, Cga, Dga] = ssdata(Ga11_d);
fprintf('  Ga11: order=%d, DC gain=%.1f dB\n', size(Aga,1), ...
    20*log10(max(abs(dcgain(Ga(1,1))),eps)));

% K: the H-infinity controller (1x2: [e_soc; e_h2] -> Pfc)
K_d = c2d(ss(Ktfpwr1), Ts, 'tustin');
[Ak_ctrl, Bk_ctrl, Ck_ctrl, Dk_ctrl] = ssdata(K_d);
fprintf('  K: order=%d, size=%dx%d\n', size(Ak_ctrl,1), size(Ck_ctrl,1), size(Bk_ctrl,2));

%% ========================================================================
%  6. Saturation limits
% ========================================================================
Pfc_min = 0;
Pfc_max = 106e3;


%% PAPER FIGURE from Simulink run: SOC tracking / M_H2 tracking / errors
%% =========================================================================
% Signal map (from the Simulink model):
%   out.SOC_abs_sim -> SOC reference   (SOC_T)
%   out.simout1     -> M_H2 reference  (M_H2,T)
%   out.simout      -> SOC output
%   out.H2_sim      -> M_H2 output

[t,       SOCref_sim] = local_getsig(out.SOC_abs_sim);
[~,       H2ref_sim ] = local_getsig(out.simout1);
[t_soc,   SOC_out   ] = local_getsig(out.simout);
[t_h2,    Mh2_out   ] = local_getsig(out.H2_sim);

% Put everything on the SOC-output time base (handles variable-step solvers)
if numel(t) ~= numel(t_soc) || any(abs(t(:) - t_soc(:)) > 1e-9)
    SOCref_sim = interp1(t, SOCref_sim, t_soc, 'linear', 'extrap');
end
[~, H2ref_t] = local_getsig(out.simout1);
if numel(H2ref_t) ~= numel(t_soc)
    % re-fetch its own time vector for interpolation
end
th2r = local_gettime(out.simout1);
if numel(th2r) ~= numel(t_soc) || any(abs(th2r(:) - t_soc(:)) > 1e-9)
    H2ref_sim = interp1(th2r, H2ref_sim, t_soc, 'linear', 'extrap');
end
if numel(t_h2) ~= numel(t_soc) || any(abs(t_h2(:) - t_soc(:)) > 1e-9)
    Mh2_out = interp1(t_h2, Mh2_out, t_soc, 'linear', 'extrap');
end
t = t_soc(:);

SOC_out     = SOC_out(:);      Mh2_out    = Mh2_out(:);
SOCref_sim  = SOCref_sim(:);   H2ref_sim  = H2ref_sim(:);

SOC_err_fig = SOC_out - SOCref_sim;
H2_err_fig  = Mh2_out - H2ref_sim;

fprintf('\n=== Nonlinear (Simulink) results ===\n');
fprintf('  SOC RMS error          = %.4e\n',  rms(SOC_err_fig));
fprintf('  SOC max |err|          = %.4e\n',  max(abs(SOC_err_fig)));
fprintf('  SOC end-of-cycle error = %+.4e\n', SOC_err_fig(end));
fprintf('  M_H2 RMS error         = %.4e\n',  rms(H2_err_fig));
fprintf('  M_H2 max |err|         = %.4e\n',  max(abs(H2_err_fig)));
fprintf('  M_H2 end-of-cycle err  = %+.4e\n', H2_err_fig(end));

figpaper = figure('Name','Paper figure: nonlinear tracking','Color','w', ...
    'Position',[100 100 700 750]);
tlp = tiledlayout(figpaper,3,1,'Padding','compact','TileSpacing','compact');

axp = nexttile(tlp);
plot(axp, t, SOC_out, 'b-', 'LineWidth', 1.3); hold(axp,'on');
plot(axp, t, SOCref_sim, 'r--', 'LineWidth', 1.3);
grid(axp,'on'); ylabel(axp,'$SOC$','Interpreter','latex');
legend(axp,{'$SOC$','$SOC_T$'},'Interpreter','latex','Location','best');
title(axp,'(a) SOC tracking');

axp = nexttile(tlp);
plot(axp, t, Mh2_out, 'b-', 'LineWidth', 1.3); hold(axp,'on');
plot(axp, t, H2ref_sim, 'r--', 'LineWidth', 1.3);
grid(axp,'on'); ylabel(axp,'$M_{H_2}$ (g)','Interpreter','latex');
legend(axp,{'$M_{H_2}$','$M_{H_2,T}$'},'Interpreter','latex','Location','best');
title(axp,'(b) Hydrogen consumption tracking');

axp = nexttile(tlp);
yyaxis(axp,'left');
plot(axp, t, SOC_err_fig, 'b-', 'LineWidth', 1.2);
ylabel(axp,'$SOC - SOC_T$','Interpreter','latex');
yyaxis(axp,'right');
plot(axp, t, H2_err_fig, 'r-', 'LineWidth', 1.2);
ylabel(axp,'$M_{H_2} - M_{H_2,T}$ (g)','Interpreter','latex');
grid(axp,'on'); xlabel(axp,'Time (s)');
legend(axp,{'SOC error (left)','$M_{H_2}$ error (right)'}, ...
    'Interpreter','latex','Location','best');
title(axp,'(c) Tracking errors');

set(findall(figpaper,'-property','FontSize'),'FontSize',11);
exportgraphics(figpaper,'nonlinear_SensorAlloc.png','Resolution',300);

function [tv, yv] = local_getsig(sig)
if isa(sig,'timeseries')
    tv = sig.Time(:);  yv = squeeze(sig.Data);
elseif isstruct(sig) && isfield(sig,'signals')      % struct-with-time
    tv = sig.time(:);  yv = squeeze(sig.signals.values);
elseif isa(sig,'Simulink.SimulationData.Signal')
    tv = sig.Values.Time(:); yv = squeeze(sig.Values.Data);
else                                                 % plain array
    yv = squeeze(sig);  tv = (0:numel(yv)-1).';
end
yv = yv(:);
end

function tv = local_gettime(sig)
[tv, ~] = local_getsig(sig);
end