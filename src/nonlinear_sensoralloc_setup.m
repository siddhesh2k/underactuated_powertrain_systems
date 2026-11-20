%% create_sensor_allocator_hinf_simulink_model.m
% Creates a Simulink model for:
%   Sensor-allocation control + SISO H-infinity controller + SOC integral bias
%
% 
clear; clc; close all;

%% ========================================================================
%  Run parameter file FIRST (it has clear all)
% ========================================================================
run('PTsystem_Openloop_simulink_parameter_v03.m');

%% ========================================================================
%  USER SETTINGS (defined AFTER the parameter file's clear all)
% ========================================================================
model_name = 'nonlinear_SensorAlloc';
data_file  = 'FCHEV_SDP_target_pdem_data_UDDS.mat'; %nonlinear refs
% load('FCHEV_SDP_Targets_linear_UDDS.mat');
Ts   = 1;
SOC0 = 0.6;

N = 1370;

% Engineering normalization from your script
SOC_scale = 0.01;
H2_scale  = 6.0596e3;

% Allocation settings
f_alloc = 0.3728;
w1 = 1;
w2 = 1;

% H-infinity synthesis settings from your best sweep
wc_hz   = 0.15;
M_w1    = 1.2;
A_w1    = 0.01;
W2scale = 0.05;

% Power limits and soft floor
Pfc_min = 0;
Pfc_max = 106e3;
Pfc_soft_floor = 0;

% SOC integral bias
% In the linear sim, KsocI = -2e6 works because the algebraic-loop solver
% inverts the sign path internally. In Simulink with direct error wiring:
%   e1 > 0 (SOC too low) -> Isoc accumulates positive -> bias should
%   INCREASE Pfc (to raise SOC) -> need KsocI > 0.
KsocI = 2e6;
leak  = 0.002;

%% ========================================================================
%  STEP 1: Load data
% ========================================================================

fprintf('\n=== Loading data ===\n');

if exist(data_file, 'file') ~= 2
    warning('Data file not found. Creating dummy test signals.');
    N = 200;
    t = (0:N-1).' * Ts;
    SOCref = SOC0 + 5e-4*sin(2*pi*t/80);
    H2ref  = 0.04*t;
    Pdem   = 20e3 + 10e3*sin(2*pi*t/50);
else
    D = load(data_file);

    required_vars = {'SOC_Target','H2_cons_target','P_dem_simulink'};
    for k = 1:numel(required_vars)
        if ~isfield(D, required_vars{k})
            error('Missing variable %s in %s.', required_vars{k}, data_file);
        end
    end

    Pdem   = D.P_dem_simulink(2,:).';
    SOCref = D.SOC_Target(2,:).';
    H2ref  = D.H2_cons_target(2,:).';
   

    N = min([N, numel(Pdem), numel(SOCref), numel(H2ref)]);
    Pdem   = Pdem(1:N);
    SOCref = SOCref(1:N);
    H2ref  = H2ref(1:N);
    t = (0:N-1).' * Ts;
end

Tsim = t(end);

% Simulink From Workspace signals use [time, value]
SOC_ref_fw = [t SOCref];
H2_ref_fw  = [t H2ref];
Pdem_fw    = [t Pdem];
% SOC_ref_fw = SOC_Target;
% H2_ref_fw  = H2_cons_target;
% Pdem_fw    = P_dem_simulink;
assignin('base','SOC_ref_fw',SOC_ref_fw);
assignin('base','H2_ref_fw',H2_ref_fw);
assignin('base','Pdem_fw',Pdem_fw);
assignin('base','SOC0',SOC0);
assignin('base','Tsim',Tsim);

fprintf('  N = %d, Ts = %.3g, Tsim = %.1f s\n', N, Ts, Tsim);
fprintf('  Pdem range   = [%+.3e, %+.3e] W\n', min(Pdem), max(Pdem));
fprintf('  SOCref range = [%.6f, %.6f]\n', min(SOCref), max(SOCref));
fprintf('  H2ref range  = [%.3f, %.3f] g\n', min(H2ref), max(H2ref));

%% ========================================================================
%  STEP 2: Build plant and synthesize SISO Hinf controller
% ========================================================================

fprintf('\n=== Building plant and sensor allocator Hinf controller ===\n');

s = tf('s');

num = [5.51 2.248e6 2.488e9 6.473e11];
den = [1044 1.239e10 2.034e13 8.21e15 3.67e16 0];

Gfc  = tf(num, den);        % Mh2 / Pfc
Gbat = -(3.858e-9)/s;       % SOC_delta / Pbat

% Full linear placeholder plant: inputs [Pdem; Pfc], outputs [SOC_delta; H2]
Gplant_tf = [Gbat, -Gbat;
             0*s,   Gfc];
Gplant_ct = minreal(ss(Gplant_tf));

% Allocation gains computed against open-loop plant from Pfc
w_alloc = 2*pi*f_alloc;
G_Pfc_to_SOC = -Gbat;
G_Pfc_to_H2  =  Gfc;

ell1 = squeeze(abs(freqresp(G_Pfc_to_SOC, w_alloc)));
ell2 = squeeze(abs(freqresp(G_Pfc_to_H2,  w_alloc)));

den_alloc = w1*ell1^2 + w2*ell2^2 + 1e-12;

fprintf('  ell1 = %.6g\n', ell1);
fprintf('  ell2 = %.6g\n', ell2);
fprintf('  den_alloc = %.6g\n', den_alloc);

% Nominal SISO plant from Pfc to allocated scalar error ec
%
% Sign analysis:
%   Ge_nom maps Pfc -> ec (through plant and error computation).
%   Ge_nom has NEGATIVE DC gain (increasing Pfc -> SOC rises -> e1 drops -> ec drops).
%   mixsyn designs K for standard negative feedback: S = 1/(1 + Ge_nom*K).
%   With negative Ge_nom, K comes out negative, so K*ec is negative when ec>0.
%   In Simulink: ec>0 (SOC too low) -> K output negative -> but we WANT
%   positive Pfc (to raise SOC by replacing battery discharge with FC power).
%
%   The linear sim's algebraic-loop solver handles this through the
%   Dp feedthrough terms, which effectively flip the sign internally.
%   In Simulink there is no algebraic loop (continuous plant, discrete K),
%   so we must correct the sign explicitly.
%
%   Solution: negate Ge_nom before synthesis (same as flip_sign=1).
%   This makes K positive, so ec>0 -> Pfc>0 -> SOC rises. Correct.

Ge_nom = -( w1*ell1*(-Gbat)/SOC_scale ...
          + w2*ell2* Gfc   /H2_scale ) / den_alloc;
Ge_nom = minreal(Ge_nom, 1e-10);

% Negate Ge_nom so that K comes out with positive gain
Ge_nom = -Ge_nom;

fprintf('  Ge_nom DC gain (after flip) = %+.4e\n', dcgain(Ge_nom));

% Hinf synthesis
Pfc_max_design = Pfc_max;
wb = 2*pi*wc_hz;

W1 = (s/M_w1 + wb) / (s + wb*A_w1);
W2 = tf(W2scale / Pfc_max_design);
W3 = [];

[K_hinf, ~, gamma_ach] = mixsyn(Ge_nom, W1, W2, W3);
K_ct = minreal(K_hinf, 1e-10);
K_ct = ss(K_ct);
K_d = c2d(K_ct, Ts, 'tustin');

[Ak, Bk, Ck, Dk] = ssdata(K_d);

% K is designed against -Ge_nom, so it has POSITIVE gain.
% ec > 0 -> K output > 0 -> Pfc increases -> SOC rises -> correct.
% No output negation needed.

fprintf('  gamma = %.4f\n', gamma_ach);
fprintf('  controller order = %d\n', size(Ak,1));
fprintf('  Dk = %.6e\n', Dk);

% Placeholder plant discretization
Gplant_d = c2d(Gplant_ct, Ts, 'tustin');
[Ap, Bp, Cp, Dp] = ssdata(Gplant_d);

% Assign variables used by Simulink blocks
assignin('base','Ak_hinf',Ak);
assignin('base','Bk_hinf',Bk);
assignin('base','Ck_hinf',Ck);
assignin('base','Dk_hinf',Dk);

assignin('base','Ap_placeholder',Ap);
assignin('base','Bp_placeholder',Bp);
assignin('base','Cp_placeholder',Cp);
assignin('base','Dp_placeholder',Dp);

assignin('base','SOC_scale',SOC_scale);
assignin('base','H2_scale',H2_scale);
assignin('base','ell1',ell1);
assignin('base','ell2',ell2);
assignin('base','w1_alloc',w1);
assignin('base','w2_alloc',w2);
assignin('base','den_alloc',den_alloc);
assignin('base','KsocI',KsocI);
assignin('base','leak',leak);
assignin('base','Pfc_min',Pfc_min);
assignin('base','Pfc_max',Pfc_max);
assignin('base','Pfc_soft_floor',Pfc_soft_floor);
assignin('base','Ts_controller',Ts);


%% PAPER FIGURE from Simulink run: SOC tracking / M_H2 tracking / errors
%% =========================================================================
% Signal map (from the Simulink model):
%   out.simout  -> SOC reference   (SOC_T)
%   out.simout1 -> M_H2 reference  (M_H2,T)
%   out.simout2 -> SOC output
%   out.simout3 -> M_H2 output

[t_socref, SOCref_sim] = local_getsig(out.simout);
[t_h2ref,  H2ref_sim ] = local_getsig(out.simout1);
[t_soc,    SOC_out   ] = local_getsig(out.simout2);
[t_h2,     Mh2_out   ] = local_getsig(out.simout3);

% Put everything on the SOC-output time base (handles variable-step solvers)
t = t_soc(:);
if numel(t_socref) ~= numel(t) || any(abs(t_socref(:) - t) > 1e-9)
    SOCref_sim = interp1(t_socref, SOCref_sim, t, 'linear', 'extrap');
end
if numel(t_h2ref) ~= numel(t) || any(abs(t_h2ref(:) - t) > 1e-9)
    H2ref_sim = interp1(t_h2ref, H2ref_sim, t, 'linear', 'extrap');
end
if numel(t_h2) ~= numel(t) || any(abs(t_h2(:) - t) > 1e-9)
    Mh2_out = interp1(t_h2, Mh2_out, t, 'linear', 'extrap');
end

SOC_out    = SOC_out(:);      Mh2_out   = Mh2_out(:);
SOCref_sim = SOCref_sim(:);   H2ref_sim = H2ref_sim(:);

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

% ---- helper (place at end of script, or save as a separate .m file) ------
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