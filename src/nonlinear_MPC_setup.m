%% setup_MPC_workspace_FCHEV.m


clear; clc; close all;

%% ========================================================================
%  1. File names
% ========================================================================

param_file = 'PTsystem_Openloop_simulink_parameter_v03.m';
data_file  = 'FCHEV_SDP_target_pdem_data_UDDS.mat';

%% ========================================================================
%  2. Run nonlinear plant parameter file FIRST
% ========================================================================

fprintf('\n=== Running nonlinear plant parameter file ===\n');

if exist(param_file, 'file') ~= 2
    error('Could not find parameter file: %s', param_file);
end

run(param_file);

%% ========================================================================
%  3. Redefine everything after parameter file
% ========================================================================
% The parameter file may clear the workspace, so define user settings here.

param_file = 'PTsystem_Openloop_simulink_parameter_v03.m'; 
data_file  = 'FCHEV_SDP_target_pdem_data_UDDS.mat';

Ts   = 1;
SOC0 = 0.6;

% MPC tuning
q1_over_q2 = 1e8;     % SOC priority. Try 1e8, 1e10, 1e12.
Np = 30;
Nc = 10;

% Fuel-cell power constraints
Pfc_min = 0;
Pfc_max = 106e3;

% Optional fuel-cell slew-rate limits. Use Inf to disable.
Pfc_rate_min = -Inf;   % e.g., -5e3
Pfc_rate_max =  Inf;   % e.g.,  5e3

% Move suppression weight
mv_rate_weight = 0.1;

%% ========================================================================
%  4. Load target and demand data
% ========================================================================

fprintf('\n=== Loading target and demand data ===\n');

if exist(data_file, 'file') ~= 2
    error('Could not find data file: %s', data_file);
end

D = load(data_file);

required_vars = {'SOC_Target','H2_cons_target','P_dem_simulink'};
for k = 1:numel(required_vars)
    if ~isfield(D, required_vars{k})
        error('Missing variable %s in %s.', required_vars{k}, data_file);
    end
end

% Use row 2, matching your MPC script convention.
SOC_ref = D.SOC_Target(2,:).';
H2_ref  = D.H2_cons_target(2,:).';
P_dem   = D.P_dem_simulink(2,:).';

N = min([numel(SOC_ref), numel(H2_ref), numel(P_dem)]);
SOC_ref = SOC_ref(1:N);
H2_ref  = H2_ref(1:N);
P_dem   = P_dem(1:N);

t = (0:N-1).' * Ts;
Tsim = t(end);

% MPC tracks SOC_delta, not absolute SOC.
SOC_ref_delta = SOC_ref - SOC0;

% Simulink From Workspace signals must be [time, value].
SOC_ref_fw       = [t SOC_ref];
SOC_ref_delta_fw = [t SOC_ref_delta];
H2_ref_fw        = [t H2_ref];
Pdem_fw          = [t P_dem];

assignin('base','SOC_ref',SOC_ref);
assignin('base','SOC_ref_delta',SOC_ref_delta);
assignin('base','H2_ref',H2_ref);
assignin('base','P_dem',P_dem);
assignin('base','t',t);
assignin('base','Tsim',Tsim);

assignin('base','SOC_ref_fw',SOC_ref_fw);
assignin('base','SOC_ref_delta_fw',SOC_ref_delta_fw);
assignin('base','H2_ref_fw',H2_ref_fw);
assignin('base','Pdem_fw',Pdem_fw);

fprintf('  N = %d samples\n', N);
fprintf('  Ts = %.3g s\n', Ts);
fprintf('  Tsim = %.1f s\n', Tsim);
fprintf('  SOC_ref range       = [%.6f, %.6f]\n', min(SOC_ref), max(SOC_ref));
fprintf('  SOC_ref_delta range = [%+.4e, %+.4e]\n', min(SOC_ref_delta), max(SOC_ref_delta));
fprintf('  H2_ref range        = [%.3f, %.3f] g\n', min(H2_ref), max(H2_ref));
fprintf('  P_dem range         = [%.1f, %.1f] W\n', min(P_dem), max(P_dem));

%% ========================================================================
%  5. Build MPC prediction plant
% ========================================================================
% Inputs:
%   input 1 = P_dem measured disturbance
%   input 2 = P_fc manipulated variable
%
% Outputs:
%   output 1 = SOC_delta
%   output 2 = H2
%
% Important for Simulink MPC block:
%   Direct feedthrough from MV to output is not allowed.
%   So use ZOH discretization and force D(:,MV) = 0.

fprintf('\n=== Building MPC prediction plant ===\n');

s = tf('s');

Gfc_num = [5.51 2.248e6 2.488e9 6.473e11];
Gfc_den = [1044 1.239e10 2.034e13 8.21e15 3.67e16 0];
Gfc = tf(Gfc_num, Gfc_den);

% SOC sign convention:
%   positive P_dem lowers SOC
%   positive P_fc raises SOC
Gbat = -(3.858e-9)/s;

Gplant_tf = [Gbat, -Gbat;
             0*s,   Gfc];

Gplant_ct = minreal(ss(Gplant_tf));
Gplant_d  = c2d(Gplant_ct, Ts, 'zoh');

[Ampc, Bmpc, Cmpc, Dmpc] = ssdata(Gplant_d);

MD_index = 1;      % P_dem
MV_index = 2;      % P_fc
MO_index = [1 2];  % [SOC_delta; H2]

% Remove direct MV feedthrough for MPC block compatibility.
Dmpc(:,MV_index) = 0;

Gplant_d = ss(Ampc, Bmpc, Cmpc, Dmpc, Ts);
Gplant_d.InputName  = {'Pdem','Pfc'};
Gplant_d.OutputName = {'SOC_delta','H2'};

fprintf('  Direct MV feedthrough norm = %.3e\n', norm(Dmpc(:,MV_index), inf));

plant_mpc = setmpcsignals(Gplant_d, ...
    'MD', MD_index, ...
    'MV', MV_index, ...
    'MO', MO_index);

assignin('base','Gplant_d',Gplant_d);
assignin('base','plant_mpc',plant_mpc);

%% ========================================================================
%  6. Create MPC controller object
% ========================================================================

fprintf('\n=== Creating MPC controller object ===\n');

mpcobj = mpc(plant_mpc, Ts);

mpcobj.PredictionHorizon = Np;
mpcobj.ControlHorizon = Nc;

% MV constraints: P_fc
mpcobj.MV.Min = Pfc_min;
mpcobj.MV.Max = Pfc_max;

if isfinite(Pfc_rate_min)
    mpcobj.MV.RateMin = Pfc_rate_min;
end
if isfinite(Pfc_rate_max)
    mpcobj.MV.RateMax = Pfc_rate_max;
end

% Weights
mpcobj.Weights.OutputVariables = [q1_over_q2, 1];
mpcobj.Weights.ManipulatedVariables = 0;
mpcobj.Weights.ManipulatedVariablesRate = mv_rate_weight;

% Nominal values for deviation-based model
mpcobj.Model.Nominal.U = [0; 0];
mpcobj.Model.Nominal.Y = [0; 0];
mpcobj.Model.Nominal.X = zeros(size(Ampc,1),1);

assignin('base','mpcobj',mpcobj);
assignin('base','Ts',Ts);
assignin('base','SOC0',SOC0);
assignin('base','Pfc_min',Pfc_min);
assignin('base','Pfc_max',Pfc_max);
assignin('base','Np',Np);
assignin('base','Nc',Nc);
assignin('base','q1_over_q2',q1_over_q2);

fprintf('  MPC object created as variable: mpcobj\n');
fprintf('  Prediction horizon = %d\n', Np);
fprintf('  Control horizon    = %d\n', Nc);
fprintf('  Output weights     = [%.3e, %.3e]\n', q1_over_q2, 1);
fprintf('  MV limits          = [%.1f, %.1f] W\n', Pfc_min, Pfc_max);

%% ========================================================================
%  7. Summary for Simulink wiring
% ========================================================================

fprintf('\n=== Simulink variables ready ===\n');
fprintf('Use these From Workspace variables:\n');
fprintf('  SOC reference absolute: SOC_ref_fw\n');
fprintf('  SOC reference delta:    SOC_ref_delta_fw\n');
fprintf('  H2 reference:           H2_ref_fw\n');
fprintf('  Power demand:           Pdem_fw\n');
fprintf('\nMPC block should use controller object: mpcobj\n');
fprintf('\nCorrect MPC block wiring:\n');
fprintf('  mo  = [SOC_actual - SOC0; H2_actual]\n');
fprintf('  ref = [SOC_ref_delta;     H2_ref]\n');
fprintf('  md  = P_dem\n');
fprintf('  mv  = P_fc command -> Saturation -> Divide6 -> PEMFC\n');
fprintf('\nSet model StopTime to Tsim or %.1f.\n', Tsim);

fprintf('\nSetup complete. You can now run your MPC Simulink model.\n');


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
exportgraphics(figpaper,'nonlinear_MPC.png','Resolution',300);

%---- helper (place at end of script, or save as a separate .m file) ------
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