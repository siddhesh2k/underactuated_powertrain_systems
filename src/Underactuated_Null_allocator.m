%% End-to-end time-domain simulation:
%   H-infinity controller K + null-space prefilter Ga(1,1) + H prefilter at w0=0.04
%
% Architecture (closed loop, with disturbance Pdem on plant):
%
%      [SOC_T]   [Ga11 0]   [H11 0]
%      [H2_T ] - [  0  1] - [H21 1] -> r ->  e = r - y  -> K(s) -> Pfc -> Plant -> y
%                                                                          ^
%                                                                         Pdem
%
% Plant (disturbance + actuator):
%    [SOC_delta]   [Gbat]            [-Gbat]
%    [   M_h2  ] = [ 0  ] * Pdem  +  [ Gfc  ] * Pfc
%
% (note: SOC_abs = SOC0 + SOC_delta with SOC0 = 0.6)
%
% This script assumes the following are in workspace from your earlier scripts:
%   - Tpwr1  (closed-loop reference-to-output transfer matrix, 2x2)
%   - Ktfpwr1 (the H-infinity controller K(s))
%   - Gppwr2 (the 2x1 plant from Pfc -> [SOC_delta; M_h2]) -- note this is
%             the plant the controller was designed against
%   - Ga    (null-space allocator output, we use Ga(1,1))
%
% If they aren't, run your H-infinity synthesis script first.

clear; close all; clc;

%% =========================================================================
%% STEP 1: Run controller-design script if needed

run('Hinf4_null_allocator.m');     % <-- adjust path/name if needed

% Sanity check
needed_vars = {'Tpwr1','Ktfpwr1','Gppwr2','Ga'};
for k = 1:length(needed_vars)
    assert(evalin('base', sprintf('exist(''%s'',''var'')', needed_vars{k})) == 1, ...
        '%s not in workspace. Run your H-infinity design script first.', needed_vars{k});
end

% Pull from base workspace
fprintf('=== Workspace check ===\n');
fprintf('  Tpwr1   present: %d\n', exist('Tpwr1','var'));
fprintf('  Ktfpwr1 present: %d\n', exist('Ktfpwr1','var'));
fprintf('  Gppwr2  present: %d\n', exist('Gppwr2','var'));
fprintf('  Ga      present: %d\n', exist('Ga','var'));

%% =========================================================================
%% STEP 2: Load reference data
%% =========================================================================
load('FCHEV_SDP_target_pdem_data_UDDS.mat');
% load('FCHEV_DP_target_UDDS.mat');
% load('FCHEV_pdem_data_UDDS.mat');
SOC_Target_full     = SOC_Target;          % 2 x 1370
H2_cons_target_full = H2_cons_target;      % 2 x 1370
P_dem_full          = P_dem_simulink;      % 2 x 1370

% Use row 2
SOCref = SOC_Target_full(2,:).';
H2ref  = H2_cons_target_full(2,:).';
Pdem   = P_dem_full(2,:).';

N  = length(SOCref);
Ts = 1;
t  = (0:N-1).' * Ts;

SOC0 = 0.6;       % SOC offset (so SOC_abs = SOC0 + SOC_delta)
SOCref_delta = SOCref - SOC0;   % delta-form for tracking against plant output

fprintf('\n=== Reference data ===\n');
fprintf('  N = %d, Ts = %.3g, T = %.1f s\n', N, Ts, N*Ts);
fprintf('  SOCref     range: [%.4f, %.4f]   (delta range: [%+.4e, %+.4e])\n', ...
    min(SOCref), max(SOCref), min(SOCref_delta), max(SOCref_delta));
fprintf('  H2ref      range: [%.3g, %.3g]\n', min(H2ref), max(H2ref));
fprintf('  Pdem       range: [%+.3e, %+.3e]\n', min(Pdem), max(Pdem));

%% =========================================================================
%% STEP 3: Build the H prefilter at omega0 = 0.04 rad/s (Variant B)
%% =========================================================================
omega0 = 0.04;
Ty_w0 = evalfr(Tpwr1, 1j*omega0);
Ty_mag = abs(Ty_w0);
det_T_mag = Ty_mag(1,1)*Ty_mag(2,2) - Ty_mag(1,2)*Ty_mag(2,1);

H_mat = [ Ty_mag(1,1)*Ty_mag(2,2)/det_T_mag,         0;
         -Ty_mag(1,1)*Ty_mag(2,1)/det_T_mag,         1 ];

fprintf('\n=== H prefilter at omega0 = %.3g rad/s ===\n', omega0);
disp(H_mat);

%% =========================================================================
%% STEP 4: Assemble the plant + controller as state-space (single, big SS)
%% =========================================================================
% We need a discrete-time simulation. The plant we control was designed
% against Gppwr2 (2x1, only Pfc input, no disturbance).  But the *real* plant
% includes the Pdem disturbance.  We build a 2-input state-space:
%
%     x_p[k+1] = Ap*x_p[k] + Bp_full * [Pdem; Pfc]
%        y[k]  = Cp*x_p[k] + Dp_full * [Pdem; Pfc]
%
% with y = [SOC_delta; M_h2].
%
% From your physical model:
%   SOC_delta_dot ~ ( Pbat ) * (1/(3600*Em*Q))     and Pbat = Pdem - Pfc
%                 = (Pdem - Pfc)/(3600*Em*Q)
%   so dG/dPdem = +k_soc,  dG/dPfc = -k_soc  with k_soc = 1/(3600*Em*Q)
%
%   M_h2_dot ~ Pfc / (eta_fc * Q_lhv)
%   so dG/dPdem = 0,       dG/dPfc = +k_h2  with k_h2 = 1/(eta_fc*Q_lhv)
%
% The single-input plant Gppwr2 we used for the H-infinity design is
% [Gbat_to_SOC * (-1) ; Gfc].  i.e., it's the Pfc->y transfer (with the right
% sign baked in for the SOC channel because SOC = (Pdem - Pfc)*k -> -Pfc*k
% when we think of Pfc as the only input).
%
% For the closed-loop sim we use:
%   G_full(s) = [ +Gbat_SOC,  -Gbat_SOC ;
%                 0        ,  +Gfc       ]
% where Gbat_SOC = (3.858e-9)/s  (SOC channel sign convention: positive Pdem
% drains battery -> SOC_delta becomes negative; pick signs consistent with
% your H-infinity plant).
%
% Easiest: reuse Gppwr2 (Pfc-only plant) for the y = -Pfc-driven SOC channel,
% and add a parallel Pdem-driven channel for the disturbance.

s = tf('s');
Gbat = (3.858e-9)/s;      % positive: dSOC/dPbat with Pbat = Pdem - Pfc

% Verify Gppwr2 = [Gp1 ; Gp2] with Gp1 = -Gbat (so Pfc -> SOC_delta has correct sign)
% From your code: G1 = tf(num1, den1) with num1 = 3.858e-9, den1 = [1 0]
%   so G1 = +3.858e-9 / s = +Gbat
% Convention: y_SOC = G1 * Pfc means Pfc actuates SOC with POSITIVE gain.
% Physically Pbat = Pdem - Pfc, so dSOC/dPfc should be NEGATIVE.
% Your H-infinity controller was designed with Gppwr2 having +Gbat in row 1.
% We must keep that exact convention here for K to act correctly.
%
% Therefore the disturbance channel uses -Gbat (so Pdem - Pfc has the right
% relationship in the closed loop).  Let's verify by checking sum:
%   y_SOC = G1*Pfc + Gd_SOC*Pdem
%   Physically:  y_SOC ~ k * (Pdem - Pfc)
%   With G1 = +k/s (your convention), need Gd_SOC = -k/s ... wait
%   That would give y_SOC = +k/s * Pfc - k/s * Pdem = -k/s*(Pdem - Pfc)
%   which is the wrong sign physically.
%
% So your H-infinity controller was actually designed with a SIGN-FLIPPED
% SOC channel.  That's fine -- the controller compensates internally.
% For a consistent closed-loop sim that uses YOUR Tpwr1, we should use
% the SAME plant convention.  So:
Gd_full = [ -Gbat ;     % d SOC_delta / d Pdem  (sign flipped to match design)
             0    ];    % d M_h2      / d Pdem  = 0

% Full disturbance + actuator plant
Gplant_full_tf = [Gd_full, Gppwr2];   % 2 outputs, 2 inputs [Pdem, Pfc]
Gplant_full_ss = ss(minreal(Gplant_full_tf));
Gplant_d = c2d(Gplant_full_ss, Ts, 'tustin');

% Discrete controller
K_d = c2d(ss(Ktfpwr1), Ts, 'tustin');

[Ap, Bp, Cp, Dp] = ssdata(Gplant_d);
[Ak, Bk, Ck, Dk] = ssdata(K_d);
nxp = size(Ap, 1);
nxk = size(Ak, 1);

fprintf('\n=== Discrete-time models ===\n');
fprintf('  nxp (plant states)      = %d\n', nxp);
fprintf('  nxk (controller states) = %d\n', nxk);
fprintf('  ||Dp||_inf = %.3e   ||Dk||_inf = %.3e   (feedthrough check)\n', ...
    norm(Dp,inf), norm(Dk,inf));

%% =========================================================================
%% STEP 5: Discretize Ga(1,1) for the simulation
%% =========================================================================
Ga11_d = c2d(ss(Ga(1,1)), Ts, 'tustin');
[Aga, Bga, Cga, Dga] = ssdata(Ga11_d);
nxga = size(Aga, 1);

fprintf('\n=== Ga(1,1) ===\n');
fprintf('  Ga11 states = %d\n', nxga);
fprintf('  Ga11 DC gain = %.3e (%.1f dB)\n', dcgain(Ga(1,1)), ...
    20*log10(max(abs(dcgain(Ga(1,1))),eps)));

%% =========================================================================
%% STEP 6: Pfc saturation
%% =========================================================================
Pfc_min = 0;
Pfc_max = 106e3;
use_sat = true;

%% =========================================================================
%% STEP 7: Simulate the closed loop
%% =========================================================================
% Per-step computation:
%   1. r_soc_pre = Ga11_d * SOCref_delta[k]      (filter SOC ref through Ga11)
%   2. r_pre     = [r_soc_pre ; H2ref[k]]
%   3. r         = H_mat * r_pre                  (apply static H prefilter)
%   4. y         = Cp*xp + Dp*[Pdem; Pfc]         (read plant output)
%   5. e         = r - y                          (error)
%   6. Pfc       = Ck*xk + Dk*e                   (controller output)
%   7. saturate Pfc
%   8. advance states xk, xp, xga

% Storage
SOC_delta_log = zeros(N,1);
SOC_abs_log   = zeros(N,1);
Mh2_log       = zeros(N,1);
Pfc_log       = zeros(N,1);
Pbat_log      = zeros(N,1);
e_log         = zeros(N,2);
r_log         = zeros(N,2);
r_pre_log     = zeros(N,2);
r_soc_pre_log = zeros(N,1);

xp  = zeros(nxp, 1);
xk  = zeros(nxk, 1);
xga = zeros(nxga, 1);

% Algebraic-loop terms (Dk and Dp may both be nonzero after Tustin):
%   y = Cp*xp + Dp*[Pdem; Pfc]
%   Pfc = Ck*xk + Dk*e = Ck*xk + Dk*(r - y) = Ck*xk + Dk*r - Dk*Cp*xp - Dk*Dp*[Pdem; Pfc]
% Solve for Pfc explicitly:
%   (I + Dk*Dp(:,2)) * Pfc = Ck*xk + Dk*r - Dk*Cp*xp - Dk*Dp(:,1)*Pdem
% where Dp(:,1) is the Pdem column and Dp(:,2) is the Pfc column.

for ii = 1:N

    % --- 1. Ga(1,1) on SOC reference (delta) ---
    r_soc_pre = Cga*xga + Dga*SOCref_delta(ii);

    % --- 2. Pre-H reference vector ---
    r_pre = [r_soc_pre; H2ref(ii)];

    % --- 3. H prefilter (static 2x2) ---
    r = H_mat * r_pre;

    % --- 4-6. Algebraic-loop closure for Pfc ---
    %   y = Cp*xp + Dp(:,1)*Pdem + Dp(:,2)*Pfc
    %   e = r - y = r - Cp*xp - Dp(:,1)*Pdem - Dp(:,2)*Pfc
    %   Pfc = Ck*xk + Dk*e
    %   => (1 + Dk*Dp(:,2)) * Pfc = Ck*xk + Dk*(r - Cp*xp - Dp(:,1)*Pdem)
    rhs = Ck*xk + Dk*(r - Cp*xp - Dp(:,1)*Pdem(ii));
    M_loop = 1 + Dk*Dp(:,2);

    if abs(M_loop) < 1e-12
        warning('Algebraic-loop term near zero at step %d', ii);
        Pfc_k = rhs;     % fall back
    else
        Pfc_k = M_loop \ rhs;
    end

    % --- 7. Saturation ---
    if use_sat
        Pfc_k = min(max(Pfc_k, Pfc_min), Pfc_max);
    end

    % --- Read plant output with this Pfc ---
    uplant = [Pdem(ii); Pfc_k];
    y = Cp*xp + Dp*uplant;
    SOC_delta_k = y(1);
    Mh2_k       = y(2);
    SOC_abs_k   = SOC0 + SOC_delta_k;

    % Error (post-saturation, used for state update)
    e = r - y;

    % --- Log ---
    SOC_delta_log(ii) = SOC_delta_k;
    SOC_abs_log(ii)   = SOC_abs_k;
    Mh2_log(ii)       = Mh2_k;
    Pfc_log(ii)       = Pfc_k;
    Pbat_log(ii)      = Pdem(ii) - Pfc_k;
    e_log(ii,:)       = e.';
    r_log(ii,:)       = r.';
    r_pre_log(ii,:)   = r_pre.';
    r_soc_pre_log(ii) = r_soc_pre;

    % --- 8. Advance states ---
    xp  = Ap*xp + Bp*uplant;
    xk  = Ak*xk + Bk*e;
    xga = Aga*xga + Bga*SOCref_delta(ii);
end

%% =========================================================================
%% STEP 8: Metrics
%% =========================================================================
SOC_err = SOC_abs_log - SOCref;
H2_err  = Mh2_log     - H2ref;

fprintf('\n=========================================\n');
fprintf(' Time-domain results\n');
fprintf('=========================================\n');
fprintf('  SOC RMS error           = %.4e\n',     rms(SOC_err));
fprintf('  SOC max |err|           = %.4e\n',     max(abs(SOC_err)));
fprintf('  SOC end-of-cycle error  = %+.4e\n',    SOC_err(end));
fprintf('  H2 RMS error            = %.4e\n',     rms(H2_err));
fprintf('  H2 max |err|            = %.4e\n',     max(abs(H2_err)));
fprintf('  H2 end-of-cycle error   = %+.4e\n',    H2_err(end));
fprintf('  Pfc range used          = [%.3e, %.3e]\n', min(Pfc_log), max(Pfc_log));
fprintf('  Pfc saturation fraction = %.3f\n', ...
    mean(Pfc_log <= Pfc_min + eps | Pfc_log >= Pfc_max - eps));
fprintf('  Energy-bal RMS resid    = %.3e (Pdem - Pfc - Pbat, should be ~0)\n', ...
    rms(Pdem - Pfc_log - Pbat_log));

%% =========================================================================
%% STEP 9: Plots
%% =========================================================================
figure('Name','SOC and H2 tracking');
subplot(3,1,1);
plot(t, SOC_abs_log, 'b', t, SOCref, 'r--', 'LineWidth', 1.2); grid on;
legend('SOC','SOC ref','Location','best');
ylabel('SOC'); title('SOC tracking');

subplot(3,1,2);
plot(t, Mh2_log, 'b', t, H2ref, 'r--', 'LineWidth', 1.2); grid on;
legend('M_{H2}','H2 ref','Location','best');
ylabel('M_{H2}'); title('H_2 tracking');

subplot(3,1,3);
plot(t, SOC_err, 'b', t, H2_err, 'r', 'LineWidth', 1.1); grid on;
legend('SOC err','H_2 err','Location','best');
xlabel('Time (s)'); ylabel('Error');
title('Tracking errors');

figure('Name','Power split');
plot(t, Pdem, 'k', t, Pfc_log, 'b', t, Pbat_log, 'r', 'LineWidth', 1.2);
hold on;
yline(Pfc_max, 'b:'); yline(Pfc_min, 'b:');
grid on;
legend('P_{dem}','P_{fc}','P_{bat}','Location','best');
xlabel('Time (s)'); ylabel('Power (W)');
title('Power split (with P_{fc} saturation limits)');


%% =========================================================================
%% Closed-loop Bode plots: full reference-to-output map T_eff = Ty * H * diag(Ga11, 1)
M_tf = Tpwr1 * H_mat;     % 2x2 transfer matrix Ty*H

% Final closed-loop entries
T_eff_11 = minreal(M_tf(1,1) * Ga(1,1), 1e-2);
T_eff_12 = minreal(M_tf(1,2),           1e-2);
T_eff_21 = minreal(M_tf(2,1) * Ga(1,1), 1e-2);
T_eff_22 = minreal(M_tf(2,2),           1e-2);

% Attempt 1 (no H, just Ga11)
T_a1_11 = minreal(Tpwr1(1,1) * Ga(1,1), 1e-2);
T_a1_12 = Tpwr1(1,2);
T_a1_21 = minreal(Tpwr1(2,1) * Ga(1,1), 1e-2);
T_a1_22 = Tpwr1(2,2);

w_bode = logspace(-3, 3, 500);

% Helper to compute |H(jw)| in dB
magdb = @(sys, w) 20*log10(squeeze(abs(freqresp(sys, w))) + eps);

% Pre-compute all magnitudes
Ty11_mag    = magdb(Tpwr1(1,1), w_bode);
Ty12_mag    = magdb(Tpwr1(1,2), w_bode);
Ty21_mag    = magdb(Tpwr1(2,1), w_bode);
Ty22_mag    = magdb(Tpwr1(2,2), w_bode);

T_a1_11_mag = magdb(T_a1_11, w_bode);
T_a1_12_mag = magdb(T_a1_12, w_bode);
T_a1_21_mag = magdb(T_a1_21, w_bode);
T_a1_22_mag = magdb(T_a1_22, w_bode);

T_eff_11_mag = magdb(T_eff_11, w_bode);
T_eff_12_mag = magdb(T_eff_12, w_bode);
T_eff_21_mag = magdb(T_eff_21, w_bode);
T_eff_22_mag = magdb(T_eff_22, w_bode);

% Use tiledlayout for proper spacing
fig = figure('Name','Closed-loop reference-to-output Bode', ...
             'Position', [100 100 1200 800]);
tl = tiledlayout(fig, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

% --- (1,1): SOCref -> SOC ---
ax1 = nexttile;
semilogx(ax1, w_bode, Ty11_mag, 'LineWidth', 1.2); hold(ax1,'on');
semilogx(ax1, w_bode, T_a1_11_mag, 'LineWidth', 1.2);
semilogx(ax1, w_bode, T_eff_11_mag, 'LineWidth', 1.5);
xline(ax1, omega0, 'm--');
title(ax1, 'SOC_{ref} \rightarrow SOC');
ylabel(ax1, 'Magnitude (dB)');
grid(ax1, 'on');
legend(ax1, 'Ty(1,1) only','Ty(1,1)\cdotG_{a11}', ...
       'Final: (TyH)(1,1)\cdotG_{a11}','Location','best','FontSize',8);

% --- (1,2): H2ref -> SOC ---
ax2 = nexttile;
semilogx(ax2, w_bode, Ty12_mag, 'LineWidth', 1.2); hold(ax2,'on');
semilogx(ax2, w_bode, T_a1_12_mag, 'LineWidth', 1.2);
semilogx(ax2, w_bode, T_eff_12_mag, 'LineWidth', 1.5);
xline(ax2, omega0, 'm--');
title(ax2, 'H_{2,ref} \rightarrow SOC');
grid(ax2, 'on');
legend(ax2, 'Ty(1,2)','Attempt 1: Ty(1,2)','Final: (TyH)(1,2)', ...
       'Location','best','FontSize',8);

% --- (2,1): SOCref -> H2  [the problem path] ---
ax3 = nexttile;
semilogx(ax3, w_bode, Ty21_mag, 'LineWidth', 1.2); hold(ax3,'on');
semilogx(ax3, w_bode, T_a1_21_mag, 'LineWidth', 1.2);
semilogx(ax3, w_bode, T_eff_21_mag, 'LineWidth', 1.5);
xline(ax3, omega0, 'm--');
title(ax3, 'SOC_{ref} \rightarrow H_2   [problem path; want low]');
ylabel(ax3, 'Magnitude (dB)');
xlabel(ax3, 'Frequency (rad/s)');
grid(ax3, 'on');
legend(ax3, 'Ty(2,1) only','Ty(2,1)\cdotG_{a11} (Attempt 1)', ...
       'Final: (TyH)(2,1)\cdotG_{a11}','Location','best','FontSize',8);

% --- (2,2): H2ref -> H2 ---
ax4 = nexttile;
semilogx(ax4, w_bode, Ty22_mag, 'LineWidth', 1.2); hold(ax4,'on');
semilogx(ax4, w_bode, T_a1_22_mag, 'LineWidth', 1.2);
semilogx(ax4, w_bode, T_eff_22_mag, 'LineWidth', 1.5);
xline(ax4, omega0, 'm--');
title(ax4, 'H_{2,ref} \rightarrow H_2');
xlabel(ax4, 'Frequency (rad/s)');
grid(ax4, 'on');
legend(ax4, 'Ty(2,2)','Attempt 1: Ty(2,2)','Final: (TyH)(2,2)', ...
       'Location','best','FontSize',8);

title(tl, 'Closed-loop reference-to-output Bode plots');

%% Numerical comparison at DC and at omega0
fprintf('\n========================================================================\n');
fprintf(' Final closed-loop T_eff = Ty * H * diag(Ga11, 1)\n');
fprintf('========================================================================\n');
fprintf('  Path                    |  DC (dB)   | @ omega0=%.3g (dB)\n', omega0);
fprintf('  -------------------------------------------------------\n');
fprintf('  SOCref -> SOC           |  %+8.1f  |  %+8.1f\n', ...
    20*log10(max(abs(dcgain(T_eff_11)),eps)), ...
    20*log10(abs(evalfr(T_eff_11, 1j*omega0))));
fprintf('  H2ref  -> SOC           |  %+8.1f  |  %+8.1f\n', ...
    20*log10(max(abs(dcgain(T_eff_12)),eps)), ...
    20*log10(abs(evalfr(T_eff_12, 1j*omega0))));
fprintf('  SOCref -> H2 (PROBLEM)  |  %+8.1f  |  %+8.1f\n', ...
    20*log10(max(abs(dcgain(T_eff_21)),eps)), ...
    20*log10(abs(evalfr(T_eff_21, 1j*omega0))));
fprintf('  H2ref  -> H2            |  %+8.1f  |  %+8.1f\n', ...
    20*log10(max(abs(dcgain(T_eff_22)),eps)), ...
    20*log10(abs(evalfr(T_eff_22, 1j*omega0))));

fprintf('\n  For comparison, Attempt 1 (Ga11 only, no H):\n');
fprintf('  SOCref -> H2 (Attempt 1) |  %+8.1f  |  %+8.1f\n', ...
    20*log10(max(abs(dcgain(T_a1_21)),eps)), ...
    20*log10(abs(evalfr(T_a1_21, 1j*omega0))));

%% =========================================================================
%% PAPER FIGURE: SOC tracking / M_H2 tracking / errors
%% =========================================================================
SOC_out = SOC_abs_log;   Mh2_out = Mh2_log;
SOC_err_fig = SOC_out - SOCref;
H2_err_fig  = Mh2_out - H2ref;

figpaper = figure('Name','Paper figure: tracking','Color','w', ...
    'Position',[100 100 700 750]);
tlp = tiledlayout(figpaper,3,1,'Padding','compact','TileSpacing','compact');

axp = nexttile(tlp);
plot(axp, t, SOC_out, 'b-', 'LineWidth', 1.3); hold(axp,'on');
plot(axp, t, SOCref, 'r--', 'LineWidth', 1.3);
grid(axp,'on'); ylabel(axp,'$SOC$','Interpreter','latex');
legend(axp,{'$SOC$','$SOC_T$'},'Interpreter','latex','Location','best');
title(axp,'(a) SOC tracking');

axp = nexttile(tlp);
plot(axp, t, Mh2_out, 'b-', 'LineWidth', 1.3); hold(axp,'on');
plot(axp, t, H2ref, 'r--', 'LineWidth', 1.3);
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
exportgraphics(figpaper,'linear_HGa.png','Resolution',300);