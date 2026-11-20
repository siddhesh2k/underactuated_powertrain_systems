%% Sensor-allocation control with H-infinity SISO controller
% Uses the best-performing hyperparameters from the sweep:
%   wc_hz   = 0.15
%   M       = 1.2
%   A       = 0.01
%   W2scale = 0.05
%   flip    = 1
%   gamma_achieved = 0.8367

clear; close all; clc;

%% Load data
load('FCHEV_SDP_target_pdem_data_UDDS.mat')
s = tf('s');

%% Time base + reference extraction
N  = 1370;
Ts = 1;
t  = (0:N-1).';

Pdem   = P_dem_simulink(2,:).';
SOCref = SOC_Target(2,:).';
H2ref  = H2_cons_target(2,:).';

SOC0 = 0.6;

%% Sanity checks
assert(length(Pdem)   >= N, 'Pdem shorter than N');
assert(length(SOCref) >= N, 'SOCref shorter than N');
assert(length(H2ref)  >= N, 'H2ref shorter than N');

Pdem   = Pdem(1:N);
SOCref = SOCref(1:N);
H2ref  = H2ref(1:N);

fprintf('=== Data ranges ===\n');
fprintf('  Pdem:   [%+.3e, %+.3e] W\n', min(Pdem),   max(Pdem));
fprintf('  SOCref: [%+.4f,  %+.4f]\n',   min(SOCref), max(SOCref));
fprintf('  H2ref:  [%+.3e, %+.3e]\n',    min(H2ref),  max(H2ref));

if any(diff(H2ref) < 0)
    warning('H2ref is not monotonically non-decreasing -- check sign convention.');
end

%% Plant models
num = [5.51 2.248e6 2.488e9 6.473e11];
den = [1044 1.239e10 2.034e13 8.21e15 3.67e16 0];

Gfc  = tf(num,den);        % Mh2 / Pfc
Gbat = -(3.858e-9)/s;      % SOC_delta / Pbat

Gplant_tf = [ Gbat, -Gbat ;
              0*s,   Gfc  ];
Gplant_ct = minreal(ss(Gplant_tf));

%% =========================================================================
%% Engineering normalization
%% =========================================================================
SOC_scale = 0.01;
H2_scale  = 6.0596e+03;

fprintf('\n=== Engineering scales ===\n');
fprintf('  SOC_scale = %.6g\n', SOC_scale);
fprintf('  H2_scale  = %.6g\n', H2_scale);

%% =========================================================================
%% Allocation gains (computed against open-loop plant from Pfc)
%% =========================================================================
f_alloc = 0.3728;
w_alloc = 2*pi*f_alloc;

G_Pfc_to_SOC = -Gbat;
G_Pfc_to_H2  =  Gfc;

ell1 = squeeze(abs(freqresp(G_Pfc_to_SOC, w_alloc)));
ell2 = squeeze(abs(freqresp(G_Pfc_to_H2,  w_alloc)));

fprintf('\n=== Allocation gains at f_alloc = %.4g Hz ===\n', f_alloc);
fprintf('  ell1 = %.6g\n', ell1);
fprintf('  ell2 = %.6g\n', ell2);

%% Allocator weights
w1 = 1;
w2 = 1;
den_alloc = w1*ell1^2 + w2*ell2^2 + 1e-12;

%% =========================================================================
%% Build nominal SISO plant Ge_nom = Pfc -> ec  (with sign flip)
%% =========================================================================
Ge_nom = -( w1*ell1*(-Gbat)/SOC_scale ...
          + w2*ell2* Gfc   /H2_scale ) / den_alloc;
Ge_nom = minreal(Ge_nom, 1e-10);

% Apply best-found sign flip
flip_sign = 1;
if flip_sign
    Ge_nom = -Ge_nom;
end

fprintf('\n=== Nominal SISO plant Ge_nom ===\n');
fprintf('  DC gain (after flip) = %+.4e\n', dcgain(Ge_nom));
fprintf('  Order                = %d\n',     order(Ge_nom));

%% =========================================================================
%% H-infinity synthesis with BEST settings
%% =========================================================================
wc_hz   = 0.15;
M_w1    = 1.2;
A_w1    = 0.01;
W2scale = 0.05;

Pfc_max_design = 106e3;
wb = 2*pi*wc_hz;

W1 = (s/M_w1 + wb) / (s + wb*A_w1);
W2 = tf(W2scale / Pfc_max_design);
W3 = [];

[K_hinf, ~, gamma_ach] = mixsyn(Ge_nom, W1, W2, W3);
K_ct = minreal(K_hinf, 1e-10);
K_ct = ss(K_ct);

fprintf('\n=== H-inf synthesis ===\n');
fprintf('  wc_hz   = %.4g\n', wc_hz);
fprintf('  M       = %.4g\n', M_w1);
fprintf('  A       = %.4g\n', A_w1);
fprintf('  W2scale = %.4g\n', W2scale);
fprintf('  flip    = %d\n',   flip_sign);
fprintf('  gamma   = %.4f  (target from sweep: 0.8367)\n', gamma_ach);
fprintf('  controller order = %d\n', order(K_ct));

%% Discretize
Gplant = c2d(Gplant_ct, Ts, 'tustin');
K_d    = c2d(K_ct,      Ts, 'tustin');

[Ap,Bp,Cp,Dp] = ssdata(Gplant);
[Ak,Bk,Ck,Dk] = ssdata(K_d);

nxp = size(Ap,1);
nxk = size(Ak,1);

fprintf('\n=== Discrete feedthrough check ===\n');
fprintf('  ||Dp||_inf = %.3e   ||Dk|| = %.3e\n', norm(Dp,inf), abs(Dk));

%% =========================================================================
%% Pfc limits + soft floor
%% =========================================================================
Pfc_min = 0;
Pfc_max = 106e3;
Pfc_soft_floor = 2e3;
use_Pfc_sat = true;

%% =========================================================================
%% SOC integral bias
%% =========================================================================
KsocI = -2e6;
leak  = 0.002;
Isoc  = 0;

%% =========================================================================
%% MAIN SIMULATION
%% =========================================================================
fprintf('\n=== Running main simulation ===\n');
sim_out = run_allocator_sim_revised(Ap,Bp,Cp,Dp, Ak,Bk,Ck,Dk, ...
    Pdem, SOCref, H2ref, SOC0, ...
    ell1, ell2, w1, w2, SOC_scale, H2_scale, ...
    Pfc_min, Pfc_max, Pfc_soft_floor, use_Pfc_sat, ...
    KsocI, leak, Isoc, Ts);

[metrics_nom] = report_metrics(sim_out, SOCref, H2ref, Pdem, Pfc_min, Pfc_max, 'H-inf nominal');

%% =========================================================================
%% Narrow f_alloc SWEEP
%% =========================================================================
fprintf('\n=== Sweep: f_alloc (narrow region) ===\n');
f_sweep = logspace(log10(0.12), log10(0.6), 9);
sweep_f_metrics = struct();

fprintf('   f_alloc (Hz)  |  ell1/ell2    |  SOC RMS      |  H2 RMS       |  SOC end err  |  H2 end err\n');
fprintf('   ----------------------------------------------------------------------------------------\n');

for k = 1:length(f_sweep)
    fk = f_sweep(k);
    wk = 2*pi*fk;

    e1k = squeeze(abs(freqresp(G_Pfc_to_SOC, wk)));
    e2k = squeeze(abs(freqresp(G_Pfc_to_H2,  wk)));

    out = run_allocator_sim_revised(Ap,Bp,Cp,Dp, Ak,Bk,Ck,Dk, ...
        Pdem, SOCref, H2ref, SOC0, ...
        e1k, e2k, w1, w2, SOC_scale, H2_scale, ...
        Pfc_min, Pfc_max, Pfc_soft_floor, use_Pfc_sat, ...
        KsocI, leak, 0, Ts);

    soc_rms = rms(out.SOC_abs - SOCref);
    h2_rms  = rms(out.Mh2     - H2ref);
    soc_end = out.SOC_abs(end) - SOCref(end);
    h2_end  = out.Mh2(end)     - H2ref(end);

    sweep_f_metrics(k).f_alloc = fk;
    sweep_f_metrics(k).ell1    = e1k;
    sweep_f_metrics(k).ell2    = e2k;
    sweep_f_metrics(k).soc_rms = soc_rms;
    sweep_f_metrics(k).h2_rms  = h2_rms;
    sweep_f_metrics(k).soc_end = soc_end;
    sweep_f_metrics(k).h2_end  = h2_end;

    fprintf('   %10.4g  |  %8.2e   |  %10.4g   |  %10.4g   |  %+10.4g   |  %+10.4g\n', ...
        fk, e1k/max(e2k,eps), soc_rms, h2_rms, soc_end, h2_end);
end

figure('Name','f_alloc sweep (H-inf)');
subplot(2,1,1);
loglog(f_sweep, [sweep_f_metrics.soc_rms], '-o', ...
       f_sweep, [sweep_f_metrics.h2_rms],  '-s', 'LineWidth', 1.2);
xlabel('f_{alloc} (Hz)'); ylabel('RMS tracking error'); grid on;
legend('SOC RMS','H2 RMS','Location','best');
title('Tracking error vs allocation frequency');

subplot(2,1,2);
semilogx(f_sweep, [sweep_f_metrics.soc_end], '-o', ...
         f_sweep, [sweep_f_metrics.h2_end],  '-s', 'LineWidth', 1.2);
xlabel('f_{alloc} (Hz)'); ylabel('End error (signed)'); grid on;
legend('SOC end err','H2 end err','Location','best');
title('End-of-cycle error vs allocation frequency');

%% =========================================================================
%% Narrow weight-ratio SWEEP
%% =========================================================================
fprintf('\n=== Sweep: weight ratio r = w1/w2 (narrow useful region) ===\n');
r_sweep = logspace(-1, 1.5, 9);
sweep_w_metrics = struct();

fprintf('    w1/w2    |  SOC RMS      |  H2 RMS       |  SOC end err  |  H2 end err   |  Pfc sat frac\n');
fprintf('   ----------------------------------------------------------------------------------------\n');

for k = 1:length(r_sweep)
    r = r_sweep(k);
    w1k = 2*r/(1+r);
    w2k = 2  /(1+r);

    out = run_allocator_sim_revised(Ap,Bp,Cp,Dp, Ak,Bk,Ck,Dk, ...
        Pdem, SOCref, H2ref, SOC0, ...
        ell1, ell2, w1k, w2k, SOC_scale, H2_scale, ...
        Pfc_min, Pfc_max, Pfc_soft_floor, use_Pfc_sat, ...
        KsocI, leak, 0, Ts);

    soc_rms = rms(out.SOC_abs - SOCref);
    h2_rms  = rms(out.Mh2     - H2ref);
    soc_end = out.SOC_abs(end) - SOCref(end);
    h2_end  = out.Mh2(end)     - H2ref(end);
    sat_frac = mean(out.Pfc <= Pfc_min + eps | out.Pfc >= Pfc_max - eps);

    sweep_w_metrics(k).ratio    = r;
    sweep_w_metrics(k).w1       = w1k;
    sweep_w_metrics(k).w2       = w2k;
    sweep_w_metrics(k).soc_rms  = soc_rms;
    sweep_w_metrics(k).h2_rms   = h2_rms;
    sweep_w_metrics(k).soc_end  = soc_end;
    sweep_w_metrics(k).h2_end   = h2_end;
    sweep_w_metrics(k).sat_frac = sat_frac;

    fprintf('   %8.3g  |  %10.4g   |  %10.4g   |  %+10.4g   |  %+10.4g   |  %7.3f\n', ...
        r, soc_rms, h2_rms, soc_end, h2_end, sat_frac);
end

figure('Name','w1/w2 sweep (H-inf)');
subplot(2,1,1);
loglog(r_sweep, [sweep_w_metrics.soc_rms], '-o', ...
       r_sweep, [sweep_w_metrics.h2_rms],  '-s', 'LineWidth', 1.2);
xlabel('w_1/w_2'); ylabel('RMS tracking error'); grid on;
legend('SOC RMS','H2 RMS','Location','best');
title('Tracking error vs weight ratio');

subplot(2,1,2);
loglog([sweep_w_metrics.soc_rms], [sweep_w_metrics.h2_rms], '-o', 'LineWidth', 1.2);
hold on;
[~, idx_nom] = min(abs(r_sweep - 1));
plot(sweep_w_metrics(idx_nom).soc_rms, sweep_w_metrics(idx_nom).h2_rms, ...
     'rp', 'MarkerFaceColor','r','MarkerSize',12);
xlabel('SOC RMS'); ylabel('H2 RMS'); grid on;
title('Pareto front: SOC vs H2 tracking (red marker = w1=w2)');

%% =========================================================================
%% Main plots (nominal sim)
%% =========================================================================
figure('Name','Tracking (H-inf nominal)');
subplot(4,1,1);
plot(t, sim_out.SOC_abs, t, SOCref, 'LineWidth', 1.1); grid on;
legend('SOC','SOCref','Location','best');
title('SOC tracking');

subplot(4,1,2);
plot(t, sim_out.Mh2, t, H2ref, 'LineWidth', 1.1); grid on;
legend('Mh2','H2ref','Location','best');
title('H2 tracking');

subplot(4,1,3);
plot(t, sim_out.ec, 'LineWidth', 1.1); grid on;
title('Allocated error command e_c');

subplot(4,1,4);
plot(t, sim_out.Isoc_log, 'LineWidth', 1.1); grid on;
title('SOC integral state');
xlabel('Time (s)');

figure('Name','Power split (H-inf nominal)');
plot(t, Pdem, 'k', t, sim_out.Pfc, 'b', t, sim_out.Pbat, 'r', 'LineWidth', 1.1);
grid on;
legend('P_{dem}','P_{fc}','P_{bat}','Location','best');
title('Power split');
xlabel('Time (s)'); ylabel('Power (W)');

figure('Name','Errors (H-inf)');
subplot(2,1,1); plot(t, sim_out.e1, 'LineWidth', 1.1); grid on;
title('SOC error e_1'); ylabel('SOC err');
subplot(2,1,2); plot(t, sim_out.e2, 'LineWidth', 1.1); grid on;
title('H2 error e_2'); xlabel('Time (s)'); ylabel('H2 err');

%% =========================================================================
%% PAPER FIGURE: SOC tracking / M_H2 tracking / errors
%% =========================================================================
SOC_out = sim_out.SOC_abs;   Mh2_out = sim_out.Mh2;
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
exportgraphics(figpaper,'linear_SensorAlloc.png','Resolution',300);

%% =========================================================================
%% Helper: revised simulation with exact algebraic-loop solve + SOC bias
%% =========================================================================
function out = run_allocator_sim_revised(Ap,Bp,Cp,Dp, Ak,Bk,Ck,Dk, ...
                                         Pdem, SOCref, H2ref, SOC0, ...
                                         ell1, ell2, w1, w2, SOC_scale, H2_scale, ...
                                         Pfc_min, Pfc_max, Pfc_soft_floor, use_Pfc_sat, ...
                                         KsocI, leak, Isoc0, Ts)

    N = length(Pdem);
    nxp = size(Ap,1);
    nxk = size(Ak,1);

    xp = zeros(nxp,1);
    xk = zeros(nxk,1);
    Isoc = Isoc0;

    SOC_delta = zeros(N,1);
    SOC_abs   = zeros(N,1);
    Mh2       = zeros(N,1);
    Pfc       = zeros(N,1);
    Pbat      = zeros(N,1);
    e1        = zeros(N,1);
    e2        = zeros(N,1);
    e1n       = zeros(N,1);
    e2n       = zeros(N,1);
    ec        = zeros(N,1);
    Isoc_log  = zeros(N,1);

    den_alloc = w1*ell1^2 + w2*ell2^2 + 1e-12;
    a = Dp(1,2);
    b = Dp(2,2);

    for ii = 1:N
        s0 = SOC0 + Cp(1,:)*xp + Dp(1,1)*Pdem(ii);
        m0 =        Cp(2,:)*xp + Dp(2,1)*Pdem(ii);

        alpha = ( w1*ell1*(SOCref(ii) - s0)/SOC_scale ...
                + w2*ell2*(H2ref(ii)  - m0)/H2_scale ) / den_alloc;
        beta  = ( w1*ell1*a/SOC_scale ...
                + w2*ell2*b/H2_scale ) / den_alloc;

        denom_loop = 1 + Dk*beta;
        if abs(denom_loop) < 1e-12
            warning('Algebraic-loop denominator near zero at step %d', ii);
            denom_loop = sign(denom_loop)*1e-12 + (denom_loop==0)*1e-12;
        end

        Pfc_k = (Ck*xk + Dk*alpha + KsocI*Isoc) / denom_loop;
        Pfc_k = max(Pfc_k, Pfc_soft_floor);
        if use_Pfc_sat
            Pfc_k = min(max(Pfc_k, Pfc_min), Pfc_max);
        end

        uplant = [Pdem(ii); Pfc_k];
        yplant = Cp*xp + Dp*uplant;
        SOCd_k   = yplant(1);
        Mh2_k    = yplant(2);
        SOCabs_k = SOC0 + SOCd_k;

        e1_k  = SOCref(ii) - SOCabs_k;
        e2_k  = H2ref(ii)  - Mh2_k;
        e1n_k = e1_k / SOC_scale;
        e2n_k = e2_k / H2_scale;
        ec_k  = (w1*ell1*e1n_k + w2*ell2*e2n_k) / den_alloc;

        Isoc = (1-leak)*Isoc + Ts*e1_k;

        SOC_delta(ii) = SOCd_k;
        SOC_abs(ii)   = SOCabs_k;
        Mh2(ii)       = Mh2_k;
        Pfc(ii)       = Pfc_k;
        Pbat(ii)      = Pdem(ii) - Pfc_k;
        e1(ii)        = e1_k;
        e2(ii)        = e2_k;
        e1n(ii)       = e1n_k;
        e2n(ii)       = e2n_k;
        ec(ii)        = ec_k;
        Isoc_log(ii)  = Isoc;

        xk = Ak*xk + Bk*ec_k;
        xp = Ap*xp + Bp*[Pdem(ii); Pfc_k];
    end

    out.SOC_delta = SOC_delta;
    out.SOC_abs   = SOC_abs;
    out.Mh2       = Mh2;
    out.Pfc       = Pfc;
    out.Pbat      = Pbat;
    out.e1        = e1;
    out.e2        = e2;
    out.e1n       = e1n;
    out.e2n       = e2n;
    out.ec        = ec;
    out.Isoc_log  = Isoc_log;
end

%% =========================================================================
%% Helper: metrics report
%% =========================================================================
function m = report_metrics(out, SOCref, H2ref, Pdem, Pfc_min, Pfc_max, label)
    fprintf('\n=== Metrics [%s] ===\n', label);

    m.ec_rms        = rms(out.ec);
    m.SOC_end_err   = out.SOC_abs(end) - SOCref(end);
    m.H2_end_err    = out.Mh2(end)     - H2ref(end);
    m.SOC_rms       = rms(out.SOC_abs - SOCref);
    m.H2_rms        = rms(out.Mh2     - H2ref);
    m.SOC_max_dev   = max(abs(out.SOC_abs - SOCref));
    m.H2_max_dev    = max(abs(out.Mh2     - H2ref));
    m.Pfc_sat_frac  = mean(out.Pfc <= Pfc_min + eps | out.Pfc >= Pfc_max - eps);
    m.Pfc_min_used  = min(out.Pfc);
    m.Pfc_max_used  = max(out.Pfc);
    m.energy_balance_rms = rms(Pdem - out.Pfc - out.Pbat);
    m.H2_violations = sum(diff(out.Mh2) < 0);

    fprintf('  ec_rms            = %.6g\n', m.ec_rms);
    fprintf('  SOC_end_err       = %+.6g\n', m.SOC_end_err);
    fprintf('  H2_end_err        = %+.6g\n', m.H2_end_err);
    fprintf('  SOC_rms           = %.6g\n', m.SOC_rms);
    fprintf('  H2_rms            = %.6g\n', m.H2_rms);
    fprintf('  SOC_max_dev       = %.6g\n', m.SOC_max_dev);
    fprintf('  H2_max_dev        = %.6g\n', m.H2_max_dev);
    fprintf('  Pfc_sat_frac      = %.3f\n',  m.Pfc_sat_frac);
    fprintf('  Pfc range used    = [%.3e, %.3e]\n', m.Pfc_min_used, m.Pfc_max_used);
    fprintf('  Energy-bal RMS    = %.3e   (should be ~0)\n', m.energy_balance_rms);
    fprintf('  H2 non-monotonic  = %d steps\n', m.H2_violations);
end