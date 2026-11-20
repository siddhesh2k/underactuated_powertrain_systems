%% Tunable MPC for the FCHEV underactuated tracking problem
%
% Single-run MPC with one main knob: q1/q2 ratio.
% Edit the parameters in the "TUNE HERE" block, run the script, see the
% tracking result. Designed for fast iteration.
%
% Other knobs you might want to play with:
%   Np, Nc                 -- prediction / control horizon
%   r_du                   -- rate (smoothness) penalty on Pfc
%   Pfc_min / Pfc_max      -- saturation limits
%   Pfc_rate_limit         -- slew-rate limit (set to inf to disable)

clear; close all; clc;
s = tf('s');

%% =========================================================================
%% TUNE 
%% =========================================================================
q1_over_q2     = 1e8;       % SOC vs H2 priority (1e4 = H2 priority, 1e12 = SOC priority)
Np             = 30;         % prediction horizon (steps)
Nc             = 10;         % control horizon (steps)
r_du           = 1e-12;      % rate penalty on Pfc (smoothness)
Pfc_min        = 0;          % min fuel-cell power (W)
Pfc_max        = 106e3;      % max fuel-cell power (W)
Pfc_rate_limit = inf;        % slew-rate limit (W per step); inf = disabled

%% =========================================================================
%% Load data
%% =========================================================================
load('FCHEV_SDP_target_pdem_data_UDDS.mat')

N  = 1370;
Ts = 1;
t  = (0:N-1).' * Ts;

Pdem   = P_dem_simulink(2,:).';   Pdem   = Pdem(1:N);
SOCref = SOC_Target(2,:).';       SOCref = SOCref(1:N);
H2ref  = H2_cons_target(2,:).';   H2ref  = H2ref(1:N);
SOC0 = 0.6;

%% =========================================================================
%% Plant: inputs [Pdem; Pfc] -> outputs [SOC_delta; M_h2]
%% =========================================================================
num = [5.51 2.248e6 2.488e9 6.473e11];
den = [1044 1.239e10 2.034e13 8.21e15 3.67e16 0];
Gfc  = tf(num, den);
Gbat = -(3.858e-9)/s;

Gplant_tf = [ Gbat,  -Gbat ;
              0*s,    Gfc  ];
Gplant_ct = minreal(ss(Gplant_tf));
Gplant_d  = c2d(Gplant_ct, Ts, 'tustin');

[Ap, Bp, Cp, Dp] = ssdata(Gplant_d);
Bd = Bp(:,1);   Bu = Bp(:,2);
Dd = Dp(:,1);   Du = Dp(:,2);

nx = size(Ap,1);
ny = size(Cp,1);

%% =========================================================================
%% Build prediction matrices
%% =========================================================================
Phi      = zeros(Np*ny, nx);
Gamma_u  = zeros(Np*ny, Nc);
Gamma_d  = zeros(Np*ny, Np);

A_pow = eye(nx);
for j = 1:Np
    A_pow = Ap * A_pow;
    Phi((j-1)*ny+1 : j*ny, :) = Cp * A_pow;
end

% Gamma_u: indirect contribution Cp*A^k*Bu
for col = 1:Nc
    for j = col:Np
        idx_y = (j-1)*ny+1 : j*ny;
        if j == col
            Gamma_u(idx_y, col) = Cp * Bu;
        else
            Gamma_u(idx_y, col) = Cp * Ap^(j-col) * Bu;
        end
    end
end
% Hold beyond Nc
for j = Nc+1:Np
    idx_y = (j-1)*ny+1 : j*ny;
    add = zeros(ny, 1);
    for held = Nc:j
        if held == j
            add = add + Du;
        else
            add = add + Cp * Ap^(j-held) * Bu;
        end
    end
    Gamma_u(idx_y, Nc) = Gamma_u(idx_y, Nc) + add;
end
% Direct feedthrough
for j = 1:Nc
    idx_y = (j-1)*ny+1 : j*ny;
    Gamma_u(idx_y, j) = Gamma_u(idx_y, j) + Du;
end

% Gamma_d (constant Pdem assumption)
for col = 1:Np
    for j = col:Np
        idx_y = (j-1)*ny+1 : j*ny;
        if j == col
            Gamma_d(idx_y, col) = Cp * Bd;
        else
            Gamma_d(idx_y, col) = Cp * Ap^(j-col) * Bd;
        end
    end
end
for j = 1:Np
    idx_y = (j-1)*ny+1 : j*ny;
    Gamma_d(idx_y, j) = Gamma_d(idx_y, j) + Dd;
end
gamma_d_const = sum(Gamma_d, 2);

%% =========================================================================
%% Cost weights
%% =========================================================================
q1 = q1_over_q2;
q2 = 1;
Q_per_step = diag([q1, q2]);
Q_blk = kron(eye(Np), Q_per_step);

% Rate matrix (penalty on U[k] - U[k-1])
D_diff = eye(Nc);
for ii = 2:Nc
    D_diff(ii, ii-1) = -1;
end
R_du = r_du * (D_diff' * D_diff);

H_qp = Gamma_u' * Q_blk * Gamma_u + R_du;
H_qp = (H_qp + H_qp')/2;

%% =========================================================================
%% Slew-rate constraint matrix (if enabled)
%% =========================================================================
if isfinite(Pfc_rate_limit)
    % |U[k] - U[k-1]| <= rate_limit  for k = 1..Nc-1
    % Plus |U[1] - U_prev| <= rate_limit (handled in loop with current U_prev)
    A_rate = zeros(2*(Nc-1), Nc);
    b_rate = Pfc_rate_limit * ones(2*(Nc-1), 1);
    for ii = 1:Nc-1
        A_rate(2*ii-1, ii)   = -1;
        A_rate(2*ii-1, ii+1) =  1;
        A_rate(2*ii,   ii)   =  1;
        A_rate(2*ii,   ii+1) = -1;
    end
else
    A_rate = [];
    b_rate = [];
end



%% =========================================================================
%% Closed-loop simulation
%% =========================================================================
SOC_abs_log = zeros(N,1);
Mh2_log     = zeros(N,1);
Pfc_log     = zeros(N,1);
Pbat_log    = zeros(N,1);

xp = zeros(nx,1);
U_prev = 0;

lb = Pfc_min * ones(Nc,1);
ub = Pfc_max * ones(Nc,1);
opts_qp = optimoptions('quadprog','Display','off','MaxIterations',200);

fprintf('=== MPC Simulation ===\n');
fprintf('  q1/q2 = %.0e,  Np = %d,  Nc = %d,  r_du = %.0e\n', ...
    q1_over_q2, Np, Nc, r_du);
fprintf('  Pfc bounds: [%.0f, %.0f] W,  Pfc rate limit: %s W/step\n', ...
    Pfc_min, Pfc_max, num2str(Pfc_rate_limit));

tic;
for ii = 1:N
    % Reference vector over horizon
    Yref = zeros(Np*ny, 1);
    for j = 1:Np
        kk = min(ii + j - 1, N);
        Yref((j-1)*ny + 1) = SOCref(kk) - SOC0;
        Yref((j-1)*ny + 2) = H2ref(kk);
    end
    
    % ====== Constant P_dem ======
    Pdem_now = Pdem(ii);
    D_contrib = gamma_d_const * Pdem_now;


    % Online Pdem prediction: simple local trend predictor
    % ============================================================
    % 
    % Pdem_now = Pdem(ii);
    % 
    % if ii == 1
    %     dPdem = 0;
    % else
    %     dPdem = Pdem(ii) - Pdem(ii-1);
    % end
    % 
    % Pdem_pred = zeros(Np,1);
    % 
    % for j = 1:Np
    %     Pdem_pred(j) = Pdem_now + (j-1)*dPdem;
    % end
    % 
    % % Engineering clipping limits
    % % Adjust these if your vehicle power limits are different
    % Pdem_min_pred = -78e3;
    % Pdem_max_pred = 106e3 + 78e3;
    % 
    % Pdem_pred = min(max(Pdem_pred, Pdem_min_pred), Pdem_max_pred);
    % 
    % % Predicted disturbance contribution over horizon
    % D_contrib = Gamma_d * Pdem_pred;
    % ==========================================================

    f_qp = 2 * Gamma_u' * Q_blk * (Phi*xp + D_contrib - Yref);
    f_qp = f_qp - 2 * R_du(:,1) * U_prev;

    % Slew-rate constraint on first step (relative to U_prev)
    if isfinite(Pfc_rate_limit)
        A_first = zeros(2, Nc);
        A_first(1, 1) =  1;
        A_first(2, 1) = -1;
        b_first = [Pfc_rate_limit + U_prev; Pfc_rate_limit - U_prev];
        A_ineq = [A_rate; A_first];
        b_ineq = [b_rate; b_first];
    else
        A_ineq = [];
        b_ineq = [];
    end

    [U_opt, ~, exitflag] = quadprog(H_qp, f_qp, A_ineq, b_ineq, ...
        [], [], lb, ub, U_prev*ones(Nc,1), opts_qp);

    if exitflag ~= 1 && exitflag ~= 2
        warning('QP at step %d returned exitflag = %d', ii, exitflag);
        U_opt = U_prev * ones(Nc,1);
    end

    Pfc_k = min(max(U_opt(1), Pfc_min), Pfc_max);

    uplant = [Pdem(ii); Pfc_k];
    y = Cp*xp + Dp*uplant;
    SOC_abs_log(ii) = SOC0 + y(1);
    Mh2_log(ii)     = y(2);
    Pfc_log(ii)     = Pfc_k;
    Pbat_log(ii)    = Pdem(ii) - Pfc_k;

    xp = Ap*xp + Bp*uplant;
    U_prev = Pfc_k;
end
sim_time = toc;
fprintf('  Simulation time: %.2f s  (%.1f ms / step)\n', sim_time, 1000*sim_time/N);

%% =========================================================================
%% Metrics
%% =========================================================================
SOC_err = SOC_abs_log - SOCref;
H2_err  = Mh2_log     - H2ref;
sat_frac = mean(Pfc_log <= Pfc_min + eps | Pfc_log >= Pfc_max - eps);

fprintf('\n=== Results ===\n');
fprintf('  SOC RMS error           = %.4e\n', rms(SOC_err));
fprintf('  SOC max |err|           = %.4e\n', max(abs(SOC_err)));
fprintf('  SOC end-of-cycle error  = %+.4e\n', SOC_err(end));
fprintf('  H2 RMS error            = %.4e\n', rms(H2_err));
fprintf('  H2 max |err|            = %.4e\n', max(abs(H2_err)));
fprintf('  H2 end-of-cycle error   = %+.4e\n', H2_err(end));
fprintf('  Pfc range used          = [%.3e, %.3e]\n', min(Pfc_log), max(Pfc_log));
fprintf('  Pfc saturation fraction = %.3f\n', sat_frac);

% Slew rate check
dPfc = diff(Pfc_log);
fprintf('  Pfc max |slew|/step     = %.3e W\n', max(abs(dPfc)));

%% =========================================================================
%% Plots
%% =========================================================================
fig = figure('Name', sprintf('MPC tracking: q1/q2 = %.0e', q1_over_q2), ...
    'Position', [100 100 1200 800]);
tl = tiledlayout(fig, 4, 1, 'Padding','compact','TileSpacing','compact');

ax = nexttile;
plot(ax, t, SOC_abs_log, 'b', t, SOCref, 'r--', 'LineWidth', 1.2);
grid(ax,'on'); ylabel(ax,'SOC');
title(ax, sprintf('SOC tracking (RMS = %.3e, end err = %+.3e)', ...
    rms(SOC_err), SOC_err(end)));
legend(ax,'SOC','SOCref','Location','best');

ax = nexttile;
plot(ax, t, Mh2_log, 'b', t, H2ref, 'r--', 'LineWidth', 1.2);
grid(ax,'on'); ylabel(ax,'M_{H_2}');
title(ax, sprintf('H_2 tracking (RMS = %.3g, end err = %+.3g)', ...
    rms(H2_err), H2_err(end)));
legend(ax,'M_{H_2}','H_2ref','Location','best');

ax = nexttile;
plot(ax, t, Pdem, 'k', t, Pfc_log, 'b', t, Pbat_log, 'g', 'LineWidth', 1.1);
hold(ax,'on'); yline(ax, Pfc_max, 'b:'); yline(ax, Pfc_min, 'b:');
grid(ax,'on'); ylabel(ax,'Power (W)');
title(ax, sprintf('Power split (sat frac = %.3f)', sat_frac));
legend(ax,'P_{dem}','P_{fc}','P_{bat}','Location','best');

ax = nexttile;
plot(ax, t, SOC_err, 'b', 'LineWidth', 1.0); hold(ax,'on');
yyaxis(ax, 'right');
plot(ax, t, H2_err, 'r', 'LineWidth', 1.0);
ylabel(ax, 'H_2 error');
yyaxis(ax, 'left');
ylabel(ax, 'SOC error');
xlabel(ax, 'Time (s)'); grid(ax, 'on');
title(ax, 'Tracking errors');
legend(ax, 'SOC err (left)', 'H_2 err (right)', 'Location', 'best');

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
exportgraphics(figpaper,'MPC_cons.png','Resolution',300);