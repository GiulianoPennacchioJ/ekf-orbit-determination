% =========================================================================
% MASTER SCRIPT: GENERATE ALL EKF DIAGNOSTIC FIGURES (FULL 6 FIGURES)
% =========================================================================

close all; clc;

% --- 1. SETUP OUTPUT DIRECTORY ---
fig_dir = fullfile(pwd, 'figures');
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

% --- 2. LOAD SIMULATION DATA FROM WORKSPACE ---
if ~exist('out_x_true', 'var') || ~exist('out_x_hat', 'var') || ~exist('out_P_cov', 'var')
    error('Simulation outputs missing. Please run Simulink model first.');
end

x_true  = out_x_true.Data;
x_hat   = out_x_hat.Data;
y_innov = out_y_innov.Data;
S_innov = out_S_innov.Data;
t_vec   = out_x_true.Time;
P_data  = out_P_cov.Data;

N_samples  = length(t_vec);
comp_names = {'X', 'Y', 'Z'};

err_pos_xyz  = (x_hat(1:3, :) - x_true(1:3, :)) * 1000; % [m]
err_vel_xyz  = (x_hat(4:6, :) - x_true(4:6, :)) * 1000; % [mm/s]
err_pos_norm = sqrt(sum(err_pos_xyz.^2, 1));           % [m]

sigma3_p = [3 * sqrt(squeeze(P_data(1,1,:)))' * 1000; ...
            3 * sqrt(squeeze(P_data(2,2,:)))' * 1000; ...
            3 * sqrt(squeeze(P_data(3,3,:)))' * 1000]; % [m]
sigma3_v = [3 * sqrt(squeeze(P_data(4,4,:)))' * 1000; ...
            3 * sqrt(squeeze(P_data(5,5,:)))' * 1000; ...
            3 * sqrt(squeeze(P_data(6,6,:)))' * 1000]; % [mm/s]

% -------------------------------------------------------------------------
% FIGURE 1: TRANSITORIO INIZIALE (0 - 5 s ZOOM)
% -------------------------------------------------------------------------
fig1 = figure('Name', 'EKF Initial Convergence', 'Color', 'w', 'Position', [100, 100, 900, 500]);
idx_5s = t_vec <= 5;
for i = 1:3
    subplot(3, 1, i);
    plot(t_vec(idx_5s), err_pos_xyz(i, idx_5s), 'b-o', 'LineWidth', 1.2, 'MarkerFaceColor', 'b', 'MarkerSize', 3);
    grid on; ylabel(sprintf('Error %s [m]', comp_names{i}));
    if i == 1, title('EKF Initial Transient Step Response (0 - 5s Zoom)'); end
    if i == 3, xlabel('Time [s]'); end
end
set(findall(fig1, 'Type', 'axes'), 'Toolbar', []);
exportgraphics(fig1, fullfile(fig_dir, '01_initial_convergence.png'), 'Resolution', 300, 'ContentType', 'image');

% -------------------------------------------------------------------------
% FIGURE 2: POSIZIONE A REGIME VS 3-SIGMA
% -------------------------------------------------------------------------
fig2 = figure('Name', 'Position 3-Sigma Consistency', 'Color', 'w', 'Position', [100, 100, 900, 600]);
idx_ss = t_vec >= 100;
t_ss   = t_vec(idx_ss);
for i = 1:3
    subplot(3, 1, i);
    sig_i = sigma3_p(i, idx_ss);
    fill([t_ss; flipud(t_ss)], [sig_i(:); flipud(-sig_i(:))], [1 0.85 0.85], ...
         'EdgeColor', [0.8 0 0], 'LineStyle', '--', 'LineWidth', 1.0); hold on;
    plot(t_ss, err_pos_xyz(i, idx_ss), 'b', 'LineWidth', 0.6);
    grid on; ylabel(sprintf('Error %s [m]', comp_names{i})); ylim([-20, 20]); xlim([100, max(t_vec)]);
    if i == 1
        title('Steady-State Position Estimation Errors vs 3\sigma Envelope');
        legend('3\sigma Uncertainty Bound', 'Position Error', 'Location', 'northeast');
    end
    if i == 3, xlabel('Time [s]'); end
end
set(findall(fig2, 'Type', 'axes'), 'Toolbar', []);
exportgraphics(fig2, fullfile(fig_dir, '02_3sigma_consistency.png'), 'Resolution', 300, 'ContentType', 'image');

% -------------------------------------------------------------------------
% FIGURE 3: VELOCITÀ A REGIME VS 3-SIGMA
% -------------------------------------------------------------------------
fig3 = figure('Name', 'Velocity 3-Sigma Consistency', 'Color', 'w', 'Position', [100, 100, 900, 600]);
for i = 1:3
    subplot(3, 1, i);
    sig_v = sigma3_v(i, idx_ss);
    fill([t_ss; flipud(t_ss)], [sig_v(:); flipud(-sig_v(:))], [0.85 0.95 0.85], ...
         'EdgeColor', [0 0.5 0], 'LineStyle', '--', 'LineWidth', 1.0); hold on;
    plot(t_ss, err_vel_xyz(i, idx_ss), 'b', 'LineWidth', 0.6);
    grid on; ylabel(sprintf('v_%s Error [mm/s]', comp_names{i})); xlim([100, max(t_vec)]);
    if i == 1
        title('Steady-State Velocity Estimation Errors vs 3\sigma Envelope');
        legend('3\sigma Uncertainty Bound', 'Velocity Error', 'Location', 'northeast');
    end
    if i == 3, xlabel('Time [s]'); end
end
set(findall(fig3, 'Type', 'axes'), 'Toolbar', []);
exportgraphics(fig3, fullfile(fig_dir, '03_velocity_3sigma.png'), 'Resolution', 300, 'ContentType', 'image');

% -------------------------------------------------------------------------
% FIGURE 4: MAHALANOBIS GATING (OUTLIER REJECTION)
% -------------------------------------------------------------------------
D_m2 = zeros(1, N_samples);
for k = 1:N_samples
    y_k = y_innov(1:6, k);
    if ndims(S_innov) == 3, S_k = S_innov(1:6, 1:6, k); else, S_k = diag(S_innov(1:6, k)); end
    if norm(y_k) > 0 && det(S_k) > 0, D_m2(k) = y_k' * (S_k \ y_k); end
end

fig4 = figure('Name', 'EKF Outlier Rejection', 'Color', 'w', 'Position', [100, 100, 950, 550]);
subplot(2, 1, 1);
plot(t_vec, D_m2, 'b-', 'LineWidth', 1.0); hold on;
yline(18.55, 'r--', 'Threshold \gamma = 18.55 (6 DOF, 3\sigma)', 'LineWidth', 1.2);
grid on; ylabel('Squared Distance D_M^2'); title('Mahalanobis Distance Metric (GNSS Fault Test)');
xlim([0, max(t_vec)]);

subplot(2, 1, 2);
plot(t_vec, err_pos_norm, 'b-', 'LineWidth', 1.0); hold on;
xline(3000, 'k:', 'Injected Outlier (t = 3000s)', 'LineWidth', 1.2);
grid on; xlabel('Time [s]'); ylabel('Position Error [m]');
title('EKF Absolute Position Estimation Error'); ylim([0, 20]); xlim([0, max(t_vec)]);

set(findall(fig4, 'Type', 'axes'), 'Toolbar', []);
exportgraphics(fig4, fullfile(fig_dir, '04_outlier_rejection.png'), 'Resolution', 300, 'ContentType', 'image');

% -------------------------------------------------------------------------
% FIGURE 5: INNOVAZIONE / RESIDUI DI MISURA
% -------------------------------------------------------------------------
fig5 = figure('Name', 'Innovation Residuals', 'Color', 'w', 'Position', [100, 100, 900, 600]);
sigma3_S = zeros(3, N_samples);
for k = 1:N_samples
    if ndims(S_innov) == 3, S_k = S_innov(1:3, 1:3, k); else, S_k = diag(S_innov(1:3, k)); end
    sigma3_S(:, k) = 3 * sqrt(diag(S_k));
end

for i = 1:3
    subplot(3, 1, i);
    sig_s = sigma3_S(i, :);
    plot(t_vec, sig_s, 'r--', 'LineWidth', 1.0); hold on;
    plot(t_vec, -sig_s, 'r--', 'LineWidth', 1.0);
    plot(t_vec, y_innov(i, :), 'b', 'LineWidth', 0.5);
    grid on; ylabel(sprintf('y_%s [m]', comp_names{i})); xlim([0, max(t_vec)]);
    if i == 1
        title('Measurement Innovation Residuals vs $\pm 3\sqrt{S_k}$ Bounds', 'Interpreter', 'latex');
        legend('$3\sqrt{S_k}$ Uncertainty Bound', '', 'Innovation Residue', 'Location', 'northeast', 'Interpreter', 'latex');
    end
    if i == 3, xlabel('Time [s]'); end
end
set(findall(fig5, 'Type', 'axes'), 'Toolbar', []);
exportgraphics(fig5, fullfile(fig_dir, '05_innovation_residuals.png'), 'Resolution', 300, 'ContentType', 'image');

% -------------------------------------------------------------------------
% FIGURE 6: TEST STATISTICO NIS
% -------------------------------------------------------------------------
fig6 = figure('Name', 'NIS Test', 'Color', 'w', 'Position', [100, 100, 900, 400]);
plot(t_vec, D_m2, 'b', 'LineWidth', 0.8); hold on;
yline(12.59, 'r--', '\chi^2 Upper Threshold (\alpha=0.05, DOF=6)', 'LineWidth', 1.2);
grid on; xlabel('Time [s]'); ylabel('NIS Stat');
title('Normalized Innovation Squared (NIS) Consistency Test');
xlim([100, max(t_vec)]); ylim([0, 30]);

set(findall(fig6, 'Type', 'axes'), 'Toolbar', []);
exportgraphics(fig6, fullfile(fig_dir, '06_nis_consistency_test.png'), 'Resolution', 300, 'ContentType', 'image');

fprintf('\nSUCCESS: All 6 figures generated, cleaned and exported successfully in "%s"\n\n', fig_dir);