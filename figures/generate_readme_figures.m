% =========================================================================
% MASTER SCRIPT: GENERATE README FIGURES FOR EKF ORBIT ESTIMATION
% =========================================================================

close all; clc; % Chiude solo i grafici aperti, NON cancella il Workspace!

% --- 1. SETUP OUTPUT DIRECTORY ---
fig_dir = 'figures';
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

% --- 2. LOAD SIMULATION DATA FROM WORKSPACE ---
if ~exist('out_x_true', 'var') || ~exist('out_x_hat', 'var')
    error('Simulation outputs missing. Please run Simulink model (ekf_orbit_sim.slx) first.');
end

x_true  = out_x_true.Data;
x_hat   = out_x_hat.Data;
y_innov = out_y_innov.Data;
S_innov = out_S_innov.Data;
t_vec   = out_x_true.Time;

% Compute Errors
err_pos_xyz  = (x_hat(1:3, :) - x_true(1:3, :)) * 1000; % [m]
err_pos_norm = sqrt(sum(err_pos_xyz.^2, 1));           % [m]

% Retrieve 3-sigma bounds from covariance P if available
if exist('out_P', 'var')
    P_data   = out_P.Data;
    sigma3_x = 3 * sqrt(squeeze(P_data(1, 1, :))) * 1000;
    sigma3_y = 3 * sqrt(squeeze(P_data(2, 2, :))) * 1000;
    sigma3_z = 3 * sqrt(squeeze(P_data(3, 3, :))) * 1000;
end

% --- FIGURE 1: INITIAL CONVERGENCE (0 - 30 s) ---
fig1 = figure('Name', 'EKF Initial Convergence', 'Color', 'w', 'Position', [100, 100, 900, 500]);
idx_30s = t_vec <= 30;

subplot(3, 1, 1);
plot(t_vec(idx_30s), err_pos_xyz(1, idx_30s), 'b-o', 'MarkerSize', 3); grid on;
ylabel('Error X [m]'); title('EKF Initial Transient Response (0 - 30s)');

subplot(3, 1, 2);
plot(t_vec(idx_30s), err_pos_xyz(2, idx_30s), 'b-o', 'MarkerSize', 3); grid on;
ylabel('Error Y [m]');

subplot(3, 1, 3);
plot(t_vec(idx_30s), err_pos_xyz(3, idx_30s), 'b-o', 'MarkerSize', 3); grid on;
ylabel('Error Z [m]'); xlabel('Time [s]');

exportgraphics(fig1, fullfile(fig_dir, '01_initial_convergence.png'), 'Resolution', 300);

% --- FIGURE 2: STEADY-STATE 3-SIGMA CONSISTENCY (t > 100 s) ---
if exist('out_P', 'var')
    fig2 = figure('Name', 'EKF 3-Sigma Bounds', 'Color', 'w', 'Position', [100, 100, 950, 600]);
    idx_ss = t_vec >= 100;
    t_ss   = t_vec(idx_ss);
    sigs   = {sigma3_x(idx_ss), sigma3_y(idx_ss), sigma3_z(idx_ss)};
    labels = {'Error X [m]', 'Error Y [m]', 'Error Z [m]'};

    for i = 1:3
        subplot(3, 1, i);
        err_i = err_pos_xyz(i, idx_ss);
        sig_i = sigs{i};
        
        fill([t_ss; flipud(t_ss)], [sig_i(:); flipud(-sig_i(:))], [1 0.8 0.8], ...
             'EdgeColor', [0.8 0 0], 'LineStyle', '--', 'LineWidth', 1.0); hold on;
        plot(t_ss, err_i, 'b', 'LineWidth', 0.6);
        grid on; ylabel(labels{i}); ylim([-25, 25]); xlim([100, max(t_vec)]);
        
        if i == 1
            title('EKF Steady-State Consistency vs 3\sigma Envelope');
            legend('3\sigma Uncertainty Bounds', 'Estimation Error', 'Location', 'northeast');
        end
        if i == 3, xlabel('Time [s]'); end
    end
    exportgraphics(fig2, fullfile(fig_dir, '02_3sigma_consistency.png'), 'Resolution', 300);
end

% --- FIGURE 3: MAHALANOBIS GATING & OUTLIER REJECTION ---
N_samples = length(t_vec);
D_m2      = zeros(1, N_samples);
for k = 1:N_samples
    y_k = y_innov(1:6, k);
    if ndims(S_innov) == 3
        S_k = S_innov(1:6, 1:6, k);
    else
        S_k = diag(S_innov(1:6, k));
    end
    if norm(y_k) > 0 && det(S_k) > 0
        D_m2(k) = y_k' * (S_k \ y_k);
    end
end

fig3 = figure('Name', 'EKF Outlier Rejection', 'Color', 'w', 'Position', [100, 100, 950, 550]);
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

exportgraphics(fig3, fullfile(fig_dir, '03_outlier_rejection.png'), 'Resolution', 300);

fprintf('\n SUCCESS: All 3 high-resolution figures saved in "%s/" folder.\n\n', fig_dir);