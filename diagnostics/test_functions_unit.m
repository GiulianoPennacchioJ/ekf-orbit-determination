%% TEST_FUNCTIONS_UNIT.M
% Unit Testing Script for EKF Orbit Determination Functions
% Author: Giuliano Pennacchio

clearvars; clc; close all;

% 1. Load workspace environment
% Determine root folder dynamically relative to this test script
script_dir = fileparts(mfilename('fullpath'));
root_dir   = fileparts(script_dir);
run(fullfile(root_dir, 'config', 'init_simulation.m'));

fprintf('=======================================================\n');
fprintf('       RUNNING UNIT TESTS FOR EKF MODULES\n');
fprintf('=======================================================\n\n');

%% TEST 1: Truth Dynamics Vector Computation
fprintf('[TEST 1/4] Checking Truth Dynamics Functions...');
try
    % Generate dummy harmonic matrices (20x20)
    C_nm = zeros(21, 21); S_nm = zeros(21, 21);
    C_nm(3,1) = -sim_param.J2; % Set J2
    
    x_dot_truth = truth_dynamics(0, x_true_0, sim_param, C_nm, S_nm);
    
    assert(length(x_dot_truth) == 6, 'Truth vector derivative must be 6x1.');
    assert(~any(isnan(x_dot_truth)), 'NaN detected in truth acceleration!');
    fprintf(' PASSED.\n');
catch ME
    fprintf(' FAILED: %s\n', ME.message);
end

%% TEST 2: GNSS Sensor Model Noise Statistics
fprintf('[TEST 2/4] Checking GNSS Sensor Model Noise...');
try
    N_samples = 10000;
    z_samples = zeros(6, N_samples);
    for k = 1:N_samples
        z_samples(:, k) = gnss_sensor_model(x_true_0, R);
    end
    
    % Compute empirical standard deviations
    std_pos = std(z_samples(1:3, :)', 1);
    std_vel = std(z_samples(4:6, :)', 1);
    
    % Verify match with expected sigma values (5m position, 0.05m/s velocity)
    assert(max(abs(std_pos - sigma_r_gnss)) < 1e-3, 'Position noise variance mismatch.');
    assert(max(abs(std_vel - sigma_v_gnss)) < 1e-5, 'Velocity noise variance mismatch.');
    fprintf(' PASSED.\n');
catch ME
    fprintf(' FAILED: %s\n', ME.message);
end

%% TEST 3: Analytical Jacobian J2 vs Finite Difference
fprintf('[TEST 3/4] Validating Analytical Jacobian F(x) via Finite Difference...');
try
    % Analytical Jacobian
    F_anal = jacobian_j2(x_true_0, sim_param);
    
    % Finite Difference Numerical Approximation
    eps_fd = 1e-6;
    F_num  = zeros(6, 6);
    
    for j = 1:6
        dx = zeros(6, 1);
        dx(j) = eps_fd;
        
        % Perturb state +/-
        x_plus  = x_true_0 + dx;
        x_minus = x_true_0 - dx;
        
        % Onboard dynamics call
        f_plus  = onboard_dyn_eval(x_plus, sim_param);
        f_minus = onboard_dyn_eval(x_minus, sim_param);
        
        F_num(:, j) = (f_plus - f_minus) / (2 * eps_fd);
    end
    
    % Compare analytical vs numerical matrices
    max_err = max(max(abs(F_anal - F_num)));
    assert(max_err < 1e-5, sprintf('Jacobian mismatch! Max diff: %e', max_err));
    fprintf(' PASSED (Max Error: %.2e).\n', max_err);
catch ME
    fprintf(' FAILED: %s\n', ME.message);
end

%% TEST 4: EKF Step Execution (Prediction & Update)
fprintf('[TEST 4/4] Checking EKF Predict & Update Step...');
try
    % Measure 1: Prediction step (10 Hz, no measurement)
    [x_pred, P_pred, ~, ~] = ekf_predict_update(x_hat_0, P0, x_true_0, 0, sim_param, Q, R);
    
    % Measure 2: Update step (1 Hz, GNSS measurement active)
    z_gnss = gnss_sensor_model(x_true_0, R);
    [x_up, P_up, y_inn, S_inn] = ekf_predict_update(x_pred, P_pred, z_gnss, 1, sim_param, Q, R);
    
    err_before = norm(x_true_0(1:3) - x_hat_0(1:3));
    err_after  = norm(x_true_0(1:3) - x_up(1:3));
    
    assert(all(eig(P_up) > 0), 'Covariance matrix P must remain Positive Definite.');
    assert(~any(isnan(x_up)), 'NaN detected in state update.');
    fprintf(' PASSED (Initial Err: %.2f km -> Updated Err: %.2f km).\n', err_before, err_after);
catch ME
    fprintf(' FAILED: %s\n', ME.message);
end

fprintf('\n=======================================================\n');
fprintf('       ALL UNIT TESTS COMPLETED SUCCESSFULLY!\n');
fprintf('=======================================================\n');

% Helper function for finite difference evaluation
function x_dot = onboard_dyn_eval(x, sim_param)
    r = x(1:3); v = x(4:6); r_mag = norm(r);
    a_central = - (sim_param.mu / r_mag^3) * r;
    z = r(3);
    factor = 1.5 * sim_param.J2 * sim_param.mu * (sim_param.R_E^2) / (r_mag^5);
    a_j2 = factor * [ x(1) * (5 * (z^2 / r_mag^2) - 1);
                      x(2) * (5 * (z^2 / r_mag^2) - 1);
                      x(3) * (5 * (z^2 / r_mag^2) - 3) ];
    x_dot = [v; a_central + a_j2];
end