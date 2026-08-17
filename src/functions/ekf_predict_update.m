function [x_hat_out, P_out, y_innov, S_innov] = ekf_predict_update(x_hat_in, P_in, z_gnss, meas_available, sim_param, Q, R)
%% EKF_PREDICT_UPDATE Discrete-Continuous Extended Kalman Filter Engine
% Inputs:
%   x_hat_in       - Previous state estimate [6x1]
%   P_in           - Previous error covariance matrix [6x6]
%   z_gnss         - Current GNSS measurement vector [6x1]
%   meas_available - Boolean flag (1 when 1 Hz GNSS measurement is active, 0 otherwise)
%   sim_param      - Simulation parameters struct
%   Q              - Process noise covariance matrix [6x6]
%   R              - Measurement noise covariance matrix [6x6]
% Outputs:
%   x_hat_out      - Updated state estimate [6x1]
%   P_out          - Updated error covariance matrix [6x6]
%   y_innov        - Innovation vector [6x1] (zeros if no measurement)
%   S_innov        - Innovation covariance matrix [6x6] (zeros if no measurement)

    dt = sim_param.dt_sim; % Prediction step size (0.1 s = 10 Hz)
    
    %% 1. PREDICTION STEP (Always executed at 10 Hz)
    % 1.1 State Propagation using On-Board Dynamics (mu + J2)
    x_dot = onboard_dynamics_j2(x_hat_in, sim_param);
    x_hat_priori = x_hat_in + x_dot * dt; % Euler integration (or RK2/RK4)
    
    % 1.2 Covariance Propagation
    F = jacobian_j2(x_hat_in, sim_param);
    Phi = eye(6) + F * dt; % First-order State Transition Matrix approximation
    P_priori = Phi * P_in * Phi' + Q;
    
    %% 2. CORRECTION / UPDATE STEP (Executed at 1 Hz when meas_available == 1)
    if meas_available
        % Measurement Matrix (Direct state observation)
        H = eye(6);
        
        % Innovation computation
        y_innov = z_gnss - x_hat_priori;
        S_innov = H * P_priori * H' + R;
        
        % Kalman Gain computation
        K = (P_priori * H') / S_innov;
        
        % State and Covariance Update (Joseph form for numerical stability)
        x_hat_out = x_hat_priori + K * y_innov;
        I_KH = eye(6) - K * H;
        P_out = I_KH * P_priori * I_KH' + K * R * K';
    else
        % No measurement available: Priori state and covariance passed through
        x_hat_out = x_hat_priori;
        P_out     = P_priori;
        y_innov   = zeros(6, 1);
        S_innov   = zeros(6, 6);
    end
end

function x_dot = onboard_dynamics_j2(x, sim_param)
    % On-board simplified dynamics: Central Gravity + J2
    r = x(1:3);
    v = x(4:6);
    r_mag = norm(r);
    
    % Central Gravity
    a_central = - (sim_param.mu / r_mag^3) * r;
    
    % J2 Perturbation
    z = r(3);
    factor = 1.5 * sim_param.J2 * sim_param.mu * (sim_param.R_E^2) / (r_mag^5);
    a_j2 = factor * [ x(1) * (5 * (z^2 / r_mag^2) - 1);
                      x(2) * (5 * (z^2 / r_mag^2) - 1);
                      x(3) * (5 * (z^2 / r_mag^2) - 3) ];
                  
    a_total = a_central + a_j2;
    x_dot = [v; a_total];
end