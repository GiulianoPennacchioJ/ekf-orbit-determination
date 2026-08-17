function z_gnss = gnss_sensor_model(x_true, R)
%% GNSS_SENSOR_MODEL Emulates Spaceborne GNSS Receiver Observations
% Inputs:
%   x_true - True Cartesian state vector [r_eci; v_eci] [6x1] (km, km/s)
%   R      - Measurement noise covariance matrix [6x6]
% Output:
%   z_gnss - Noisy GNSS measurement vector [r_gnss; v_gnss] [6x1] (km, km/s)

    % Extract standard deviations from diagonal of R matrix
    sigma_vector = sqrt(diag(R)); % [3x pos std dev; 3x vel std dev]
    
    % Generate Gaussian zero-mean white noise vector
    measurement_noise = sigma_vector .* randn(6, 1);
    
    % Corrupt true state with measurement noise
    z_gnss = x_true + measurement_noise;
end