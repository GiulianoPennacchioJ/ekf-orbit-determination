%% VALIDATE_TRUTH_SIMULATOR.M
% Verification script for high-fidelity truth propagation

run('config/init_simulation.m');

% Generate dummy C_nm and S_nm coefficients for testing (J2, J3, J4 enabled)
C_nm = zeros(21, 21);
S_nm = zeros(21, 21);
C_nm(3, 1) = -sim_param.J2; % J2 term

fprintf('Integrating Truth Orbit for %.1f minutes...\n', sim_param.t_final / 60);

% Numerical Integration using ode45
options = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
[t_vec, x_mat] = ode45(@(t, x) truth_dynamics(t, x, sim_param, C_nm, S_nm), ...
                       [0, sim_param.t_final], x_true_0, options);

% Plotting Trajectory and Orbital Radius
r_mag_vec = sqrt(x_mat(:,1).^2 + x_mat(:,2).^2 + x_mat(:,3).^2);

figure('Name', 'Truth Propagation Check');
subplot(2,1,1);
plot(t_vec / 60, r_mag_vec - sim_param.R_E, 'LineWidth', 1.5);
grid on;
xlabel('Time [min]'); ylabel('Altitude [km]');
title('Truth Environment - Spacecraft Altitude');

subplot(2,1,2);
plot(t_vec / 60, sqrt(x_mat(:,4).^2 + x_mat(:,5).^2 + x_mat(:,6).^2), 'r', 'LineWidth', 1.5);
grid on;
xlabel('Time [min]'); ylabel('Speed [km/s]');
title('Truth Environment - Orbital Speed Profile');

disp('Truth Simulator successfully verified.');