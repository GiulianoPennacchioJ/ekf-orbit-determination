%% INIT_SIMULATION.M
% =========================================================================
% Project: Extended Kalman Filter (EKF) for LEO Orbit Determination
% Author: Giuliano Pennacchio - Flight Dynamics & GNC Engineer
% Description: Master initialization script for the Simulink model 'ekf_orbit_sim.slx'.
%              Defines astrodynamic constants, Sun-Synchronous Orbit (SSO)
%              initial conditions, truth vs. estimated states, discrete-continuous
%              execution rates, and EKF covariance matrices (P0, Q, R).
% =========================================================================

clearvars; clc; close all;

% Add project root and all subdirectories to MATLAB search path
project_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(project_root));

%% 1. ASTRODYNAMIC CONSTANTS & REFERENCE FRAME (WGS-84 / EGM96)
sim_param.mu      = 398600.4418;       % Standard gravitational parameter [km^3/s^2]
sim_param.R_E     = 6378.137;          % Equatorial Earth radius [km]
sim_param.J2      = 1.08262668e-3;     % Second zonal harmonic coefficient (-)
sim_param.omega_E = 7.2921159e-5;      % Earth rotation rate [rad/s]

%% 2. TARGET ORBITAL PARAMETERS (Sun-Synchronous Orbit - SSO)
sim_param.alt     = 600.0;                        % Altitude above sea level [km]
sim_param.a       = sim_param.R_E + sim_param.alt; % Semi-major axis [km]
sim_param.e       = 0.0001;                       % Near-circular eccentricity (-)

% Correct analytical SSO inclination computation for h = 600 km (i ≈ 97.79 deg)
dOmega_sun        = 0.9856 * (pi / 180) / 86400;  % Mean solar precession rate [rad/s]
cos_i             = -2/3 * (dOmega_sun / sim_param.J2) * ...
                    (1 / (sim_param.R_E^2 * sqrt(sim_param.mu))) * ...
                    (sim_param.a^(7/2));
sim_param.inc     = acos(cos_i);                  % Orbital inclination [rad] (~1.7067 rad = 97.79 deg)

sim_param.RAAN    = deg2rad(45.0);                % Right Ascension of Ascending Node [rad]
sim_param.argp    = deg2rad(0.0);                 % Argument of perigee [rad]
sim_param.nu      = deg2rad(0.0);                 % Initial true anomaly [rad]

% Orbital period and circular orbital speed
sim_param.T_orbit = 2 * pi * sqrt(sim_param.a^3 / sim_param.mu); % Orbital period [s]
sim_param.v_circ  = sqrt(sim_param.mu / sim_param.a);          % Circular velocity [km/s]

%% 3. TRUTH INITIAL STATE COMPUTATION (x_true_0 in ECI Frame)
% Transformation from Keplerian orbital elements to Cartesian ECI state vector
[r_eci_0, v_eci_0] = oe2eci(sim_param.a, sim_param.e, sim_param.inc, ...
                             sim_param.RAAN, sim_param.argp, sim_param.nu, sim_param.mu);

x_true_0 = [r_eci_0; v_eci_0]; % Initial true state vector [6x1] (km, km/s)

%% 4. ESTIMATED INITIAL STATE (x_hat_0 for EKF Initialization)
% Representative On-Board Computer (OBC) boot/injection offset:
% - Position offset: 1.0 km per axis
% - Velocity offset: 0.001 km/s (1 m/s) per axis
delta_r_0 = [1.0; -1.0; 0.5];          % Position estimation error vector [km]
delta_v_0 = [0.001; -0.001; 0.0005];   % Velocity estimation error vector [km/s]

x_hat_0 = x_true_0 + [delta_r_0; delta_v_0]; % Initial estimated state vector [6x1]

%% 5. SIMULATION TIMING & SAMPLE RATES
sim_param.dt_sim   = 0.1;              % Continuous integration step / EKF Prediction rate [s] (10 Hz)
sim_param.dt_gnss  = 1.0;              % GNSS measurement period / EKF Update rate [s] (1 Hz)
sim_param.t_final  = 2 * sim_param.T_orbit; % Total simulation duration (2 full orbits ~ 3.2 hours) [s]

%% 6. EKF COVARIANCE MATRICES CONFIGURATION (P0, Q, R)

% 6.1 Initial Error Covariance Matrix (P0)
% Configured strictly consistent with initial injection state offsets (sigma^2)
sigma_r0 = 1.0;                        % Initial position standard deviation [km]
sigma_v0 = 0.001;                      % Initial velocity standard deviation [km/s] (1 m/s)

P0 = diag([sigma_r0^2 * ones(1,3), sigma_v0^2 * ones(1,3)]);

% 6.2 GNSS Measurement Noise Covariance Matrix (R)
% Nominal spaceborne LEO GNSS receiver performance metrics:
% - Position 1-sigma RMS: 5.0 meters = 0.005 km -> sigma_r_gnss^2 = 2.5e-5 km^2
% - Velocity 1-sigma RMS: 0.05 m/s = 5e-5 km/s -> sigma_v_gnss^2 = 2.5e-9 (km/s)^2
sigma_r_gnss = 0.005;                  % Position noise std dev [km] (5 m)
sigma_v_gnss = 5.0e-5;                 % Velocity noise std dev [km/s] (0.05 m/s)

R = diag([sigma_r_gnss^2 * ones(1,3), sigma_v_gnss^2 * ones(1,3)]);

% 6.3 Process Noise Covariance Matrix (Q) - RETUNED FOR STEP 6
% Compensates for unmodeled high-degree harmonics (EGM96 20x20 vs J2 EKF model)
% Unit: Position [km^2], Velocity [(km/s)^2]

q_pos = 1.0e-8;  % Kinematic process noise spectral density [km^2]
q_vel = 1.0e-5;  % Dynamic noise: accounts for ~3.16 mm/s^2 unmodeled acceleration [(km/s)^2]

Q = diag([q_pos * ones(1,3), q_vel * ones(1,3)]);

% Pre-allocate 21x21 matrices for gravity field model
C_nm = zeros(21, 21);
S_nm = zeros(21, 21);

% Set J2 coefficient directly into C(2,0) position (degree 2, order 0)
C_nm(3,1) = -sim_param.J2; % In MATLAB indexing: row = degree+1, col = order+1

%% 7. AUXILIARY FUNCTION: OE2ECI (Keplerian Elements to ECI Cartesian Transformation)
function [r_eci, v_eci] = oe2eci(a, e, inc, RAAN, argp, nu, mu)
    % Compute magnitude and perifocal state vector
    p = a * (1 - e^2);
    r_val = p / (1 + e * cos(nu));
    
    r_perifocal = [r_val * cos(nu); r_val * sin(nu); 0];
    v_perifocal = sqrt(mu / p) * [-sin(nu); e + cos(nu); 0];
    
    % Direction Cosine Matrix (DCM) from Perifocal to ECI frame (3-1-3 Euler Angles)
    R3_RAAN = [cos(RAAN) -sin(RAAN) 0; sin(RAAN) cos(RAAN) 0; 0 0 1];
    R1_inc  = [1 0 0; 0 cos(inc) -sin(inc); 0 sin(inc) cos(inc)];
    R3_argp = [cos(argp) -sin(argp) 0; sin(argp) cos(argp) 0; 0 0 1];
    
    Q_p2eci = R3_RAAN * R1_inc * R3_argp;
    
    % Frame transformation
    r_eci = Q_p2eci * r_perifocal;
    v_eci = Q_p2eci * v_perifocal;
end

disp('===========================================================');
disp(' Initialization completed successfully.');
disp([' Target Orbit: LEO SSO at h = ', num2str(sim_param.alt), ' km']);
disp([' Calculated Inclination: ', num2str(rad2deg(sim_param.inc)), ' deg']);
disp([' Orbital Period: ', num2str(sim_param.T_orbit / 60), ' minutes']);
disp(' Base Workspace updated: sim_param, x_true_0, x_hat_0, P0, Q, R');
disp('===========================================================');