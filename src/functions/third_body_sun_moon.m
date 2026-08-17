function a_3b = third_body_sun_moon(r_eci, t)
%% THIRD_BODY_SUN_MOON Solar and Lunar Gravitational Perturbations
% Inputs:
%   r_eci - Spacecraft position in ECI frame [km] [3x1]
%   t     - Current simulation time [s]
% Output:
%   a_3b  - Net third-body acceleration in ECI [km/s^2] [3x1]

    % Constants
    mu_sun  = 1.32712440018e11; % Gravitational parameter Sun [km^3/s^2]
    mu_moon = 4902.8000;        % Gravitational parameter Moon [km^3/s^2]
    
    % Analytical Ephemerides Approximation for Sun and Moon ECI vectors
    r_sun  = get_sun_position(t);   % [km]
    r_moon = get_moon_position(t);  % [km]
    
    % Differential acceleration Sun
    r_s_sat = r_sun - r_eci;
    a_sun   = mu_sun * (r_s_sat / norm(r_s_sat)^3 - r_sun / norm(r_sun)^3);
    
    % Differential acceleration Moon
    r_m_sat = r_moon - r_eci;
    a_moon  = mu_moon * (r_m_sat / norm(r_m_sat)^3 - r_moon / norm(r_moon)^3);
    
    a_3b = a_sun + a_moon;
end

function r_sun = get_sun_position(t)
    % Simplified analytical solar position model
    omega_sun = 1.990986e-7; % [rad/s]
    AU = 149597870.7;        % [km]
    obliquity = deg2rad(23.439);
    
    lambda_s = omega_sun * t;
    r_sun = AU * [cos(lambda_s); sin(lambda_s) * cos(obliquity); sin(lambda_s) * sin(obliquity)];
end

function r_moon = get_moon_position(t)
    % Simplified analytical lunar position model
    omega_moon = 2.661866e-6; % [rad/s]
    r_m_mag    = 384400.0;     % [km]
    inc_moon   = deg2rad(5.145);
    
    lambda_m = omega_moon * t;
    r_moon = r_m_mag * [cos(lambda_m); sin(lambda_m) * cos(inc_moon); sin(lambda_m) * sin(inc_moon)];
end