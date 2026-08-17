function a_drag = atmospheric_drag(r_eci, v_eci, sim_param)
%% ATMOSPHERIC_DRAG Exponential Atmospheric Drag Model
% Inputs:
%   r_eci, v_eci - Spacecraft position [km] and velocity [km/s] in ECI
%   sim_param    - Struct containing vehicle parameters (Cd, A, mass)
% Output:
%   a_drag       - Drag acceleration in ECI frame [km/s^2] [3x1]

    % Spacecraft ballistic properties
    Cd   = 2.2;                  % Drag coefficient (-)
    Area = 2.0e-6;               % Cross-sectional area [km^2] (2 m^2)
    mass = 500.0;                % Mass [kg]
    
    % Altitude computation
    r_mag = norm(r_eci);
    alt   = r_mag - sim_param.R_E; % [km]
    
    % Exponential Density Model (Base parameters for 600 km)
    rho_0 = 1.454e-4;            % Reference density at 600 km [kg/km^3]
    h_0   = 600.0;               % Reference altitude [km]
    H     = 71.8;                % Scale height [km]
    
    rho = rho_0 * exp(-(alt - h_0) / H); % Atmospheric density [kg/km^3]
    
    % Relative velocity considering Earth rotation
    v_rel = v_eci - cross([0; 0; sim_param.omega_E], r_eci);
    v_rel_mag = norm(v_rel);
    
    % Drag Acceleration
    a_drag = -0.5 * rho * (Cd * Area / mass) * v_rel_mag * v_rel;
end