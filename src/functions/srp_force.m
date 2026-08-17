function a_srp = srp_force(r_eci, t, sim_param)
%% SRP_FORCE Solar Radiation Pressure Cannonball Model
% Inputs:
%   r_eci     - Spacecraft position [km]
%   t         - Time [s]
%   sim_param - Parameters struct
% Output:
%   a_srp     - SRP acceleration in ECI [km/s^2] [3x1]

    Cr    = 1.3;               % Radiation pressure coefficient (-)
    Area  = 2.0e-6;            % Cross-sectional area [km^2]
    mass  = 500.0;             % Mass [kg]
    P_sun = 4.56e-6;           % Solar radiation pressure at 1 AU [N/m^2] -> [kg/(km*s^2)]
    P_sun = P_sun * 1.0e-3;    % Convert to [kg/(km*s^2)]
    
    r_sun     = get_sun_position(t);
    r_sun_sat = r_eci - r_sun;
    u_sun     = r_sun_sat / norm(r_sun_sat);
    
    % Cylindrical Shadow Function (1 = Sunlit, 0 = Eclipse)
    nu_shadow = shadow_function(r_eci, r_sun, sim_param.R_E);
    
    a_srp = -nu_shadow * P_sun * Cr * (Area / mass) * u_sun;
end

function nu = shadow_function(r_sat, r_sun, R_E)
    % Simple cylindrical Earth eclipse check
    u_sun = r_sun / norm(r_sun);
    proj  = dot(r_sat, u_sun);
    
    if proj > 0
        nu = 1.0; % Sunlit
    else
        d_perp = norm(r_sat - proj * u_sun);
        if d_perp < R_E
            nu = 0.0; % Eclipse
        else
            nu = 1.0; % Sunlit
        end
    end
end

function r_sun = get_sun_position(t)
    omega_sun = 1.990986e-7;
    AU        = 149597870.7;
    obliquity = deg2rad(23.439);
    lambda_s  = omega_sun * t;
    r_sun     = AU * [cos(lambda_s); sin(lambda_s) * cos(obliquity); sin(lambda_s) * sin(obliquity)];
end