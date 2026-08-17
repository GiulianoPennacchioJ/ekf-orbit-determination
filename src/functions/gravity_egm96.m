function a_grav = gravity_egm96(r_eci, t, sim_param, C_nm, S_nm, max_degree)
%% GRAVITY_EGM96 High-Degree Spherical Harmonic Gravity Model
% Inputs:
%   r_eci      - Spacecraft position in ECI frame [km] [3x1]
%   t          - Current simulation time [s]
%   sim_param  - Struct containing mu, R_E, omega_E
%   C_nm, S_nm - Spherical harmonic coefficients matrices (normalized or unnormalized)
%   max_degree - Maximum degree/order truncation (e.g., 20)
% Output:
%   a_grav     - Gravitational acceleration in ECI frame [km/s^2] [3x1]

    % Convert ECI position to ECEF frame via Greenwich Mean Sidereal Time (GMST)
    theta_gmst = sim_param.omega_E * t;
    R_eci2ecef = [ cos(theta_gmst),  sin(theta_gmst), 0;
                  -sin(theta_gmst),  cos(theta_gmst), 0;
                                 0,                0, 1];
                             
    r_ecef = R_eci2ecef * r_eci;
    r_mag  = norm(r_ecef);
    
    % Spherical coordinates in ECEF
    phi   = asin(r_ecef(3) / r_mag);        % Geocentric latitude [rad]
    lambda = atan2(r_ecef(2), r_ecef(1));   % Longitude [rad]
    
    % Legendre Polynomials P_nm(sin(phi))
    P = legendre_polynomials(max_degree, sin(phi));
    
    % Spherical harmonic potential derivatives (Radial, Latitudinal, Longitudinal)
    dU_dr     = -sim_param.mu / r_mag^2;
    dU_dphi   = 0;
    dU_dlambda = 0;
    
    for n = 2:max_degree
        r_ratio = (sim_param.R_E / r_mag)^n;
        term_r  = 0;
        term_phi = 0;
        term_lam = 0;
        
        for m = 0:n
            C = C_nm(n+1, m+1);
            S = S_nm(n+1, m+1);
            trig_sum = C * cos(m * lambda) + S * sin(m * lambda);
            trig_diff = -C * sin(m * lambda) + S * cos(m * lambda);
            
            term_r   = term_r + P(n+1, m+1) * trig_sum;
            term_phi = term_phi + (P(n+1, m+2) - m * tan(phi) * P(n+1, m+1)) * trig_sum;
            term_lam = term_lam + m * P(n+1, m+1) * trig_diff;
        end
        
        dU_dr      = dU_dr - (sim_param.mu / r_mag^2) * (n + 1) * r_ratio * term_r;
        dU_dphi    = dU_dphi + (sim_param.mu / r_mag) * r_ratio * term_phi;
        dU_dlambda = dU_dlambda + (sim_param.mu / r_mag) * r_ratio * term_lam;
    end
    
    % Convert spherical gradient components to Cartesian ECEF acceleration
    a_ecef = [ (dU_dr * cos(phi) * cos(lambda) - dU_dphi * sin(phi) * cos(lambda) - dU_dlambda * sin(lambda) / cos(phi));
               (dU_dr * cos(phi) * sin(lambda) - dU_dphi * sin(phi) * sin(lambda) + dU_dlambda * cos(lambda) / cos(phi));
               (dU_dr * sin(phi) + dU_dphi * cos(phi)) ];
           
    % Transform acceleration back to ECI frame
    a_grav = R_eci2ecef' * a_ecef;
end

function P = legendre_polynomials(n_max, x)
    % Compute Associated Legendre Polynomials P_nm(x)
    P = zeros(n_max + 2, n_max + 2);
    P(1, 1) = 1; % P_0,0
    
    for n = 1:n_max
        P(n+1, 1) = legendre_p_single(n, 0, x);
        for m = 1:n
            P(n+1, m+1) = legendre_p_single(n, m, x);
        end
    end
end

function val = legendre_p_single(n, m, x)
    P_all = legendre(n, x);
    val = P_all(m + 1);
end