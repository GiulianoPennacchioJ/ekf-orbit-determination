function a_grav = gravity_egm96(r_eci, t, sim_param, C_nm, S_nm, max_degree)
%% GRAVITY_EGM96 Spherical Harmonics Gravity Model (Code-Gen Compatible)
% Computes gravitational acceleration usingAssociated Legendre Functions via direct recursion.

    mu = sim_param.mu;
    RE = sim_param.R_E;
    
    r_mag = norm(r_eci);
    x = r_eci(1); y = r_eci(2); z = r_eci(3);
    
    % Latitude and Longitude calculation
    phi = asin(z / r_mag);       % Geocentric latitude
    lambda = atan2(y, x);        % Longitude
    
    u = sin(phi);                % sin(phi)
    cos_phi = cos(phi);
    if abs(cos_phi) < 1e-12
        cos_phi = 1e-12;        % Avoid singularity at poles
    end
    
    % Precompute powers of (R_E / r)
    ratio = RE / r_mag;
    
    % Compute Associated Legendre Polynomials P_nm(u) via pure recursion
    P = compute_legendre_pure(max_degree, u);
    
    % Precompute cos(m*lambda) and sin(m*lambda)
    m_vec = (0:max_degree)';
    cos_m_lambda = cos(m_vec * lambda);
    sin_m_lambda = sin(m_vec * lambda);
    
    % Radial, latitudinal, and longitudinal acceleration components in Spherical Coords
    dU_dr     = -mu / r_mag^2;
    dU_dphi   = 0;
    dU_dlambda = 0;
    
    term_r = 0;
    term_phi = 0;
    term_lambda = 0;
    
    for n = 2:max_degree
        ratio_n = ratio^n;
        
        sum_r = 0;
        sum_phi = 0;
        sum_lambda = 0;
        
        for m = 0:n
            C = C_nm(n+1, m+1);
            S = S_nm(n+1, m+1);
            
            if C == 0 && S == 0
                continue;
            end
            
            P_nm = P(n+1, m+1);
            
            % Derivative dP_nm/dphi using recursion relation
            if m == n
                P_n_mplus1 = 0;
            else
                P_n_mplus1 = P(n+1, m+2);
            end
            dP_dphi = m * tan(phi) * P_nm + P_n_mplus1;
            
            % Harmonic combinations
            A_nm = C * cos_m_lambda(m+1) + S * sin_m_lambda(m+1);
            B_nm = S * cos_m_lambda(m+1) - C * sin_m_lambda(m+1);
            
            sum_r      = sum_r + (n + 1) * P_nm * A_nm;
            sum_phi    = sum_phi + dP_dphi * A_nm;
            sum_lambda = sum_lambda + m * P_nm * B_nm;
        end
        
        term_r      = term_r + ratio_n * sum_r;
        term_phi    = term_phi + ratio_n * sum_phi;
        term_lambda = term_lambda + ratio_n * sum_lambda;
    end
    
    dU_dr      = dU_dr + (mu / r_mag^2) * term_r;
    dU_dphi    = (mu / r_mag) * term_phi;
    dU_dlambda = (mu / r_mag) * term_lambda;
    
    % Transform Spherical Gradient (dU/dr, dU/dphi, dU/dlambda) to ECI Cartesian Frame
    r_xy = sqrt(x^2 + y^2);
    if r_xy < 1e-12
        r_xy = 1e-12;
    end
    
    % Unit Vectors
    e_r = [x; y; z] / r_mag;
    e_phi = [-x*z / (r_mag*r_xy); -y*z / (r_mag*r_xy); r_xy / r_mag];
    e_lambda = [-y / r_xy; x / r_xy; 0];
    
    a_grav = dU_dr * e_r + (1 / r_mag) * dU_dphi * e_phi + (1 / (r_mag * cos_phi)) * dU_dlambda * e_lambda;
end

function P = compute_legendre_pure(N, u)
%% COMPUTE_LEGENDRE_PURE Associated Legendre Polynomials Recursion (No 'legendre' call)
% Returns (N+1)x(N+1) matrix P where P(n+1, m+1) = P_nm(u)

    P = zeros(N+1, N+1);
    
    % Base values
    P(1,1) = 1.0; % P_0,0
    
    if N == 0
        return;
    end
    
    cos_phi = sqrt(max(0, 1 - u^2));
    
    % Sectorial terms P_m,m (m = n)
    for m = 1:N
        P(m+1, m+1) = (2*m - 1) * cos_phi * P(m, m);
    end
    
    % Terms P_m+1,m
    for m = 0:N-1
        P(m+2, m+1) = (2*m + 1) * u * P(m+1, m+1);
    end
    
    % Full recursion for remaining P_n,m
    for m = 0:N-2
        for n = m+2:N
            P(n+1, m+1) = ((2*n - 1) * u * P(n, m+1) - (n + m - 1) * P(n-1, m+1)) / (n - m);
        end
    end
end