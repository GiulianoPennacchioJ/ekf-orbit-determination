function F = jacobian_j2(x_hat, sim_param)
%% JACOBIAN_J2 Analytical Jacobian Matrix for Central Gravity + J2 Perturbation
% Inputs:
%   x_hat     - Estimated state vector [r_eci; v_eci] [6x1] (km, km/s)
%   sim_param - Struct containing mu, R_E, J2
% Output:
%   F         - Analytical Jacobian Matrix [6x6] (df/dx evaluated at x_hat)

    r = x_hat(1:3);
    x = r(1); y = r(2); z = r(3);
    r_mag = norm(r);
    
    mu  = sim_param.mu;
    RE  = sim_param.R_E;
    J2  = sim_param.J2;
    
    % Powers of radial distance
    r2 = r_mag^2;
    r3 = r_mag^3;
    r5 = r_mag^5;
    r7 = r_mag^7;
    r9 = r_mag^9;
    z2 = z^2;
    
    % --- 1. Central Gravity Part (da/dr_central) ---
    G_central = (3 * mu / r5) * (r * r') - (mu / r3) * eye(3);
    
    % --- 2. J2 Perturbation Part (da/dr_J2) ---
    % Pre-factor for J2 acceleration partial derivatives
    k = 1.5 * J2 * mu * RE^2;
    
    % Partial derivatives of J2 acceleration components w.r.t [x, y, z]
    daJ2_dx = k * [ (1 - 5*z2/r2)*(1/r5 - 5*x^2/r7) + 10*x^2*z2/r9 ;
                    (1 - 5*z2/r2)*(-5*x*y/r7) + 10*x*y*z2/r9 ;
                    (3 - 5*z2/r2)*(-5*x*z/r7) + 10*x*z^3/r9 - 10*x*z/r7 ];
                
    daJ2_dy = k * [ (1 - 5*z2/r2)*(-5*x*y/r7) + 10*x*y*z2/r9 ;
                    (1 - 5*z2/r2)*(1/r5 - 5*y^2/r7) + 10*y^2*z2/r9 ;
                    (3 - 5*z2/r2)*(-5*y*z/r7) + 10*y*z^3/r9 - 10*y*z/r7 ];
                
    daJ2_dz = k * [ -10*x*z/r7 + (1 - 5*z2/r2)*(-5*x*z/r7) + 10*x*z^3/r9 ;
                    -10*y*z/r7 + (1 - 5*z2/r2)*(-5*y*z/r7) + 10*y*z^3/r9 ;
                    -10*z/r7 + (3 - 5*z2/r2)*(1/r5 - 5*z2/r7) + 10*z^4/r9 - 10*z2/r7 ];
                
    G_J2 = [daJ2_dx, daJ2_dy, daJ2_dz];
    
    % --- 3. Assemble Full Jacobian Matrix (6x6) ---
    F = [zeros(3,3), eye(3,3);
         G_central + G_J2, zeros(3,3)];
end