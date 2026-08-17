function x_dot = truth_dynamics(t, x_true, sim_param, C_nm, S_nm)
%% TRUTH_DYNAMICS Master Truth Simulator Vector Differential Equation
% Inputs:
%   t          - Current simulation time [s]
%   x_true     - State vector [r_eci; v_eci] [6x1]
%   sim_param  - Simulation parameters struct
%   C_nm, S_nm - EGM96 Spherical Harmonic matrices
% Output:
%   x_dot      - State derivative [v_eci; a_total] [6x1]

    r_eci = x_true(1:3);
    v_eci = x_true(4:6);
    
    % 1. High-Degree Gravity Acceleration (EGM96 20x20)
    max_degree = 20;
    a_grav = gravity_egm96(r_eci, t, sim_param, C_nm, S_nm, max_degree);
    
    % 2. Third-Body Acceleration (Sun & Moon)
    a_3b = third_body_sun_moon(r_eci, t);
    
    % 3. Atmospheric Drag Acceleration
    a_drag = atmospheric_drag(r_eci, v_eci, sim_param);
    
    % 4. Solar Radiation Pressure (SRP)
    a_srp = srp_force(r_eci, t, sim_param);
    
    % Net Acceleration
    a_total = a_grav + a_3b + a_drag + a_srp;
    
    x_dot = [v_eci; a_total];
end