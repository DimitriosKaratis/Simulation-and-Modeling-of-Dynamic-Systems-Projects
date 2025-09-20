%% SIMULATION AND MODELING OF DYNAMIC SYSTEMS
% Assignment 2 - April 2025
% KARATIS DIMITRIOS 10775

%% Exercise 2 - b
% System Dynamics Function
function dx = roll_angle_dynamics(t, x, a1, a2, a3, b, ...
    phi_0, phi_inf, lambda, rho, k1, k2, r_d)

    % Extract current state variables
    r = x(1);           % Position
    r_dot = x(2);       % Velocity

    % Desired trajectory at time t
    r_des = r_d(t);

    % Compute time-varying gain phi(t)
    phi = (phi_0 - phi_inf) * exp(-lambda * t) + phi_inf;

    % First tracking error (normalized by phi)
    z1 = (r - r_des) / phi;

    % Feedback law based on first error
    a_term = -k1 * log((1 + z1) / (1 - z1));

    % Second tracking error (based on velocity deviation)
    z2 = (r_dot - a_term) / rho;

    % Final control input using second feedback term
    u = -k2 * log((1 + z2) / (1 - z2));
    
    % Zero disturbance  
    dist = 0; 

    % Initialize derivative of state vector
    dx = zeros(2,1);
    
    % First state derivative: dr/dt = r_dot
    dx(1) = r_dot;

    % Second state derivative: full nonlinear dynamics with input u and disturbance
    dx(2) = -a1 * r_dot - a2 * sin(r) + a3 * r_dot^2 * sin(2 * r) + b * u + dist;
end