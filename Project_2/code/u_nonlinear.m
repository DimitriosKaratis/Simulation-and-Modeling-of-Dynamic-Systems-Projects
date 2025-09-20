%% SIMULATION AND MODELING OF DYNAMIC SYSTEMS
% Assignment 2 - April 2025
% KARATIS DIMITRIOS 10775

%% Exercise 2 - a
% Controller Function
function u = u_nonlinear(t, x, r_d_fun, phi0, phi_inf, lambda, rho, k1, k2, T)
    % Compute phi(t)
    phi = (phi0 - phi_inf)*exp(-lambda*t) + phi_inf;

    % Normalized position error
    z1 = (x(1) - r_d_fun(t)) / phi;

    % Virtual control law alpha
    alpha = -k1 * T(z1);

    % Normalized velocity error
    z2 = (x(2) - alpha) / rho;

    % Final control input
    u = -k2 * T(z2);
end