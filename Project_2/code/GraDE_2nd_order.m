%% SIMULATION AND MODELING OF DYNAMIC SYSTEMS
% Assignment 2 - April 2025
% KARATIS DIMITRIOS 10775

%% Exercise 1 - a
% Gradient Descent Estimator for mx'' + bx' + kx = u
function xxdot = GraDE_2nd_order(t, xx)
    global gamma m b k c d f p1 p2

    % States
    x = xx(1);
    x_dot = xx(2);
    j1 = xx(3);
    j2 = xx(4);
    j3 = xx(5);
    dj1 = xx(6);
    dj2 = xx(7);
    dj3 = xx(8);
    theta1_est = xx(9);
    theta2_est = xx(10);
    theta3_est = xx(11);
    
    lambda1 = p1 + p2;
    lambda2 = p1 * p2;

    % Input
    u = c * sin(d * t) + f;

    theta_true = [-b/m, -k/m, 1/m];

    % Compute estimated acceleration
    x_est = theta1_est * j1 + theta2_est * j2 + theta3_est * j3;
    x_ddot_true = theta_true(1) * x_dot + theta_true(2) * x + theta_true(3) * u;

    % Modeling error
    e = x - x_est;

    % xxdot = zeros(size(xx));
    xxdot = zeros(11,1);

    % System dynamics
    xxdot(1) = x_dot; 
    xxdot(2) = x_ddot_true;

    % Filters (dynamics of regressors)
    xxdot(3) = dj1;
    xxdot(4) = dj2;
    xxdot(5) = dj3;

    xxdot(6) = - dj1 * lambda1 - j1 * lambda2 - x_dot;
    xxdot(7) = - dj2 * lambda1 - j2 * lambda2 - x;
    xxdot(8) = - dj3 * lambda1 - j3 * lambda2 + u;

    % Parameter updates
    xxdot(9) = gamma * e * j1;
    xxdot(10) = gamma * e * j2;
    xxdot(11) = gamma * e * j3;
end


