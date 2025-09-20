%% SIMULATION AND MODELING OF DYNAMIC SYSTEMS
% Assignment 2 - April 2025
% KARATIS DIMITRIOS 10775

%% Exercise 1 - b
% Lyapunov Method / Series-Parallel Configuration Estimator
function xxdot = LyapSP(t,xx)

    global A B c d f gamma thetam
    
    % xx description
    x1 = xx(1); 
    x2 = xx(2); 
    x1_est = xx(3);
    x2_est = xx(4);
    a21_est = xx(5);
    a22_est = xx(6);
    b2_est = xx(7);

    % Input
    u = c * sin(d * t) + f;

    % Errors
    e1 = x1 - x1_est;
    e2 = x2 - x2_est;

    xxdot = zeros(size(xx));
    
    xxdot(1) = A(1,1)*x1 + A(1,2)*x2 + B(1)*u;                             % x1_dot
    xxdot(2) = A(2,1)*x1 + A(2,2)*x2 + B(2)*u;                             % x2_dot
    xxdot(3) = x2 + thetam(1) * e1 + thetam(2) * e2;                       % x1_hat_dot
    xxdot(4) = a21_est*x1 + a22_est*x2 + b2_est*u + thetam(3) * e1 + thetam(4) * e2; % x2_hat_dot
    xxdot(5) = gamma(1) * e2 * x1;                                         % a21_hat_dot
    xxdot(6) = gamma(1) * e2 * x2;                                         % a22_hat_dot
    xxdot(7) = gamma(2) * e2 * u;                                          % b2_hat_dot

end

