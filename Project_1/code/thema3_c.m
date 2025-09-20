%% SIMULATION AND MODELING OF DYNAMIC SYSTEMS
% Assignment 1 - March 2025
% KARATIS DIMITRIOS 10775

%% Exercise 3 - c 

clear all;
close all;
clc;

% Define the parameters of the system
m = 0.75;   % kg
L = 1.25;   % m
c = 0.15;   % N*m*sec
g = 9.81;   % m/sec^2
omega = 2;  % rad/sec

% Define the simulation time and fixed sampling time
T_sim = 20; 
dt = 1e-4;  
tspan = 0:dt:T_sim;

% Initial conditions
x0 = [0; 0];

% Define A0 values to test
A0_values = [0.01, 0.1, 1, 10, 100]; 

% Preallocate error storage
error_values_L = zeros(1, length(A0_values));
error_values_m = zeros(1, length(A0_values));
error_values_c = zeros(1, length(A0_values));

% Poles of the system 
p1 = 0.5;
p2 = 0.5;

% Loop over each A0 value
for i = 1:length(A0_values)
    A0 = A0_values(i);
    
    % Input function
    u = @(t) A0 * sin(omega * t);
    
    % Simulate the real system
    dynamics = @(t, x) [x(2); (1/(m*L^2)) * (u(t) - c*x(2) - m*g*L*x(1))];
    [t, X] = ode45(dynamics, tspan, x0);
    q = X(:,1);
    q_dot = X(:,2);
    
    % Downsampling to 0.1-second intervals
    dt_sampled = 0.1;  
    t_sampled = 0:dt_sampled:T_sim;  
    
    % Interpolate the high-resolution data to match the new sampling points
    q_sampled = interp1(t, q, t_sampled);   
    u_sampled = A0 * sin(omega * t_sampled);  
        
    % Least-squares estimation
    phi = zeros(length(t_sampled'), 3);
    phi(:,1) = lsim(tf([-1 0], [1 (p1+p2) p1*p2]), q_sampled', t_sampled');
    phi(:,2) = lsim(tf(-1, [1 (p1+p2) p1*p2]), q_sampled', t_sampled');
    phi(:,3) = lsim(tf(1, [1 (p1+p2) p1*p2]), u_sampled', t_sampled');
    
    % Phi squared (or ΦΤΦ)
    phiTphi = (phi') * phi;                              
    
    theta0 = inv(phiTphi) * (phi') * (q_sampled');
    
    % Estimated parameters
    L_estim = g / (theta0(2) + (p1 * p2));
    m_estim = 1 / ((L_estim^2) * theta0(3));
    c_estim = (m_estim * L_estim^2) * (theta0(1) + (p1 + p2));
    
    % Calculate errors
    error_values_L(i) = abs(L - L_estim);
    error_values_m(i) = abs(m - m_estim);
    error_values_c(i) = abs(c - c_estim);
end

%% Plot the errors vs A0 values
figure;
plot(A0_values, error_values_L, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 8);
hold on;
plot(A0_values, error_values_m, 'rs-', 'LineWidth', 1.5, 'MarkerSize', 8);
plot(A0_values, error_values_c, 'gd-', 'LineWidth', 1.5, 'MarkerSize', 8);
hold off;
grid on;
xlabel('Input Amplitude (A0) [N*m]');
ylabel('Absolute Error');
title('Parameter Estimation Errors vs Input Amplitude');
legend('L error', 'm error', 'c error', 'Location', 'best');

%% Plot the error for L vs A0 values
figure;
plot(A0_values, error_values_L, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 8);
grid on;
xlabel('Input Amplitude (A0) [N*m]');
ylabel('Absolute Error in L');
title('Error in L vs Input Amplitude');

%% Plot the error for m vs A0 values
figure;
plot(A0_values, error_values_m, 'rs-', 'LineWidth', 1.5, 'MarkerSize', 8);
grid on;
xlabel('Input Amplitude (A0) [N*m]');
ylabel('Absolute Error in m');
title('Error in m vs Input Amplitude');

%% Plot the error for c vs A0 values
figure;
plot(A0_values, error_values_c, 'gd-', 'LineWidth', 1.5, 'MarkerSize', 8);
grid on;
xlabel('Input Amplitude (A0) [N*m]');
ylabel('Absolute Error in c');
title('Error in c vs Input Amplitude');
