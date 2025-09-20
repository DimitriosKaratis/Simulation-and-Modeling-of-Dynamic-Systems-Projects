%% SIMULATION AND MODELING OF DYNAMIC SYSTEMS
% Assignment 1 - March 2025
% KARATIS DIMITRIOS 10775

%% Exercise 3 - b

clear all;
close all;
clc;


% Define the parameters of the system
m = 0.75;  % kg
L = 1.25; % m
c = 0.15; % N*m*sec
g = 9.81; % m/sec^2
A0 = 4;   % Ν*m
omega = 2; % rad/sec

%% Simulation of the REAL system
% Define the simulation time
T_sim = 20; 
dt = 1e-4;  
tspan = 0:dt:T_sim;

% Initial conditions
x0 = [0; 0];

% Input function
u = @(t) A0 * sin(omega * t);

% Define the system
dynamics = @(t, x) [x(2); (1/(m*L^2)) * (u(t) - c*x(2) - m*g*L*x(1))];

% Solve the system using ode45
[t, X] = ode45(dynamics, tspan, x0);

q = X(:,1);
q_dot = X(:,2);

%% Estimating with Least-Squares method

% Poles of the system
p1 = 0.5;
p2 = 0.5;

%% Initialize vector for dt_sampled values

dt_sampled_values = [0.001, 0.01, 0.1];  
error_values_L = zeros(1, length(dt_sampled_values));
error_values_m = zeros(1, length(dt_sampled_values));
error_values_c = zeros(1, length(dt_sampled_values));

for i = 1:length(dt_sampled_values)
    dt_sampled = dt_sampled_values(i); 
    t_sampled = 0:dt_sampled:T_sim;    

    % Interpolate the high-resolution data to match the new sampled time
    q_sampled = interp1(t, q, t_sampled);   
    u_sampled = A0 * sin(omega * t_sampled); 

    % Create the phi matrix 
    phi = zeros(length(t_sampled'), 3);

    phi1 = lsim(tf([-1 0], [1 (p1+p2) p1*p2]), q_sampled', t_sampled');
    phi2 = lsim(tf(-1, [1 (p1+p2) p1*p2]), q_sampled', t_sampled');
    phi3 = lsim(tf(1, [1 (p1+p2) p1*p2]), u_sampled', t_sampled');

    phi(:,1) = phi1;
    phi(:,2) = phi2;
    phi(:,3) = phi3;

    % Phi squared (or ΦΤΦ)
    phiTphi = (phi') * phi;                              

    theta0 = inv(phiTphi) * (phi') * (q_sampled');

    % Find the estimated parameters
    L_estim = g / (theta0(2) + (p1 * p2));
    m_estim = 1 / ((L_estim^2) * (theta0(3)));
    c_estim = (m_estim * L_estim^2) * (theta0(1) + (p1 + p2));

    % Calculate errors
    error_values_L(i) = abs(L - L_estim);
    error_values_m(i) = abs(m - m_estim);
    error_values_c(i) = abs(c - c_estim);
    
end

%% Plot the errors vs dt_sampled_values
figure;
plot(dt_sampled_values, error_values_L, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 8);
hold on;
plot(dt_sampled_values, error_values_m, 'rs-', 'LineWidth', 1.5, 'MarkerSize', 8);
plot(dt_sampled_values, error_values_c, 'gd-', 'LineWidth', 1.5, 'MarkerSize', 8);
hold off;
grid on;
xlabel('Sampling Time (dt) [s]');
ylabel('Absolute Error');
title('Parameter Estimation Errors vs Sampling Time');
legend('L error', 'm error', 'c error', 'Location', 'best');
