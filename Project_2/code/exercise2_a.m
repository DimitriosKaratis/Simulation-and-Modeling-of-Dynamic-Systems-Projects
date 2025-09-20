%% SIMULATION AND MODELING OF DYNAMIC SYSTEMS
% Assignment 2 - April 2025
% KARATIS DIMITRIOS 10775

%% Exercise 2 - a

clear;
close all;
clc;

% System Parameters
a1 = 1.315;
a2 = 0.725;
a3 = 0.225;
b = 1.175;

% Controller Parameters
phi_0 = 1.0;
phi_inf = 0.01;
lambda = 1;
rho = 1.5;
k1 = 2;
k2 = 4;

% Desired Trajectory
t_final = 20;
r_d_fun = @(t) (pi/20)*(1 - cos(2*pi*t/t_final)); % Desired position: 0 -> pi/10 -> 0
dr_d_fun = @(t) (pi/20)*(2*pi/t_final)*sin(2*pi*t/t_final); % First derivative (velocity)
d2r_d_fun = @(t) (pi/20)*(2*pi/t_final)^2*cos(2*pi*t/t_final); % Second derivative (acceleration)

% Nonlinear transformation T(z)
T = @(z) log((1 + z) ./ (1 - z)); 

% System Dynamics with Nonlinear Controller
system = @(t, x) [
    % dx1/dt = x2 (first state derivative)
    x(2);
    % Second derivative
    -a1*x(2) - a2*sin(x(1)) + a3*x(2)^2*sin(2*x(1)) + ...
    b * u_nonlinear(t, x, r_d_fun, phi_0, phi_inf, lambda, rho, k1, k2, T)  
];

% Initial Conditions and Simulation Time Span
x0 = [0; 0]; 
tspan = [0 t_final];

% Solve ODE using ode45
[t, x] = ode45(@(t,x) system(t, x), tspan, x0);

% Compute desired trajectory for plotting
r_d = arrayfun(r_d_fun, t);

% Plotting Results
figure;
plot(t, x(:,1), 'b', 'LineWidth', 2); hold on;      % Actual response
plot(t, r_d, 'r--', 'LineWidth', 2);                % Desired trajectory
xlabel('Time [sec]');
ylabel('Angle r(t) [rad]');
legend('System Response r(t)', 'Desired Trajectory r_d(t)');
title('System Response with Nonlinear Controller');
grid on;