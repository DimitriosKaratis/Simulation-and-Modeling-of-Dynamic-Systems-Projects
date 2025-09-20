%% SIMULATION AND MODELING OF DYNAMIC SYSTEMS
% Assignment 2 - April 2025
% KARATIS DIMITRIOS 10775

%% Exercise 1 - b
%% Real time estimator / Lyapunov Method
                           % i) Parallel Configuration
                           % ii) Series-Parallel Configuration

%% i) Lyapunov Method / Real time estimator / Parallel Configuration
%% Clearing
clear all;
close all;
clc;
format long

%% Real values of parameters and configuration parameters
global A B c d f  gamma

% True system parameters
m = 1.315;         % Mass
b = 0.225;         % Damping
k = 0.725;         % Stiffness

% A Matrix real parameters
a11 = 0;
a12 = 1;
a21 = -k/m;
a22 = -b/m;
A = [a11 a12; a21 a22];


% b vector real parameters
b1 = 0;
b2 = 1/m;
B = [b1;b2];

% Input parameters
c = 2.5;
d = 1;
f = 0;
gamma1 = 0.06;
gamma2 = 0.019;
gamma = [gamma1 gamma2]; 

% Time Span
tStart = 0;   
tStep = 0.01;
tEnd = 100;
tspan = tStart:tStep:tEnd;
% Initial conditions: [x1; x2; x1_est; x2_est; a21_est; a22_est; b2_est]
initCond = [0; 0; 0; 0; 0.0; 0.0; 0.0];
% Solving ODE 
[t,x] = ode45(@(tpar,xpar) LyapPar(tpar,xpar), tspan, initCond);

% Getting results 
x1 = x(:,1);
x2 = x(:,2);
x1_est = x(:,3);
x2_est = x(:,4);
a21_est = x(:,5);
a22_est = x(:,6);
b2_est = x(:,7);

% Errors
e1 = x1 - x1_est;
e2 = x2 - x2_est;

%% Comparison Plot: x(t) vs x̂(t)
figure
plot(t,x1,'r',t,x1_est,'b');
title('x_1 & x_1_,_e_s_t (Parallel)')
xlabel('Time [s]')
ylabel('x_1 ,x_1_,_e_s_t')
grid on
legend('x_1','x_1_,_e_s_t')

%% Error Plot: e_x(t) = x(t) - x̂(t)

figure
plot(t,e1);
title('Error of x_1 estimation (e_1) (Parallel)')
xlabel('Time [s]')
ylabel('e_1')
grid on

%% Parameter Estimation Plot: m̂(t), b̂(t), k̂(t)

m_hat = 1 ./ b2_est ;  
b_hat = -a22_est .* m_hat;
k_hat = -a21_est .* m_hat;

% Mask for t >= 0.5
mask = t >= 0.5;

% Apply mask to all relevant vectors
t_plot = t(mask);
m_hat_plot = m_hat(mask);
b_hat_plot = b_hat(mask);
k_hat_plot = k_hat(mask);

figure;
subplot(3,1,1);
plot(t_plot, m_hat_plot, 'LineWidth', 1.5); 
yline(m, '--r', 'LineWidth', 1.5); 
ylabel('m̂(t)'); 
grid on; 
title('Estimated Parameters Over Time (Parallel)');
legend('m̂(t)', 'm_{true}');

subplot(3,1,2);
plot(t_plot, b_hat_plot, 'LineWidth', 1.5); 
yline(b, '--r', 'LineWidth', 1.5); 
ylabel('b̂(t)'); 
grid on;
legend('b̂(t)', 'b_{true}');

subplot(3,1,3);
plot(t_plot, k_hat_plot, 'LineWidth', 1.5); 
yline(k, '--r', 'LineWidth', 1.5); 
ylabel('k̂(t)'); 
xlabel('Time [s]');
grid on;
legend('k̂(t)', 'k_{true}');

%% Final Estimated Values
fprintf('\n=== Final Parameter Estimates (Parallel Configuration) ===\n');
fprintf('Mass:      m_true = %.4f\t m_est = %.4f\n', m, m_hat(end));
fprintf('Damping:   b_true = %.4f\t b_est = %.4f\n', b, b_hat(end));
fprintf('Stiffness: k_true = %.4f\t k_est = %.4f\n', k, k_hat(end));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ii) Lyapunov Method / Real time estimator / Series-Parallel Configuration
global thetam

thetam = [0.2 0.2 0.2 0.2];
gamma = [0.06 1.5]; 

% Time Span
tStart = 0;   
tStep = 0.01;
tEnd = 100;
tspan = tStart:tStep:tEnd;
% Initial conditions: [x1; x2; x1_est; x2_est; a21_est; a22_est; b2_est]
initCond = [0; 0; 0; 0; 0.01; 0.01; 0.01];                                 
% Solving ODE 
[t,x] = ode45(@(tpar,xpar) LyapSP(tpar,xpar), tspan, initCond);

% Getting results 
x1 = x(:,1);
x2 = x(:,2);
x1_est = x(:,3);
x2_est = x(:,4);
a21_est = x(:,5);
a22_est = x(:,6);
b2_est = x(:,7);

A_est = [a11 a12; a21_est a22_est];
b_est = [b1; b2_est];

% Errors
e1 = x1 - x1_est;
e2 = x2 - x2_est;
e = [e1; e2];

%% Comparison Plot: x(t) vs x̂(t)
figure
plot(t,x1,'r',t,x1_est,'b');
title('x_1 & x_1_,_e_s_t (Series-Parallel)')
xlabel('Time [s]')
ylabel('x_1 ,x_1_,_e_s_t')
grid on
legend('x_1','x_1_,_e_s_t')


%% Error Plot: e_x(t) = x(t) - x̂(t)

figure
plot(t,e1);
title('Error of x_1 estimation (e_1) (Series-Parallel)')
xlabel('Time [s]')
ylabel('e_1')
grid on

%% Parameter Estimation Plot: m̂(t), b̂(t), k̂(t)

m_hat = 1 ./ b2_est ;  
b_hat = -a22_est .* m_hat;
k_hat = -a21_est .* m_hat;

% Mask for t >= 0.5
mask = t >= 0.5;

% Apply mask to all relevant vectors
t_plot = t(mask);
m_hat_plot = m_hat(mask);
b_hat_plot = b_hat(mask);
k_hat_plot = k_hat(mask);

% Apply moving average filter to reject big spikes
windowSize = 500;
m_hat_plot = movmean(m_hat_plot, windowSize);
b_hat_plot = movmean(b_hat_plot, windowSize);
k_hat_plot = movmean(k_hat_plot, windowSize);


figure;
subplot(3,1,1);
plot(t_plot, m_hat_plot, 'LineWidth', 1.5); 
yline(m, '--r', 'LineWidth', 1.5); 
ylabel('m̂(t)'); 
grid on; 
title('Estimated Parameters Over Time (Series-Parallel)');
legend('m̂(t)', 'm_{true}');

subplot(3,1,2);
plot(t_plot, b_hat_plot, 'LineWidth', 1.5); 
yline(b, '--r', 'LineWidth', 1.5); 
ylabel('b̂(t)'); 
grid on;
legend('b̂(t)', 'b_{true}');

subplot(3,1,3);
plot(t_plot, k_hat_plot, 'LineWidth', 1.5); 
yline(k, '--r', 'LineWidth', 1.5); 
ylabel('k̂(t)'); 
xlabel('Time [s]');
grid on;
legend('k̂(t)', 'k_{true}');

%% Final Estimated Values
fprintf('\n=== Final Parameter Estimates (Series-Parallel Configuration) ===\n');
fprintf('Mass:      m_true = %.4f\t m_est = %.4f\n', m, m_hat(end));
fprintf('Damping:   b_true = %.4f\t b_est = %.4f\n', b, b_hat(end));
fprintf('Stiffness: k_true = %.4f\t k_est = %.4f\n', k, k_hat(end));
                                     
