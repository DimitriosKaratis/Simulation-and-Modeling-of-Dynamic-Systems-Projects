%% SIMULATION AND MODELING OF DYNAMIC SYSTEMS
% Assignment 2 - April 2025
% KARATIS DIMITRIOS 10775

%% Exercise 1 - c
%% Real time estimator / Lyapunov Method
                           % i) Parallel Configuration with Noise
                           % ii) Series-Parallel Configuration with Noise

%% i) Lyapunov Method / Real time estimator / Parallel Configuration with Noise
%% Clearing
clear all;
close all;
clc;
format long

tic;                        % Start clock for code evaluation

%% Real values of parameters and configuration parameters
global A B c d f gamma h0 f0

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
h0 = 0.25;
f0 = 20;
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
[t,x] = ode45(@(tpar,xpar) LyapParNoise(tpar,xpar), tspan, initCond);

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
title('x_1 & x_1_,_e_s_t (Parallel with NOISE)')
xlabel('Time [s]')
ylabel('x_1 ,x_1_,_e_s_t')
grid on
legend('x_1','x_1_,_e_s_t')


%% Error Plot: e_x(t) = x(t) - x̂(t)

figure
plot(t,e1);
title('Error of x_1 estimation (e_1) (Parallel with NOISE)')
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
title('Estimated Parameters Over Time (Parallel with NOISE)');
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
fprintf('\n=== Final Parameter Estimates (Parallel Configuration with NOISE) ===\n');
fprintf('Mass:      m_true = %.4f\t m_est = %.4f\n', m, m_hat(end));
fprintf('Damping:   b_true = %.4f\t b_est = %.4f\n', b, b_hat(end));
fprintf('Stiffness: k_true = %.4f\t k_est = %.4f\n', k, k_hat(end));

toc;                                      

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ii) Lyapunov Method / Real time estimator / Series-Parallel Configuration with Noise
global thetam

thetam = [0.2 0.2 0.2 0.2];
gamma = [0.06 1.5]; 

tic;                        % Start clock for code evaluation

% Time Span
tStart = 0;   
tStep = 0.01;
tEnd = 100;
tspan = tStart:tStep:tEnd;
% Initial conditions: [x1; x2; x1_est; x2_est; a21_est; a22_est; b2_est]
initCond = [0; 0; 0; 0; 0.01; 0.01; 0.01];                                 
% Solving ODE 
[t,x] = ode45(@(tpar,xpar) LyapSPNoise(tpar,xpar), tspan, initCond);

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
title('x_1 & x_1_,_e_s_t (Series-Parallel with NOISE)')
xlabel('Time [s]')
ylabel('x_1 ,x_1_,_e_s_t')
grid on
legend('x_1','x_1_,_e_s_t')

%% Error Plot: e_x(t) = x(t) - x̂(t)

figure
plot(t,e1);
title('Error of x_1 estimation (e_1) (Series-Parallel with NOISE)')
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
m_hat_smooth = movmean(m_hat_plot, windowSize);
b_hat_smooth = movmean(b_hat_plot, windowSize);
k_hat_smooth = movmean(k_hat_plot, windowSize);

figure;
subplot(3,1,1);
plot(t_plot, m_hat_smooth, 'LineWidth', 1.5); 
yline(m, '--r', 'LineWidth', 1.5); 
ylabel('m̂(t)'); 
grid on; 
title('Estimated Parameters Over Time (Series-Parallel with NOISE)');
legend('m̂(t)', 'm_{true}');

subplot(3,1,2);
plot(t_plot, b_hat_smooth, 'LineWidth', 1.5); 
yline(b, '--r', 'LineWidth', 1.5); 
ylabel('b̂(t)'); 
grid on;
legend('b̂(t)', 'b_{true}');

subplot(3,1,3);
plot(t_plot, k_hat_smooth, 'LineWidth', 1.5); 
yline(k, '--r', 'LineWidth', 1.5); 
ylabel('k̂(t)'); 
xlabel('Time [s]');
grid on;
legend('k̂(t)', 'k_{true}');

%% Final Estimated Values
fprintf('\n=== Final Parameter Estimates (Series-Parallel Configuration with NOISE) ===\n');
fprintf('Mass:      m_true = %.4f\t m_est = %.4f\n', m, m_hat(end));
fprintf('Damping:   b_true = %.4f\t b_est = %.4f\n', b, b_hat(end));
fprintf('Stiffness: k_true = %.4f\t k_est = %.4f\n', k, k_hat(end));

toc;                                      

%% Effect of Noise Amplitude on Parameter Estimation Accuracy (Parallel Configuration)

gamma = [0.06 0.019];

noiseLevels = 0:0.05:0.5;  % Noise amplitude range
numLevels = length(noiseLevels);

% Store estimation errors
m_errors = zeros(1, numLevels);
b_errors = zeros(1, numLevels);
k_errors = zeros(1, numLevels);

% Loop over different noise amplitudes
for i = 1:numLevels
    h0 = noiseLevels(i);  
    
    % Initial conditions
    initCond = [0; 0; 0; 0; 0; 0; 0];                                 
    [t,x] = ode45(@(tpar,xpar) LyapParNoise(tpar,xpar), tspan, initCond);

    a21_est = x(:,5);
    a22_est = x(:,6);
    b2_est  = x(:,7);

    m_hat = 1 ./ b2_est;
    b_hat = -a22_est .* m_hat;
    k_hat = -a21_est .* m_hat;

    % Final estimation error
    m_errors(i) = abs(m_hat(end) - m);
    b_errors(i) = abs(b_hat(end) - b);
    k_errors(i) = abs(k_hat(end) - k);
end

% Plotting estimation errors vs noise amplitude
figure;
subplot(3,1,1);
plot(noiseLevels, m_errors, '-o', 'LineWidth', 1.5);
ylabel('|m̂ - m|');
title('Effect of Noise Amplitude on Parameter Estimation (Parallel)');
grid on;

subplot(3,1,2);
plot(noiseLevels, b_errors, '-s', 'LineWidth', 1.5);
ylabel('|b̂ - b|');
grid on;

subplot(3,1,3);
plot(noiseLevels, k_errors, '-^', 'LineWidth', 1.5);
ylabel('|k̂ - k|');
xlabel('Noise Amplitude h_0');
grid on;

%% Effect of Noise Amplitude on Parameter Estimation Accuracy (Series-Parallel Configuration)

thetam = [0.2 0.2 0.2 0.2];
gamma = [0.06 1.5];

noiseLevels = 0:0.05:0.5;  % Noise amplitude range
numLevels = length(noiseLevels);

% Store estimation errors
m_errors = zeros(1, numLevels);
b_errors = zeros(1, numLevels);
k_errors = zeros(1, numLevels);

% Loop over different noise amplitudes
for i = 1:numLevels
    h0 = noiseLevels(i);  
    
    % Initial conditions
    initCond = [0; 0; 0; 0; 0.01; 0.01; 0.01];                                 
    [t,x] = ode45(@(tpar,xpar) LyapSPNoise(tpar,xpar), tspan, initCond);

    a21_est = x(:,5);
    a22_est = x(:,6);
    b2_est  = x(:,7);

    m_hat = 1 ./ b2_est;  
    b_hat = -a22_est .* m_hat;
    k_hat = -a21_est .* m_hat;

    % Final estimation error
    m_errors(i) = abs(m_hat(end) - m);
    b_errors(i) = abs(b_hat(end) - b);
    k_errors(i) = abs(k_hat(end) - k);
end

% Plotting estimation errors vs noise amplitude
figure;
subplot(3,1,1);
plot(noiseLevels, m_errors, '-o', 'LineWidth', 1.5);
ylabel('|m̂ - m|');
title('Effect of Noise Amplitude on Parameter Estimation (Series-Parallel)');
grid on;

subplot(3,1,2);
plot(noiseLevels, b_errors, '-s', 'LineWidth', 1.5);
ylabel('|b̂ - b|');
grid on;

subplot(3,1,3);
plot(noiseLevels, k_errors, '-^', 'LineWidth', 1.5);
ylabel('|k̂ - k|');
xlabel('Noise Amplitude h_0');
grid on;