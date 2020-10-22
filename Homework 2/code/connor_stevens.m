function [V, x] = connor_stevens(t, I_ext)
% [V,X] = CONNOR_STEVENS(T, I_EXT) simulates connor_stevens neuron model
% the state variables are x = [V, n, m, h, a, b]
% the function returns
% 1. V: Voltage trace
% 2. x: state variables across time
%

% Assumes that time is given in units of seconds, and converted into milisecond
t = 1000*t;
dt = t(2) - t(1);
x = zeros(6, numel(t)); % V, n, m, h, a, b

% initial conditions
x(:, 1) = [-67.97; 0.1559; 0.01008; 0.9659; 0.5404; 0.2885];

for tt = 2:length(t)
    dx = connor_stevens_ode(x(:, tt-1), I_ext(tt));
    new_x = x(:, tt-1) + dt * dx;
    
    % the voltage should be bounded between [-80, 50] mV
    new_V = new_x(1);
    new_V(new_V>50) = 50;
    new_V(new_V<-80) = -80;
    new_x = new_x(2:end);
    
    % the state variables should be bounded between [0, 1]
    new_x(new_x<0) = 0;
    new_x(new_x>1) = 1;

    x(:, tt) = [new_V; new_x];
end
V = x(1, :);
x = x(2:end, :);
end