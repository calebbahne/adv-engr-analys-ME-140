%% In-Class 9 - Integration
% Evaluate using Romberg Integration. Use a stopping criterion of 0.5%. 
% Evaluate using the MATLAB Integral function. 
clear; clc;

a = 1; b = 2; es = 0.5;
f = @(x) (x + 1./x).^(2); % need dot operators for the integral ftn
% help romberg
q_rom = romberg(f,a,b,es)

% Matlab integral function
% help integral
q_integ = integral(f,a,b)