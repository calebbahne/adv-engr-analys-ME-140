%% In-Class 5
% Employ the Newton-Raphson method to determine a real root for the 
% polynomial below. Use an initial guess of 4.5. 
% Determine all the roots of the polynomial with MATLAB ( i.e. use roots())
%f(x) = -2+6*x-4*x^2+0.5*x^3

f = @(x) -2+6*x-4*x^2+0.5*x^3;
fplot(f,[0 10]);

%help newtraph % get the format we need

func = f; xr = 4.5;
dfunc = @(x) 6-8*x+1.5*x^2;
[root_NR,ea_NR,iter_NR]=newtraph(func,dfunc,xr)

% Custom Newton_Raphson method that internally 
% calculates the derivative
[root_NRc,ea_NRc,iter_NRc]=newtraph_autodiff(func,xr)

% Calculate all roots with MATLAB
r = roots([0.5 -4 6 -2])

% Note: there's a function I found that'll automatically
% extract coefficients
syms x;
coeff = sym2poly(f(x)); % Gets in vector form
r_coeff = roots(coeff)
