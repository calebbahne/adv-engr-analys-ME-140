function [root,ea,iter] = newtraph_autodiff(func,xr,es,maxit,varargin)
% newtraph_autodiff: Newton-Raphson root finding with automatic differentiation
% [root,ea,iter] = newtraph_autodiff(func,xr,es,maxit,p1,p2,...)
% 
% Uses symbolic differentiation to compute derivative automatically!
%
% input:
%   func  = function handle (e.g., @(x) x^3 - 6*x^2 + 11*x - 6.1)
%   xr    = initial guess
%   es    = desired relative error (default = 0.0001%)
%   maxit = maximum iterations (default = 50)
%   p1,p2,... = additional parameters for func
%
% output:
%   root  = real root
%   ea    = approximate relative error (%)
%   iter  = number of iterations

if nargin < 2, error('At least 2 input arguments required'), end
if nargin < 3 || isempty(es), es = 0.0001; end
if nargin < 4 || isempty(maxit), maxit = 50; end

% Convert function handle to symbolic
syms x
f_sym = func(x,varargin{:});       % evaluate function handle in symbolic x
df_sym = diff(f_sym, x);           % symbolic derivative

% Convert back to function handles
f = matlabFunction(f_sym, 'Vars', {x,varargin{:}});
df = matlabFunction(df_sym, 'Vars', {x,varargin{:}});

iter = 0; 
ea = 100; 

while true
    xrold = xr;
    xr = xr - f(xr,varargin{:})/df(xr,varargin{:});
    iter = iter + 1;
    if xr ~= 0
        ea = abs((xr - xrold)/xr) * 100;
    end
    if ea <= es || iter >= maxit
        break
    end
end

root = xr;
end
