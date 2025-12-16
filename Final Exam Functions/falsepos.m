function [root,froot,ea,iter] = falsepos(func,xl,xu,es,maxit,varargin)
% falsepos: root location zeroes
% [root,froot,ea,iter]=falsepos(func,xl,xu,es,maxit,p1,p2,...):
% uses false position to find the root of func
% input:
%   func = name of function
%   xl, xu = lower and upper guesses
%   es = desired relative error (default = 0.0001%)
%   maxit = maximum allowable iterations (default = 50)
%   p1,p2,... = additional parameters used by function
% output:
%   root = root estimate
%   froot = function at root estimate
%   ea = approximate relative error (%)
%   iter = number of iterations

if nargin < 3
    error('At least 3 input arguments required');
end
test = func(xl,varargin{:}) * func(xu,varargin{:});
if test > 0
    error('No sign change within bracket');
end
if nargin < 4 || isempty(es),    es    = 0.0001; end   % desired error (%)
if nargin < 5 || isempty(maxit), maxit = 50;     end

iter = 0;
ea   = 100;
xr   = xl;

while true
    xrold = xr;
    fl = func(xl,varargin{:});
    fu = func(xu,varargin{:});
    xr = xu - fu*(xl - xu)/(fl - fu);   % false-position formula
    iter = iter + 1;

    if xr ~= 0
        ea = abs((xr - xrold)/xr) * 100;
    end

    if fl * func(xr,varargin{:}) < 0
        xu = xr;
    else
        xl = xr;
    end

    if ea <= es || iter >= maxit
        break
    end
end

root  = xr;
froot = func(xr,varargin{:});
end
