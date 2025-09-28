function [root,relErr,iter] = bisect(func,xl,xu,es,maxit,varargin)
%bisect: root location zeros
% [root,relErr,iter] = bisect(func,xl,xu,es,maxit,p1,p2,...)
% input:
%   func = function handle (@)
%   xl, xu = lower and upper guesses
%   es = desired relative error (default = 0.0001%) (input %)
%   maxit = maximum allowable iterations (default 50)
%   p1,p2... = additional parameters used by func
% output:
%   root = real root
%   relErr = approximate relative error (%)
%   iter = number of iterations

if nargin<3, error('at least 3 input arguments required'), end

% See if crosses x axis
test = func(xl,varargin{:})*func(xu,varargin{:}); % Negative if sign change

if test > 0, error('no sign change'), end
if nargin < 4 || isempty(es), es = 0.0001; end % default
if nargin < 5 || isempty(maxit), maxit = 50; end

iter = 0; xr = xl;
while (1)
    xr_old = xr;
    xr = (xl+xu)/2;
    iter = iter + 1;
    if xr ~= 0, relErr = abs((xr - xr_old)/xr)*100; end
    test = func(xl,varargin{:})*func(xr,varargin{:});
    if test < 0
        xu = xr;
    elseif test > 0
        xl = xr;
    else
        relErr = 0;
    end
    if relErr <= es || iter >= maxit, break, end
end 
root = xr;

end