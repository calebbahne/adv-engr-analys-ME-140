function [x1,ea,iter] = fixpt(g,x0,es,maxit)
% fixpt: fixed point iteration root locator
% [x1,ea,iter] = fixpt(g,x0,es,maxit)
% This function determines the root of x = g(x) with fixed
% The method is repeated until either the percent relative
% is equal to or less than es (default:1.e–6) or the number
% exceeds maxit (default:50).
% Input:
%   g = the function for g(x)
%   x0 = initial guess for x
%   es = relative error stopping criterion (%)
%   maxit = maximum number of iterations
% Output:
%   x1 = solution estimate
%   ea = relative error
%   iter = number of iterations

if nargin < 2, error('at least 2 arguments required'), end
if nargin < 3||isempty(es),es = 1e–6;end %if es is blank se
if nargin < 4||isempty(maxit),maxit = 50;end %if maxit is
iter = 0; ea = 100;
while (1)
x1 = g(x0);
iter = iter + 1;
if x1 ~= 0, ea = abs((x1 − x0)/x1)*100; end
if (ea <= es || iter >= maxit),break,end
x0=x1;
end
end
