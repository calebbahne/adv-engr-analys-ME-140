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
%{
    Goal: solve f(x) = 0 on [a,b] using fixed point iteration.
    Transform (algebraically) f(x) = 0 to g(x) = x, then find where the fixed point is. 
    Requirements for g(x):
        If all of g(x) on [a,b] has a range within [a,b], there's at least one fixed point on the interval.
        If (1) and g'(x) (the slope) is always within [-1,1], then the fixed point is unique.
        We want to find situations where both 1 and 2 apply before we use g(x).
    Input: g, p0, tol, Maxit
        g(x) is input using the @ symbol (e.g. @(x) x.^2 + 3*x + 2)
        p0 is the initial guess for the solution
        tol is the minimum tolerance required for the final solution (i.e. 10^-4)
        Maxit is the maximum number of iterations
    Output: p, i
        p is the solution, within tol of the true fixed point
        i is the number of iterations required to get there.
%}

if nargin < 2, error('at least 2 arguments required'), end
if nargin < 3||isempty(es),es = 1e-6;end %if es is blank se
if nargin < 4||isempty(maxit),maxit = 50;end %if maxit is
iter = 0; ea = 100;
while (1)
x1 = g(x0);
iter = iter + 1;
if x1 ~= 0, ea = abs((x1 - x0)/x1)*100; end
if (ea <= es || iter >= maxit),break,end
x0=x1;
end
end
