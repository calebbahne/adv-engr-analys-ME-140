function [q,ea,iter] = romberg(func,a,b,es,maxit,varargin)
% romberg: Romberg integration quadrature
%   [q,ea,iter] = romberg(func,a,b,es,maxit)
%   Romberg integration.
% input:
%   func = name of function to be integrated
%   a, b = integration limits
%   es = desired relative error (default = 0.000001%)
%   maxit = maximum allowable iterations (default = 30)
%   p1,p2,... = additional parameters used by func
% output:
%   q = integral estimate
%   ea = approximate relative error (%)
%   iter = number of iterations

if nargin<3,error('at least 3 input arguments required'),end
if nargin<4|isempty(es), es=0.000001;end
if nargin<5|isempty(maxit), maxit=50;end
n = 1;
I(1,1) = trap(func,a,b,n,varargin{:});
iter = 0;
while iter<maxit
iter = iter+1;
n = 2^iter;
I(iter+1,1) = trap(func,a,b,n,varargin{:});
for k = 2:iter+1
j = 2+iter-k;
I(j,k) = (4^(k-1)*I(j+1,k-1)-I(j,k-1))/(4^(k-1)-1);
end
ea = abs((I(1,iter+1)-I(2,iter))/I(1,iter+1))*100;
if ea<=es, break; end
end
q = I(1,iter+1);

%{
ROMBERG INTEGRATION – METHOD OVERVIEW, PERFORMANCE, AND LIMITATIONS
-------------------------------------------------------------------
This function implements Romberg integration, which is a systematic 
refinement of the composite trapezoidal rule using Richardson 
extrapolation. 
  • Subsequent columns apply Richardson extrapolation to eliminate leading
    error terms and accelerate convergence toward the true integral.

At iteration n, the trapezoidal rule is recomputed using 2^n subintervals. 
These values populate the leftmost column of the Romberg table. Then,
extrapolated values are computed using:

      I(j,k) = (4^(k-1)*I(j+1,k-1) - I(j,k-1)) / (4^(k-1) - 1)

This formula cancels higher-order truncation errors, effectively upgrading 
the method from O(h^2) (trapezoid) toward O(h^4), O(h^6), etc., with each 
extrapolation step. Convergence is monitored with the approximate relative 
error between successive refined estimates.

WHEN ROMBERG WORKS WELL
------------------------
Romberg integration is extremely effective when:
  • The integrand is **smooth and sufficiently differentiable** on [a,b].
  • The function has **no sharp corners, discontinuities, or noise**.
  • The interval is finite and well-behaved.
Under these conditions, Richardson extrapolation works as intended, and 
convergence is very rapid—often geometric—yielding high accuracy with fewer 
function evaluations compared to basic composite methods.

WHEN ROMBERG PERFORMS POORLY OR FAILS
--------------------------------------
Romberg integration is less effective or inappropriate when:
  • The integrand has **discontinuities**, **cusps**, or **non-smooth** behavior.
  • The integral is **improper** (infinite limits or singular behavior at endpoints).
  • The function is **highly oscillatory** (e.g., sin(1/x) or high-frequency waves).
  • The evaluation cost of the function is extremely high, since Romberg repeatedly
    recomputes trapezoidal rules with increasing density.
  • Round-off errors accumulate for very high iteration counts, especially if 
    the function varies only slightly across the interval.

In such cases, adaptive quadrature or Gaussian quadrature often performs better.

Overall, Romberg integration is a powerful, general-purpose technique for smooth 
functions, providing high accuracy with structured error control and predictable 
convergence behavior.
%}
