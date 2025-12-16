function [d,ea,iter]=rombdiff(func,x,es,maxit,varargin)
% romberg: Romberg integration quadrature
%   [d,ea,iter]=rombdiff(func,x,es,maxit,varargin)
% Romberg differentiation. 
%   Implicitly calls dydxnew
% input:
%   func = name of function to be integrated
%   es = desired relative error (default = 0.000001%)
%   maxit = maximum allowable iterations (default = 50)
%   p1,p2,... = additional parameters used by func
% output:
%   d = integral estimate
%   ea = approximate relative error (%)
%   iter = number of iterations
if nargin<2,error('at least 2 input arguments required'),end
if nargin<3|isempty(es), es=0.000001;end
if nargin<4|isempty(maxit), maxit=50;end
n = 1;
DY(1,1) = dydxnew(func,x,n,varargin{:});
iter = 0;
ea=100;
while iter<maxit
iter = iter+1;
n = 2^iter;
DY(iter+1,1) = dydxnew(func,x,n,varargin{:});
for k = 2:iter+1
j = 2+iter-k;
DY(j,k) = (4^(k-1)*DY(j+1,k-1)-DY(j,k-1))/(4^(k-1)-1);
end
if DY(1,iter+1)~=0,ea = abs((DY(1,iter+1)-DY(2,iter))/DY(1,iter+1))*100;end
if ea<=es, break; end
end
d = DY(1,iter+1);

%{
ROMBERG DIFFERENTIATION – OVERVIEW AND LIMITATIONS
--------------------------------------------------
This function applies Romberg-style Richardson extrapolation to numerical
differentiation. It repeatedly computes symmetric finite-difference
approximations with shrinking step sizes (h = x/n) using:

    d ≈ [f(x+h) – f(x–h)] / (2h)

These base estimates fill the first column of the Romberg tableau, and
higher-order columns cancel leading truncation errors to improve accuracy.
The method stops when successive extrapolated values converge within the
desired relative error.

WHEN THIS METHOD IS LESS EFFECTIVE OR INAPPROPRIATE
----------------------------------------------------
Romberg differentiation can struggle or fail in the following cases:

• **Non-smooth functions near x**  
  Kinks, cusps, discontinuities, or nondifferentiability cause large
  errors because Richardson extrapolation assumes a smooth Taylor series.

• **Functions with significant numerical noise**  
  Subtracting values very close together (as h → 0) amplifies noise and
  floating-point roundoff, often causing divergence rather than refinement.

• **Very steep or oscillatory behavior near x**  
  Rapid changes make finite differences unstable and can destroy the
  expected error cancellation.

• **x near a domain boundary**  
  Centered differences cannot be used cleanly if x ± h leaves the domain.

• **Expensive function evaluations**  
  Romberg repeatedly evaluates f(x±h) at many shrinking step sizes,
  making it costly compared to simpler derivative approximations.

In summary, Romberg differentiation works best for smooth, well-behaved
functions near the point of interest, but becomes unreliable when the
function is noisy, non-smooth, oscillatory, or poorly scaled.
%}
