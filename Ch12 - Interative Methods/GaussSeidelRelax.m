function [x,ea,iter] = GaussSeidelRelax(A,b,lambda,es,maxit)
% GaussSeidel: Gauss Seidel method with relaxation
% [x,ea,iter] = GaussSeidelRelax(A,b,lambda,es,maxit): Gauss Seidel with relaxation
% input:
% A = coefficient matrix
% b = right hand side vector
% lambda = relation factor (default = 1)
% es = stop criterion (default = 0.00001%)
% maxit = max iterations (default = 50)
% output:
% x = solution vector
% ea = maximum relative error (%)
% iter = number of iterations

% GaussSeidelRelax: Solves a system of linear equations [A]{x}={b} using the
%                   Gauss-Seidel method enhanced with a relaxation factor (lambda).
%
% This function is a generalization of the standard Gauss-Seidel method. The
% relaxation factor allows the user to intentionally over- or under-estimate
% the new value of x(i) to control convergence rate and stability.
%
% [x,ea,iter] = GaussSeidelRelax(A,b,lambda,es,maxit):
%
% input:
% A      = coefficient matrix (must be square)
% b      = right hand side vector
% lambda = **Relaxation factor.** Must be 0 < lambda <= 2.
%          - lambda = 1: Standard Gauss-Seidel (no relaxation).
%          - 0 < lambda < 1 (Underrelaxation): Used for non-convergent or unstable systems to dampen oscillations and improve stability.
%          - 1 < lambda <= 2 (Overrelaxation/Successive Over-Relaxation - SOR): Used to accelerate convergence for systems that already converge.
% es     = stop criterion (desired approximate percent relative error, default = 0.00001%)
% maxit  = maximum allowable iterations (default = 50)
%
% output:
% x      = solution vector
% ea     = final maximum approximate percent relative error (%) across all elements of x
% iter   = number of iterations performed
%
% -------------------------------------------------------------------------
% CONVERGENCE, USAGE, AND EFFECTS OF LAMBDA
% -------------------------------------------------------------------------
% 1. Standard Gauss-Seidel (lambda=1):
%    - Works best/guaranteed to converge for Diagonally Dominant systems.
%    - May diverge for non-diagonally dominant systems.
%
% 2. Relaxation (lambda != 1):
%    - **Accelerated Convergence (Overrelaxation, 1 < lambda <= 2):** When a system
%      already converges, choosing an optimal lambda > 1 (usually between 1.2 and 1.8)
%      can drastically reduce the required number of iterations, leading to faster
%      computation.
%    - **Stabilization (Underrelaxation, 0 < lambda < 1):** For highly unstable or
%      oscillatory systems (which might diverge at lambda=1), using lambda < 1 can
%      dampen the updates, forcing the method toward a stable solution at the cost
%      of slower convergence.
%    - **When it FAILS:** If the optimal lambda is chosen poorly, or if the system is
%      fundamentally ill-conditioned, the method may still diverge or converge
%      very slowly, regardless of the relaxation factor.
% -------------------------------------------------------------------------

if nargin<2,error('at least 2 input arguments required'),end
if nargin<5|isempty(maxit),maxit=50;end
if nargin<4|isempty(es),es=0.00001;end
if nargin<3|isempty(lambda),lambda=1;end
[m,n] = size(A);
if m~=n, error('Matrix A must be square'); end
C = A;
for i = 1:n
C(i,i) = 0;
x(i) = 0;
end
x = x';
for i = 1:n
C(i,1:n) = C(i,1:n)/A(i,i);
end
for i = 1:n
d(i) = b(i)/A(i,i);
end
iter = 0;
while (1)
xold = x;
for i = 1:n
x(i) = d(i)-C(i,:)*x;
x(i) = lambda*x(i) + (1 - lambda)*xold(i);
if x(i) ~= 0
ea(i) = abs((x(i) - xold(i))/x(i)) * 100;
end
end
iter = iter+1;
if max(ea)<=es | iter >= maxit, break, end
end
ea=max(ea);
