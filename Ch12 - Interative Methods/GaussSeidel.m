function [x,ea,iter] = GaussSeidel(A,b,es,maxit)
% GaussSeidel: Gauss Seidel method
% [x,ea,iter] = GaussSeidel(A,b,es,maxit): Gauss Seidel without relaxation
% input:
% A = coefficient matrix
% b = right hand side vector
% es = stop criterion (default = 0.00001%)
% maxit = max iterations (default = 50)
% output:
% x = solution vector
% ea = maximum relative error (%)
% iter = number of iterations

% % GaussSeidel: Solves a system of linear equations [A]{x}={b} using the Gauss-Seidel iterative method.
%
% This function is suitable for solving systems that are large or sparse, where direct methods
% (like matrix inversion or Gaussian elimination) may be less efficient.
% The method iteratively refines the solution vector 'x' until the approximate
% relative error is below a specified tolerance 'es' or a maximum number of iterations 'maxit' is reached.
%
% [x,ea,iter] = GaussSeidel(A,b,es,maxit): Gauss Seidel without relaxation
%
% input:
%   A = coefficient matrix (must be square)
%   b = right hand side vector
%   es = stop criterion (desired approximate percent relative error, default = 0.00001%)
%   maxit = maximum allowable iterations (default = 50)
%
% output:
%   x = solution vector
%   ea = final maximum approximate percent relative error (%) across all elements of x
%   iter = number of iterations performed
%
% Cases where the method works well:
% - When the system is **diagonally dominant** (the absolute value of the diagonal element in each row is greater than the sum of the absolute values of the off-diagonal elements in that row) Diagonal dominance guarantees convergence.
% - When the system is large and sparse (contains many zero elements).
% - The Gauss-Seidel method generally converges faster than the Jacobi method because it uses the most recently calculated values of the unknowns in the current iteration.
%
% Potential problems/Limitations:
% - **Convergence is not guaranteed.** The method may **diverge** if the system is not diagonally dominant.
% - An initial guess is needed, and a poor guess can lead to slower convergence or divergence.
% - The method only works for a square coefficient matrix A.

% COMMENT ON REARRANGING FOR DIAGONAL DOMINANCE
% -------------------------------------------------------------------------
% The Gauss-Seidel method is only guaranteed to converge if the coefficient
% matrix A is diagonally dominant. If the original system [A]{x}={b} is NOT
% diagonally dominant, you may be able to rearrange the rows (equations)
% to achieve this condition, which is a required step before calling this function.
%
% **Diagonal Dominance Condition:**
% For every row i, the absolute value of the diagonal element |a_ii| must be
% greater than the sum of the absolute values of the other elements in that row:
%
% |a_ii| > sum(|a_ij|) for j != i
%
% **Example Rearrangement:**
% Consider a non-diagonally dominant system:
%
% Original System:
% Equation 1: [ 1x + 7y - 3z = 5] -> |1| < |7| + |-3|  (FAIL)
% Equation 2: [ 8x - 2y + 1z = 1] -> |8| > |-2| + |1|  (PASS)
% Equation 3: [-1x + 1y + 5z = 9] -> |5| > |-1| + |1|  (PASS)
%
% The problem is Row 1. We must find an equation where the 'x' coefficient
% is the largest, and swap it into Row 1. Equation 2 has the largest 'x'
% coefficient (8).
%
% Rearranged System (Swap Row 1 and Row 2):
% Equation 1: [ 8x - 2y + 1z = 1] -> |8| > |-2| + |1|  (PASS)
% Equation 2: [ 1x + 7y - 3z = 5] -> |7| > |1| + |-3|  (PASS, now for 'y')
% Equation 3: [-1x + 1y + 5z = 9] -> |5| > |-1| + |1|  (PASS)
%
% This rearranged matrix is now diagonally dominant and should be passed
% into the GaussSeidel function for reliable convergence.
% -------------------------------------------------------------------------

if nargin<2,error('at least 2 input arguments required'),end
if nargin<4|isempty(maxit),maxit=50;end
if nargin<3|isempty(es),es=0.00001;end
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
if x(i) ~= 0
ea(i) = abs((x(i) - xold(i))/x(i)) * 100;
end
end
iter = iter+1;
if max(ea)<=es | iter >= maxit, break, end
end
ea=max(ea);
