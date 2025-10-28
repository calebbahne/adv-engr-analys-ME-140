function [x,f,ea,iter]=newtmult(func,x0,es,maxit,varargin)
% newtmult: Newton-Raphson root zeroes nonlinear systems
% [x,f,ea,iter]=newtmult(func,x0,es,maxit,p1,p2,...):
% [x,f,ea,iter]=newtmult(@jfreact2,x0); [have u&v w/ ftns]
% uses the Newton-Raphson method to find the roots of
% a system of nonlinear equations
% input:
%   func = name of function that returns f and J
%   x0 = initial guess
%   es = desired percent relative error (default = 0.0001%)
%   maxit = maximum allowable iterations (default = 50)
%   p1,p2,... = additional parameters used by function
% output:
%   x = vector of roots
%   f = vector of functions evaluated at roots (residual, ~ 0)
%   ea = approximate percent relative error (%)
%   iter = number of iterations

% newtmult: Newton-Raphson method for solving systems of nonlinear equations.
%
% This method is significantly faster than single-variable methods (like
% bisection or false position) but requires the calculation of the Jacobian
% matrix (the matrix of partial derivatives).
%
% [x,f,ea,iter]=newtmult(func,x0,es,maxit,p1,p2,...):
%
% input:
%   func = name of function handle that returns both the Jacobian (J) and the function vector (f).
%          (e.g., @jfreact2)
%   x0 = initial guess vector for the roots (e.g., [x1_guess; x2_guess; ...])
%   es = desired approximate percent relative error (default = 0.0001%)
%   maxit = maximum allowable iterations (default = 50)
%   p1,p2,... = optional additional parameters passed to the function 'func'.
%
% output:
%   x = vector of calculated roots (the solution) [columns are pairs]
%   f = vector of functions evaluated at the final roots (should be close to zero)
%   ea = final maximum approximate percent relative error (%), based on the change in x
%   iter = number of iterations performed
%
% -------------------------------------------------------------------------
% CONVERGENCE AND USAGE NOTES
% -------------------------------------------------------------------------
% 1. Convergence Rate: The method exhibits **quadratic convergence** (very fast),
%    meaning the number of significant digits approximately doubles with each
%    iteration, *provided* the initial guess is close enough to the true root.
%
% 2. Potential Problems (Failures):
%    - **Poor Initial Guess:** If x0 is far from the true root, the method may
%      diverge or converge to an unintended root.
%    - **Singular Jacobian:** If the Jacobian matrix J becomes singular (non-invertible)
%      at any point during the iteration, the calculation (dx = J\f) fails, and
%      the program will halt with a division-by-zero or matrix singularity error.
%      This often happens far from the root or at a point where the function's
%      surface is flat (a local extremum).
%    - **Computational Cost:** Calculating and solving the Jacobian (J\f) at every
%      step is computationally intensive compared to simpler methods.
% -------------------------------------------------------------------------

% **Example Case: **
% System: (x-4)^2 + (y-4)^2 = 5  and  x^2 + y^2 = 16
%
% The required function files (u.m and v.m) would look like this:
%
% %--- File: u.m (for f1 = 0) ---
% function result = u(x1, x2)
%   result = (x1-4)^2 + (x2-4)^2 - 5;
% end
%
% %--- File: v.m (for f2 = 0) ---
% function result = v(x1, x2)
%   result = x1^2 + x2^2 - 16;
% end

if nargin<2,error('at least 2 input arguments required'),end
if nargin<3|isempty(es),es=0.0001;end
if nargin<4|isempty(maxit),maxit=50;end
iter = 0;
x=x0; % Initialize the current root estimate vector

while (1)
    % 1. Evaluate the Jacobian matrix J and the function vector f at the current x.
    % The user-supplied function 'func' must return both J and f.
    [J,f]=func(x,varargin{:});

    % 2. Solve the linear system J * dx = -f for the correction vector dx.
    % MATLAB's backslash operator (J\f) efficiently solves J * dx = f for dx.
    % This is the most computationally expensive step, equivalent to:
    % dx = inv(J) * f  (though backslash avoids explicit inversion).
    dx=J\f;

    % 3. Update the root estimate: x_new = x_old - dx
    x=x-dx;

    % 4. Update iteration counter
    iter = iter + 1;

    % 5. Calculate Approximate Percent Relative Error (ea)
    % This finds the maximum relative change across all elements in the root vector x.
    ea=100*max(abs(dx./x));

    % 6. Check for convergence or max iterations
    if iter>=maxit|ea<=es, break, end
end
