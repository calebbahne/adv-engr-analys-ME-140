function [J, f] = jfreact2_v2(x, u, v, varargin)
%JFREACT2   Compute the Jacobian and function vector for multivariable Newton-Raphson.
%
%   [J, f] = jfreact2(x, u, v)
%
%   Numerically estimates the Jacobian matrix J and evaluates the function vector f
%   for a system of two nonlinear equations f1(x1,x2) = 0 and f2(x1,x2) = 0.
%
%   INPUTS:
%       x  - 2×1 vector [x1; x2], current variable estimates.
%       u  - Function handle for the first equation, f1(x1,x2).
%       v  - Function handle for the second equation, f2(x1,x2).
%       varargin - (optional) extra parameters to pass to u and v.
%
%   OUTPUTS:
%       J  - 2×2 Jacobian matrix of partial derivatives.
%       f  - 2×1 vector of function values [f1; f2].
%
%   The Jacobian is computed using finite differences:
%       J(i,j) ≈ [f_i(x + Δx_j) - f_i(x)] / Δx_j
%
%   Example:
%       u = @(x,y) (x-4).^2 + (y-4).^2 - 5;
%       v = @(x,y) (x-6).^2 + (y-2).^2 - 5;
%       x0 = [3; 3];
%       [J, f] = jfreact2(x0, u, v);
%
% -------------------------------------------------------------------------

    del = 1e-6;

    df1dx1 = (u([x(1)+del*x(1); x(2)]) - u(x)) / (del*x(1));
    df1dx2 = (u([x(1); x(2)+del*x(2)]) - u(x)) / (del*x(2));

    df2dx1 = (v([x(1)+del*x(1); x(2)]) - v(x)) / (del*x(1));
    df2dx2 = (v([x(1); x(2)+del*x(2)]) - v(x)) / (del*x(2));

    J = [df1dx1 df1dx2;
         df2dx1 df2dx2];

    f = [u(x);
         v(x)];
end
