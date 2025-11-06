function yint = Lagrange(x,y,xx)
% Lagrange: Lagrange interpolating polynomial
% yint = Lagrange(x,y,xx): Uses an (n order
%   Lagrange interpolating polynomial based on n data points
%   to determine a value of the dependent variable (yint) at
%   a given value of the independent variable, xx.
% input:
%   x = independent variable
%   y = dependent variable
%   xx = value of independent variable at which the
%       interpolation is calculated
% output:
%   yint = interpolated value of dependent variable

% WHEN IT WORKS WELL:
%   - When the number of data points (n) is small (typically n ≤ 10)
%   - When data points are evenly spaced or mildly uneven
%   - When the function being modeled is smooth and slowly varying
%   - Ideal for quick, one-off interpolation when you don’t need to
%     repeatedly add new points
%
% WHEN IT DOESN'T WORK WELL:
%   - For large n (many data points), Lagrange interpolation suffers from
%     numerical instability and Runge’s phenomenon (oscillations between
%     points).
%   - Every time a new data point is added, the entire polynomial must
%     be recomputed — computationally expensive for iterative data sets.
%   - Poor performance if data are unevenly spaced or highly nonlinear.
%
% COMPARISON TO OTHER METHODS:
%   * Newton Interpolation:
%       - Better suited for adding data incrementally (can reuse divided
%         differences).
%       - Computationally more efficient when adding points.
%       - Numerically more stable for many data points.
%
%   * Polynomial Interpolation (polyfit/polyval):
%       - Uses least-squares fitting rather than exact interpolation.
%       - Better for noisy data where an approximate trend is desired.
%       - Avoids oscillation problems for large n by fitting a lower-degree
%         polynomial instead of an exact match.
%
% -------------------------------------------------------------------------
% EXAMPLE:
%   x = [1 2 4];
%   y = [1 4 16];
%   xx = 3;
%   yint = Lagrange(x,y,xx)
% -------------------------------------------------------------------------
n = length(x);
if length(y)~=n, error('x and y must be same length'); end
s = 0;
for i = 1:n
product = y(i);
for j = 1:n
if i ~= j
product = product*(xx-x(j))/(x(i)-x(j));
end
end
s = s+product;
end
yint = s;
