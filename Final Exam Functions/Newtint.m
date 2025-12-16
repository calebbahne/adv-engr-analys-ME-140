function yint = Newtint(x,y,xx)
% Newtint: Newton interpolating polynomial
% yint = Newtint(x,y,xx): Uses an (n order Newton
% interpolating polynomial based on n data points (x, y)
% to determine a value of the dependent variable (yint)
% at a given value of the independent variable, xx.
% input:
%   x = independent variable
%   y = dependent variable
%   xx = value of independent variable at which
%           interpolation is calculated
% output:
%   yint = interpolated value of dependent variable
%       compute the finite divided differences in the form of a
%       difference table
% x & y must be same length

% METHOD OVERVIEW:
%   This method constructs the interpolating polynomial incrementally
%   using finite divided differences. Unlike the Lagrange form, which
%   must be fully recomputed if a new point is added, the Newton form
%   can easily incorporate additional data by extending the divided
%   difference table — making it more efficient for iterative problems.
%
% -------------------------------------------------------------------------
% WHEN IT WORKS WELL:
%   - When you need to **add new data points** without recomputing the
%     entire polynomial.
%   - When the **data set is large**, since divided differences are
%     computationally efficient.
%   - When you need **better numerical stability** compared to the
%     Lagrange method.
%   - Suitable for **smooth, continuous functions** with well-spaced data.
%
% WHEN IT DOESN'T WORK WELL:
%   - Like all high-order interpolating polynomials, it still suffers
%     from **Runge’s phenomenon** if too many points are used.
%   - **Unevenly spaced** or rapidly varying data can lead to large
%     interpolation errors.
%   - **Noisy data** may cause unrealistic oscillations.
%
% -------------------------------------------------------------------------
% COMPARISON TO OTHER METHODS:
%   * Lagrange Interpolation:
%       - Conceptually simpler, but must recompute from scratch if
%         new points are added.
%       - Less numerically efficient for large data sets.
%
%   * Polynomial Interpolation (polyfit/polyval):
%       - Fits a polynomial in a least-squares sense, not exact at
%         all points (better for noisy data).
%       - Generally more stable for large data sets or approximate fits.
%
% -------------------------------------------------------------------------
% EXAMPLE:
%   x = [1 2 4];
%   y = [1 4 16];
%   xx = 3;
%   yint = Newtint(x,y,xx)
% -------------------------------------------------------------------------

n = length(x);
if length(y)~=n, error('x and y must be same length'); end
b = zeros(n,n);
% assign dependent variables to the first column of b.
b(:,1) = y(:); % the (:) ensures that y is a column vector.
for j = 2:n
for i = 1:n-j+1
b(i,j) = (b(i+1,j-1)-b(i,j-1))/(x(i+j-1)-x(i));
end
end
% use the finite divided differences to interpolate
xt = 1;
yint = b(1,1);
for j = 1:n-1
xt = xt*(xx-x(j));
yint = yint+b(1,j+1)*xt;
end
