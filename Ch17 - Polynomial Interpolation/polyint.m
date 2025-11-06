function yint = polyint(x,y,xx)
% polyint: Polynomial interpolation
%   yint = polyint(x,y,xx): Uses an (n - 1)-order
%   polynomial based on n data points (x, y)
%   to determine a value of the dependent variable (yint)
%   at a given value of the independent variable, xx.
% input:
%   x = independent variable
%   y = dependent variable
%   xx = value of independent variable at which
%       interpolation is calculated
% output:
%   yint = interpolated value of dependent variable

% -------------------------------------------------------------------------
% METHOD OVERVIEW:
%   polyfit(x,y,n-1) computes the coefficients of a polynomial of
%   degree (n - 1) that passes exactly through the n given data points.
%   polyval(p,xx) then evaluates that polynomial at one or more
%   specified points.
%
%   This method is a general polynomial interpolation approach that
%   provides a concise and efficient implementation in MATLAB.
%
% -------------------------------------------------------------------------
% WHEN IT WORKS WELL:
%   - When the number of data points is **moderate (n ≤ 10–12)**.
%   - When the underlying function is **smooth and continuous**.
%   - For **quick polynomial fitting** or when a closed-form polynomial
%     expression is desired.
%   - Excellent for **symbolic or analytical purposes** (e.g., to find
%     coefficients or derivatives of the fitted polynomial).
%
% WHEN IT DOESN'T WORK WELL:
%   - For large n (many data points), the polynomial can become unstable
%     due to rounding errors and **Runge’s phenomenon** (oscillations
%     between points).
%   - Performs poorly if data are **unevenly spaced** or **highly nonlinear**.
%   - Not ideal for **noisy data**, since the polynomial tries to pass
%     exactly through all points — overfitting can occur.
%
% -------------------------------------------------------------------------
% COMPARISON TO OTHER METHODS:
%   * Lagrange Interpolation:
%       - Similar mathematical foundation, but computed explicitly
%         using Lagrange basis polynomials.
%       - Must recompute from scratch if new points are added.
%
%   * Newton Interpolation:
%       - More efficient if data are added incrementally.
%       - Numerically more stable for large data sets due to the use of
%         divided differences.
%
%   * Polynomial Regression (polyfit with degree < n-1):
%       - When data are noisy, fitting a lower-degree polynomial gives a
%         smoother and more general trend rather than exact interpolation.
%
% -------------------------------------------------------------------------
% EXAMPLE:
%   x = [1 2 4];
%   y = [1 4 16];
%   xx = 3;
%   yint = polyint(x,y,xx)
% -------------------------------------------------------------------------

n = length(x);
if length(y)~=n, error('x and y must be same length'); end
p = polyfit(x,y,n-1);
yint = polyval(p,xx);
end
