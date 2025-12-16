function [a, r2] = linregr(x,y)
% linregr: linear regression curve fitting
%   [a, r2] = linregr(x,y): Least squares fit of straight
%   line to data by solving the normal equations
% input:
%   x = independent variable
%   y = dependent variable
% output:
%   a = vector of slope, a(1), and intercept, a(2)
%   r2 = coefficient of determination
%
% linearization help:
%   y = alph * exp(beta*x)
%       ln(y) = ln(alph) + beta*x 
%           Note: use log() in MATLAB
%       [a, r2] = linregr(x, log(y));
%       alpha = exp(a(2));
%       beta  = a(1);
%
%   y = alph * x^beta
%       log10(y) = log10(alph) + beta*log10(x)
%       [a, r2] = linregr(log10(x), log10(y));
%       alpha = 10^(a(2));
%       beta  = a(1);
%
%   y = alph * x/(beta + x)
%       1/y = 1/alph + (beta/alph)*(1/x)
%         (linear in 1/x)
%       [a, r2] = linregr(1./x, 1./y);
%       alpha = 1 / a(2);
%       beta  = a(1) / a(2);

n = length(x);
if length(y)~=n, error('x and y must be same length'); end
x = x(:); y = y(:); % convert to column vectors
sx = sum(x); sy = sum(y);
sx2 = sum(x.*x); sxy = sum(x.*y); sy2 = sum(y.*y);
a(1) = (n*sxy-sx*sy)/(n*sx2-sx^2);
a(2) = sy/n-a(1)*sx/n;
r2 = ((n*sxy-sx*sy)/sqrt(n*sx2-sx^2)/sqrt(n*sy2-sy^2))^2;
% create plot of data and best fit line
xp = linspace(min(x),max(x),2);
yp = a(1)*xp+a(2);
plot(x,y,'o',xp,yp)
grid on
