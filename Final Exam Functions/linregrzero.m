function [a1] = linregrzero(x,y)
% linregrzero: linear regression curve fitting with zero intercept
%   [a, r2] = linregr(x,y): Least squares fit of straight
%   line to data with zero intercept
% input:
%   x = independent variable
%   y = dependent variable
% output:
%   a1 = slope
n = length(x);
if length(y)~=n, error('x and y must be same length'); end
x = x(:); y = y(:); % convert to column vectors
sx2 = sum(x.*x); sxy = sum(x.*y);
a1 = sxy/sx2;
% create plot of data and best fit line
xp = linspace(0,max(x),2);
yp = a1*xp;
plot(x,y,'o',xp,yp)
grid on