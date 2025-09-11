function [r,th] = cart2polar(x,y)
%CART2POLAR: converts cartesian to polar form
% inputs: exactly what you think
% outputs: r and theta

if size(x) ~= size(y)
    error('x and y must be the same size');
end

r = sqrt(x.^2+y.^2);

% Initialize
th = nan(size(x));

% Q1 and Q4, nonzero y
idx1 = x > 0 & y ~= 0;
th(idx1) = atan(y(idx1)./x(idx1));

% Q2
idx2 = x < 0 & y > 0;
th(idx2) = atan(y(idx2)./x(idx2)) + pi;

% Q3
idx3 = x < 0 & y < 0;
th(idx3) = atan(y(idx3)./x(idx3)) - pi;

% Vertical up
idxV_up = x == 0 & y > 0;
th(idxV_up) = pi/2;

% Vertical down
idxV_down = x == 0 & y < 0;
th(idxV_down) = -pi/2;

% Zero
idxZ = x == 0 & y == 0;
th(idxZ) = 0;

% Left
idxL = x < 0 & y == 0;
th(idxL) = pi;

% Right
idxR = x > 0 & y == 0;
th(idxR) = 0;

th = rad2deg(th);

end