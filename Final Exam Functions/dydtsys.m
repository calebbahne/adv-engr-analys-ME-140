function dy = dydtsys(t,y)
dy = y*t.^2 - 1.1*y;

%{
% Example implementation
function dy = dydtsys(t,y)
dy = [y(2);
    9.81-0.25/68.1*y(2)^2];

tspan = [0 15];
y0 = [1000; 0];   % 1000 m height, zero initial velocity
h = 0.05;

[t,y] = rk4(@dydtsys, tspan, y0, h);

height = y(:,1);
velocity = y(:,2);
%}
