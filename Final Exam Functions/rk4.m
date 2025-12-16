function [t,y] = rk4(dydt,tspan,y0,h)
% [t,y] = rk4(dydt,tspan,y0,h):
% [t,y] = rk4(@dydtsys, tspan, y0, h);
%   uses the fourth-order Runge-Kutta method to integrate an ODE
% input:
%   dydt = name of the M-file that evaluates the ODE 
%   tspan = [ti, tf] where ti and tf = initial and
%           final values of independent variable
%   y0 = initial value of dependent variable
%   h = step size
% output:
%   t = vector of independent variable
%   y = vector of solution for dependent variable

ti = tspan(1);
tf = tspan(2);
t = (ti:h:tf)';
n = length(t);
% if necessary, add an additional value of t
% so that range goes from t = ti to tf
if t(n)<tf  
  t(n+1) = tf;
  n = n+1;
end
y = y0*ones(n,1); %preallocate y to improve efficiency
for i = 1:n-1  
  hh = t(i+1) - t(i);  
  k1 = feval(dydt,t(i),y(i));   
  ymid = y(i) + k1*hh/2;
  k2 = feval(dydt,t(i)+hh/2,ymid);   
  ymid = y(i) + k2*hh/2;
  k3 = feval(dydt,t(i)+hh/2,ymid);   
  yend = y(i) + k3*hh;
  k4 = feval(dydt,t(i)+hh,yend);   
  phi = (k1+2*(k2+k3)+k4)/6;
  y(i+1) = y(i) + phi*hh;
end
plot(t,y)

% -------------------------------------------------------------------------
% HOW TO USE THIS FUNCTION:
%
% 1) Create a separate function (M-file) that defines dy/dt.
%    This function must accept two inputs: t (time) and y (state),
%    and return dy/dt.
%
% 2) Call rk4 with:
%       - the function handle or name of your dydt function
%       - a time span [ti tf]
%       - an initial condition y0
%       - a step size h
%
% EXAMPLE:
%   Solve the ODE:
%       dy/dt = -2y
%       y(0) = 1
%       over the interval 0 ≤ t ≤ 5
%
% Step 1: Create a file called dydtsys.m:
% -------------------------------------------------
% function dy = dydtsys(t,y)
%     dy = -2*y;
% end
% -------------------------------------------------
%
% Step 2: Call rk4 from the command window or a script:
% -------------------------------------------------
% tspan = [0 5];
% y0 = 1;
% h = 0.1;
% [t,y] = rk4(@dydt_s, tspan, y0, h);
% -------------------------------------------------
%
% The function will return vectors t and y and automatically
% plot y versus t.




% -------------------------------------------------------------------------
% EXAMPLE (SYSTEM OF ODEs):
%
%   Solve the system:
%       dy1/dt = y2
%       dy2/dt = -y1
%
%   with initial conditions:
%       y1(0) = 1
%       y2(0) = 0
%
%   This represents simple harmonic motion.
%
% IMPORTANT:
%   For systems, y must be a COLUMN VECTOR and dydt must return
%   a COLUMN VECTOR of the same size.
%
% Step 1: Create a file called dydt_system.m:
% -------------------------------------------------
% function dydt = dydt_system(t,y)
%     dydt = zeros(2,1);
%     dydt(1) = y(2);
%     dydt(2) = -y(1);
% end
% -------------------------------------------------
%
% Step 2: Call rk4:
% -------------------------------------------------
% tspan = [0 10];
% y0 = [1; 0];     % initial conditions as a column vector
% h = 0.01;
% [t,y] = rk4(@dydt_system, tspan, y0, h);
%
% % To plot individual states:
% figure
% plot(t,y(:,1),'b-',t,y(:,2),'r--')
% legend('y_1','y_2')
% xlabel('t')
% ylabel('State values')
% grid on
% -------------------------------------------------
%
% NOTE:
%   If using systems, remove or comment out the line:
%       plot(t,y)
%   inside rk4.m to avoid dimension errors.
% -------------------------------------------------------------------------


% RK4 advances the solution in fixed time steps h from t_i to t_f
% At each step, four slope estimates (k1–k4) are computed within the interval
% k1 uses the slope at the start, k2 and k3 at the midpoint, k4 at the end
% These slopes are combined in a weighted average to estimate the net change
% The state y is updated and stored at the next time value
