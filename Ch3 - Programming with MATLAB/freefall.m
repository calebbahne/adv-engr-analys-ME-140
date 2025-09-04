function v = freefall(t, m, cd)
%FREEFALL: bungee velocity with second order drag
% v=freevall(t,m,cd) computes the free-fall velocity 
%                    of an object with second-order drag
% input:
%   t = time (s)
%   m = mass (kg)
%   cd = second-order drag coefficient (kg/m)
% output:
%   v = downward velocity (m/s)

g = 9.81; % accel of grav
v = sqrt(g*m./cd).*tanh(sqrt(g*cd./m).*t);
end