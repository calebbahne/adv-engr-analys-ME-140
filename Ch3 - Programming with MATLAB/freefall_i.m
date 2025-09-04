function freefall_i
%FREEFALL_I: interactive bungee velocity with second order drag
% v=freevall(t,m,cd) computes the free-fall velocity 
%                    of an object with second-order drag
%                    using user input.
% No input/output variables, purely to the screen

g = 9.81; % accel of grav

m = input('Mass (kg): ');
cd = input('Drag coefficient (kg/m): ');
t = input('Time (s): ');
disp(' '); % works like a newline
disp('Velocity (m/s): ');
disp(sqrt(g*m/cd)*tanh(sqrt(g*cd/m)*t));
end