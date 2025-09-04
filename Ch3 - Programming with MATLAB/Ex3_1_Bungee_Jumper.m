% Ex3_1_Bungee_Jumper
%   Can type that into the command window and it'll run this script, fyi
g = 9.81; m = 68.1; t = 12; cd = 0.25;
v = sqrt(g*m/cd)*tanh(sqrt(g*cd/m)*t)