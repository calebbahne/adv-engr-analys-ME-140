%% In-Class 3 - 9.9 Rocket Velocity
velo_4 = rocketVelo(4)
velo_375 = rocketVelo(37.5)

function velo = rocketVelo(t)
%ROCKETVELO: velocity of a rocket at time t (piecewise ftns)
% Inputs:
%   t = time (sec)
% Outputs:
%   velo = velocity (m/s)

if t < 0
    velo = 0; % time must be non-negative
elseif t > 26
    velo = 2136*exp(-0.1*(t-26));
elseif t >= 16
    velo = 36*t+12*(t-16)^2;
elseif t >= 8
    velo = 624 - 3*t;
else % 0 to 8
    velo = 10*t^2-5*t;
end
end