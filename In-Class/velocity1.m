function vend = velocity1(dt, ti, tf, vi)
% velocity1: Euler solution for bungee velocity
%   vend - velocity1(dt, ti, tf, vi)
%   input
%       dt = time increment
%       ti = start time
%       tf = end time
%       vi = initial velo

t = ti;
v = vi;
n = (tf-ti) / dt;
for i = 1:n             % should replace w/ a while loop to handle ~= divis
    dvdt = deriv(v);
    v = v+dvdt*dt;
    t = t+dt;
end
vend = v;
end

function dv = deriv(v)
dv=9.81-(0.25/68.1)*v*abs(v);
end
