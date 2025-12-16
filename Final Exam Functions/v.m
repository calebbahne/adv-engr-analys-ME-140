function fv = v(x1,x2)
% v: fu = v(x1,x2)
%   This is required for a poorly written newtmult code. 
%   Have to pass in like this instead of as parameters.
fv = x2+5*x1*x2-x1^2;
end