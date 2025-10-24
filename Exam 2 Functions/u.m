function fu = u(x1,x2)
% u: fu = u(x1,x2)
%   This is required for a poorly written newtmult code. 
%   Have to pass in like this instead of as parameters.
fu = -x1^2+x1+0.75-x2; % set = 0
end