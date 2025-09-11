function vol = volConeCyl(R,d)
%VOLCONECYL: CH3 P03.001 - 1 - Computes volume in a cylinder/cone thing
% inputs:
%   R = height of cone & radius of tank (vector)
%       2R = length of cylinder
%   d = fluid height                    (vector)
% outputs:
%   vol = volume of the piecewise thingy
%   Note: V=πr^2h/3 volume of a cone
% Logic
%   If d <= R, volume of cone up to d
%       45 deg angle -> r = d
%   If d > R and d <= 3R, full vol of cone + (d-R) for cyl
%   If d > 3R, output an error
%   Tricky cuz gotta handle element wise

% DOESN'T WORK: not set up to handle vectors

if d < 0
    error('d must be positive');
elseif d <= R
    r = d; % radius of cone
    h = d; % height of cone
    H = 0; % height in cylinder
elseif (d > R) & (d <= 3*R)
    r = R; 
    h = R;
    H = (d-R);
else
    r = -999; 
    h = R;
    H = (d-R);
    fprintf('Overtop.');
end

vol = pi.*r.^2.*h/3 + pi.*R.^2.*H; % cone + cyl
end