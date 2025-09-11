function vol = volConeCyl_v2(R,d)
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

% Set up for matrix > for loop
% For = easier to code, but slower

if size(d) ~= size(R)
    error('Matrices must be same size.');
end

vol = nan(size(d));

% Case 1: just in cone
idx_1 = (d < R) & (d > 0); % [1 0]s, 1 & for vector not scalar
vol(idx_1) = pi.*d(idx_1).^2.*d(idx_1)/3; % r = d, h = d

% Case 2: in the cylinder
idx_2 = (d > R) & (d <= 3*R);
vol(idx_2) = pi.*R(idx_2).^2.*R(idx_2)/3 + pi.*R(idx_2).^2.*(d(idx_2)-R(idx_2));

% Case 3: overflowing
idx_3 = (d > 3*R);
vol(idx_3) = -1; % flag
fprintf('Warning: some entries overtopped the tank.\n');

end