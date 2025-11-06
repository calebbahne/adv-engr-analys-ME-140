function [Ra, T_s, T_inf] = getRa(fluid, ra_type)
% getRa Calculate Rayleigh number using Grashof number
%   [Ra, T_s, T_inf] = getRa(fluid, ra_type)
%
%   fluid   : 'air' or 'water'
%   ra_type : 'Ra_D'/'D' (use diameter) or 'Ra_L'/'L' (use length)
%
%   Calculates: Ra = Gr * Pr
%   where Gr = g * beta * (T_s - T_inf) * L^3 / (nu^2)
%   and beta = 1 / T_K  (ideal gas approximation)
%

fluid = lower(string(fluid));
ra_type = lower(string(ra_type));

% Determine characteristic name for prompts
if contains(ra_type, 'd')
    charName = 'diameter (D)';
else
    charName = 'length (L)';
end

% Prompt user to select input method
disp('Select Ra input method:');
fprintf('   1. Compute from (%s)\n', charName);
disp('   2. Direct input (Ra)');
RaChoice = input('Ra input method: ');
clc;

if RaChoice == 2
    Ra = input('Enter Rayleigh number (Ra): ');
    return
end
clc;

% Prompt for surface and ambient temperatures
T_s = input('Enter surface temperature, T_s (C): ');
T_inf = input('Enter ambient temperature, T_inf (C): ');

T_s = T_s +273.15;
T_inf = T_inf + 273.15;

% Prompt for characteristic dimension
L = input(sprintf('Enter characteristic %s (m): ', charName));

clc;

T_f = mean([T_s T_inf]);

% Compute property values based on fluid type
if fluid == "air"
    nu = getAirProp(T_f, 'nu');       % m^2/s
    alpha = getAirProp(T_f, 'alpha'); % m^2/s
    Pr = getAirProp(T_f, 'pr');
elseif fluid == "water"
    nu = getWaterProp(T_f, 'nu');
    alpha = getWaterProp(T_f, 'alpha');
    Pr = getWaterProp(T_f, 'pr');
else
    error('Unknown fluid: use ''air'' or ''water''.');
end

% Constants
g = 9.81;                  % gravitational acceleration (m/s^2)
beta = 1 / T_f;            % thermal expansion coefficient (1/K)
deltaT = T_s - T_inf;      % temperature difference

% Compute Grashof number
Gr = g * beta * deltaT * L^3 / (nu^2);

% Compute Rayleigh number
Ra = abs(Gr * Pr);

end