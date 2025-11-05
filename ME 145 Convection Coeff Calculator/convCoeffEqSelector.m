%% Convection Coefficient Calculator
% Purpose: Calculate the avg convection coefficient h for a given condition
%   Determines the appropriate case, then does calcs from there.

clc; clear;

disp('Convection Coefficient Calculator');
disp('   - This code calculates the avg convection coefficient (h) for common heat & mass situations.');
disp('   - It selects the appropriate case, calculates the avg Nusselt (Nu), then solves for h.');
disp(' ');
input('Press ENTER to continue...','s');
clc;

fluid = getFluid();
clc;

%% Overarching Convection Case
disp('Select the convection case:');
disp('  1. Forced External Flow');
disp('  2. Forced Internal Flow'); 
disp('  3. Free Convection');
convCase = input('Convection case (1-3): ');

switch convCase 

    case 1  % External Flow
        clc;
        disp('Select an external flow:');
        disp('  1. Flat plate in parallel flow');
        disp('  2. Cylinder in cross flow');
        disp('  3. Sphere');
        disp('  4. Bank of tubes');
        extConvCase = input('External flow case (1–4): ');
        
        switch extConvCase
            case 1 
                clc;
                disp('External conv, flat plate');
                % From here, get Re and Pr
                % Then go ftn to find Nu (we'll find h @ end)
                % using ftns for this will decrease the amt of stuff we've
                %   gotta put here, making the code easier to interpret

                T_K = getTemp(); % Tfilm
                Re_L = getRe(T_K, fluid, 'Re_L');
                Pr = getPr(fluid, T_K); % make sure it's T_f

                [Nu_L, convType] = ExtConvFlatPlate(Re_L, Pr);
            case 2
                clc;
                disp('External conv, cylinder');

                T_K = getTemp(); % Tfilm
                Re_D = getRe(T_K, fluid, 'Re_D');
                Pr = getPr(fluid, T_K); % make sure it's T_f


                Nu_D = ExtConvCyl(Re_D, Pr);
            case 3
                clc;
                disp('External conv, sphere');

                T_K = getTemp(); % T_inf
                Re_D = getRe(T_K, fluid, 'Re_D');
                disp('***REMOVE*** Make sure its T_inf for Pr');
                Pr = getPr(fluid, T_K); % make sure it's T_inf
                mu = input('Dynamic viscosity (mu): ');
                mu_s = input('Dynamic viscosity @ surface (mu_s):');

                Nu_D = ExtConvSphere(Re_D, Pr, mu, mu_s);
            case 4
                clc;
                disp('External conv, Bank of tubes');

                T_K = getTemp(); % Tfilm
                Re_Dmax = getRe(T_K, fluid, 'Re_D');

                disp('***REMOVE*** Make sure its T_bar for Pr');
                Pr = getPr(fluid, T_K); % make sure it's T_bar
                disp('***REMOVE*** Make sure its T_s for Pr_s')
                Pr_s = getPr(fluid, T_K); % make sure it's T_bar **** update to take in temp too

                disp('Select a tube configuration:');
                disp('  1. Aligned');
                disp('  2. Staggered');
                tubeType = input('Tube config: ');

                N_L = input('Longitudinal pitch (N_L): ');
                mu_s = input('Transverse pitch (N_T):');

                Nu_D = ExtConvTubeBank(Re_Dmax, Pr, Pr_s, tubeType, N_L, N_T);
            otherwise
                disp('Invalid selection.');
        end

    case 2  % Internal Flow
        clc;
        disp('Select an internal flow:');
        disp('  1. Circular tube');
        disp('  2. Noncircular tube');
        disp('  3. Concentric tube annulus');
        intConvCase = input('Internal flow case (1–3): ');

        switch intConvCase
            case 1 % Circular tube
                clc;
                disp('Select a circular tube flow:');
                disp('  1. Laminar flow');
                disp('  2. Turbulent flow');
                circTubeFlow = input('Circular tube flow (1–2): ');
                if circTubeFlow == 1
                    disp('Int conv, laminar flow in circular tube');
                elseif circTubeFlow == 2
                    disp('Int conv, turbulent flow in circular tube');
                else
                    disp('Invalid selection.');
                end

            case 2 % Noncircular tube
                clc;
                disp('Select a noncircular tube:');
                disp('  1. Square/rectangle');
                disp('  2. Triangle');
                nonCircGeom = input('Select geometry (1–2): ');
                if nonCircGeom == 1
                    clc;
                    disp('Int conv, noncirc, square/rectangular duct');
                    % Use b/a to determine subcase.
                    % could be square, rect, inf rect, etc
                elseif nonCircGeom == 2
                    clc;
                    disp('Int conv, noncirc, triangular duct');
                else
                    clc;
                    disp('Invalid selection.');
                end

            case 3 % Concentric Tube Annulus
                clc;
                disp('Int conv, noncirc, concentric tube annulus');
            otherwise
                disp('Invalid selection.');
        end

    case 3  % Free Convection
        clc;
        disp('Select a free convection:');
        disp('  1. Immersed free convection');
        disp('  2. Enclosed free convection');
        freeCase = input('Free convection case (1–2): ');

        switch freeCase
            case 1  % External
                clc;
                disp('Select an immersed free convection:');
                disp('  1. Flat plate');
                disp('  2. Horizontal cylinder');
                disp('  3. Sphere');
                immFree = input('Immersed free conv case (1–3): ');

                switch immFree
                    case 1  % Flat plate
                        clc;
                        disp('Select a flat plate orientation:');
                        disp('  1. Vertical');
                        disp('  2. Horizontal');
                        disp('  3. Inclined');
                        flatCaseOrient = input('Flat plate case (1–3): ');

                        switch flatCaseOrient
                            case 1
                                clc;
                                disp('Free conv, immersed, vert flat plate');

                                % Get Ra_L, Pr
                                convType = 'Length'; % change
                                [Ra_L, T_s, T_inf] = getRa(fluid, convType);
                                T_f = mean([T_s T_inf]);
                                Pr = getPr(fluid, T_f);
                                
                                % Check precondit
                                Nu_L = freeConvExtFlatPlateVert(Ra_L, Pr);
                            case 2
                                clc;
                                disp('Select a horizontal plate:');
                                disp('  1. Hot surface up / cold surface down');
                                disp('  2. Cold surface up / hot surface down');
                                hotcold = input('Orientation (1–2): ');
                                if hotcold == 1
                                    clc;
                                    disp('Free conv, immersed, flat plate, hot surface up / cold surface down');
                                    convType = 'Ra_L';

                                    [Ra_L, T_s, T_inf] = getRa(fluid, convType);
                                    T_f = mean([T_s T_inf]);
                                    Pr = getPr(fluid, T_f);
                                    Nu_L = freeConvExtPlateHorizHotUpper(Ra_L,Pr);
                                elseif hotcold == 2
                                    clc;
                                    disp('Free conv, immersed, flat plate, cold surface up / hot surface down');
                                    convType = 'Ra_L';

                                    [Ra_L, T_s, T_inf] = getRa(fluid, convType);
                                    T_f = mean([T_s T_inf]);
                                    Pr = getPr(fluid, T_f);
                                    Nu_L = freeConvExtPlateHorizHotLower(Ra_L,Pr);
                                else
                                    clc;
                                    disp('Invalid selection.');
                                end
                            case 3
                                clc;
                                disp('Select an inclined flat plate:');
                                disp('  1. Cold surface up / hot surface down');
                                disp('  2. Hot surface up / cold surface down');
                                incline = input('Select configuration (1–2): ');
                                if incline == 1
                                    clc;
                                    disp('Free conv, immersed, inclined plate, old surface up / hot surface down');

                                    theta = input('Inclination angle theta (deg):');
                                    convType = 'Ra_L';
                                    %   Ra calc'd with g*cos(theta) *****
                                    [Ra_L, T_s, T_inf] = getRa(fluid, convType);
                                    T_f = mean([T_s T_inf]);
                                    Pr = getPr(fluid, T_K);

                                    Nu_L = freeConvExtPlateInc(Ra_L,Pr,theta);
                                else
                                    clc;
                                    disp('Not currently supported (not in textbook)');
                                    % this one's not in table 9.3
                                end
                            otherwise
                                disp('Invalid selection.');
                        end

                    case 2
                        clc;
                        disp('Free conv, immersed, horizontal cylinder');
                        convType = 'Ra_D';
                        
                        [Ra_D, T_s, T_inf] = getRa(fluid, convType);
                        T_f = mean([T_s T_inf]);
                        Pr = getPr(fluid, T_f);
                        Nu_D = freeConvExtHorizCyl(Ra_D,Pr);
                    case 3
                        clc;
                        disp('Free conv, immersed, sphere');
                        convType = 'Ra_D';
                        
                        [Ra_D, T_s, T_inf] = getRa(fluid, convType);
                        T_f = mean([T_s T_inf]);
                     
                        Pr = getPr(fluid, T_K);
                        
                        Nu_D = freeConvExtSphere(Ra_D,Pr);
                    otherwise
                        clc;
                        disp('Invalid selection.');
                end

            case 2  % Internal free convection
                clc;
                disp('Select an enclosed free convection geometry:');
                disp('  1. Rectangular cavity');
                disp('  2. Concentric cylinders');
                disp('  3. Concentric spheres');
                encFreeGeom = input('Geometry (1–3): ');

                switch encFreeGeom
                    case 1
                        clc;
                        disp('Select a rectangular cavity orientation:');
                        disp('  1. Vertical');
                        disp('  2. Horizontal');
                        disp('  3. Inclined');
                        cavity = input('Orientation (1–3): ');
                        if cavity == 1
                            clc;
                            disp('Free convect, enclosed, rectangular, vertical rectangular cavity');
                        elseif cavity == 2
                            clc;
                            disp('Free convect, enclosed, rectangular, horizontal rectangular cavity');
                        elseif cavity == 3
                            clc;
                            disp('Free convect, enclosed, rectangular, inclined rectangular cavity');
                        else
                            clc;
                            disp('Invalid selection.');
                        end
                    case 2
                        clc;
                        disp('Free convect, enclosed, concentric cylinders');
                    case 3
                        clc;
                        disp('Free convect, enclosed, concentric spheres');
                    otherwise
                        disp('Invalid selection.');
                end
            otherwise
                clc;
                disp('Invalid selection.');
        end

    otherwise
        disp('Invalid main selection.');
end

disp(' ');
disp('Get stuff for Re, Ra, Pr, etc.')
disp('Calcs for Nu go here');
disp('Then we calc h');
% h = k*Nu_L/L;
% h = k*Nu_D/D;
disp(' ');
disp('Output the following:');
disp('  Selected case type');
disp('  Re, Pr, Ra, geometry stuff, Nu');
disp('  h')

disp('----------------------------------');
disp('Convection coefficient calculation complete');
%disp('   We love you Troy, so please give us the A :)');
disp('----------------------------------');


%% Functions ============================================================

% Get user input
function fluid = getFluid()
% gets the fluid
% currently only water or air, but easy to add more if needed

clc;

disp('Select the fluid to use:');
disp('   1. Air (atmospheric)');
disp('   2. Water (saturated liquid)');
fluidSelect = input('Fluid: ');

switch fluidSelect
    case 1
        fluid = 'air';
    case 2
        fluid = 'water';
    otherwise
        fluid = 'air';
        warning('Invalid entry. Air is autoselected');
end
end

function T_K = getTemp()
% getTemp: Get the temperature from the user
% Inputs
%   Add more, this is a placeholder
T_C = input('Enter the Temperature (C): ');
T_K = T_C+273.15;
clc;
end

function Re = getRe(T_K, fluid, re_type)
% getRe Calculate Reynolds number using nu
%   Re = getRe(T_K, fluid, re_type)
%
%   T_K     : temp in Kelvin
%   fluid   : 'air' or 'water'
%   re_type : 'Re_D'/'D' (use diameter) or 'Re_L'/'L' (use length)
%

fluid = lower(string(fluid)); % turns out format as str helps out
re_type = lower(string(re_type)); % no need for strcompi

% Determine characteristic name for prompts
if contains(re_type, 'd')
    charName = 'diameter (D)';
else
    charName = 'length (L)'; % defaults to length automatically
end

% Ask user whether to compute or directly supply Re
disp('Select Re input method:');
fprintf('   1. Velocity (V), %s\n', charName);
disp('   2. Direct input (Re)');
ReChoice = input('Re input method: ');
clc;

if ReChoice == 2
    Re = input('Enter Re: ');
    return
end
clc;

% Otherwise compute from V and characteristic length
V = input('Enter velocity V (m/s): ');
L = input(sprintf('Enter characteristic %s (m): ', charName));

clc;

% Get kinematic viscosity from property functions
if fluid == "air"
    nu = getAirProp(T_K, 'nu');   % expects m^2/s
elseif fluid == "water"
    nu = getWaterProp(T_K, 'nu'); % expects m^2/s
else
    error('Unknown fluid: use ''air'' or ''water''.');
end

% Compute Reynolds number
Re = V * L / nu;
end


function Pr = getPr(fluid, T_K)
% Pr = getPr(fluid, T_K)
% Returns Prandtl number
% fluid = 'air' or 'water'
% T_K = temp of fluid

if strcmpi(fluid, 'water')
    Pr = getWaterProp(T_K, 'pr');
elseif strcmpi(fluid, 'air')
    Pr = getAirProp(T_K, 'pr');
else
    error('Unknown fluid. Use ''water'' or ''air''.');
end
end

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

%% Forced External

% Flat plate
function [Nu_L, convType] = ExtConvFlatPlate(Re_L, Pr)
% Case: External convection, flat plate
%   Pr >= 0.6
%   T_F used for Pr
%   Selects laminar/mixed based on Re
% Nu_L = ExtConvFlatPlate(Re_L, Pr)
% Inputs:
%   Re_L, Pr
% Outputs:
%   Nu_L -- avg Nusselt
%   convType - laminar, mixed, turbulent

% Get Pr from T_f

if Re_L <= 5*10^5 && Pr >= 0.6
    convType = 'laminar';
    Nu_L = 0.664*Re_L.^(1/2)*Pr.^(1/3);
elseif Re_L >= 5*10^5 && Re_L <= 10^8 && Pr >= 0.6 && Pr <= 60
    convType = 'mixed';
    Nu_L = (0.037*Re_L.^(4/5)-871) .*Pr.^(1/3);
else
    warning('Re_L or Pr outside acceptable range.');
    warning('Make sure to use T_f for Pr.');
    Nu_L = [];
    convType = 'RE_OUTSIDE';
end
end

% Cylinder
function Nu_D = ExtConvCyl(Re_D, Pr)
% Case: External convection, cylinder in cross flow
%   Re_D*Pr >= 0.2
%   T_F used for Pr
% Nu_D = ExtConvCyl(Re_D, Pr)
% Inputs:
%   Re_D, Pr
% Outputs:
%   Nu_D -- avg Nusselt

% Get Pr from T_f

if Re_D*Pr >= 0.2
    Nu_D = 0.3 + (0.62*Re_D.^(1/2).*Pr.^(1/3).* (1+(0.4/Pr).^(2/3)).^(-1/4)).* (1+(Re_D/282000).^(5/8)).^(4/5);
else
    disp('Re_D*Pr outside acceptable range.');
    disp('Make sure to use T_f for Pr.')
    Nu_D = [];
end
end

% Sphere
function Nu_D = ExtConvSphere(Re_D, Pr, mu, mu_s)
% Case: External convection, sphere
%   3.5 <= Re_D <= 7.6e4
%   0.71 <= Pr <= 380
%   1.0 <= mu/mu_s <= 3.2
%   T_F used for Pr
% Nu_D = ExtConvSphere(Re_D, Pr, mu, mu_s)
% Inputs:
%   Re_D, Pr, mu, mu_s
%   mu = dynamic viscosity
%   mu_s = dynamic viscosity at the surface
% Outputs:
%   Nu_D -- avg Nusselt

% Get Pr from T_inf

if Re_D <= 7.6e4 && Re_D >= 3.5 && Pr >= 0.71 && Pr <= 380 && mu/mu_s >= 1 && mu/mu_s <= 3.2
    Nu_D = 2 + (0.4*Re_D.^(1/2) + 0.06*Re_D.^(2/3)).*Pr.^(0.4).*(mu./mu_s).^(1/4);
else
    disp('Re_D, Pr, or mu/mu_s outside acceptable range.');
    disp('Make sure to use T_f for Pr.')
    Nu_D = [];
end
end

% Tube bank
function Nu_D = ExtConvTubeBank(Re_Dmax, Pr, Pr_s, tubeType, N_L, N_T, S_T, S_L)
% Case: External convection, tube bank
%   10 <= Re_D <= 2e6
%   0.7 <= Pr <= 500
%   Use T_bar (need a calc)
% Nu_D = ExtConvTubeBank(Re_Dmax, Pr, tubeType, N_L, N_T)
% Inputs:
%   Re_D, Pr, tubeType, N_L, N_T
%   tubeType = 1 (Aligned) or 2 (Staggered)
%   N_L = longitudinal pitch
%   N_T = transverse pitch
% Outputs:
%   Nu_D -- avg Nusselt

% Calc Re_Dmax
%   Calc Vmax with S_T, D
% Calc T_bar
%   T_bar = mean([T_i T_o])
% Calc dT_lm
%   Need T_s, T_o, T_i
% Calc C1 and C2
%   Switch statement to select C1 and m based on Re_Dmax & tubeType,
%       S_T/S_L
%   Calc C2 based on N_L

disp('****INCOMPLETE');

if Re_Dmax 
    Nu_D = C1*C2*Re_Dmax.^(m).*Pr.^0.36.*(Pr./Pr_s).^(1/4);
else
    disp('Re_D, Pr, or mu/mu_s outside acceptable range.');
    disp('Make sure to use T_f for Pr.')
    Nu_D = [];
end

% Can also output q' = N*h*pi*D*dT_lm
% Can also calc pressure drop delP
end

%% Forced Internal

%% Free Convection

% Vertical plate
function Nu_L = freeConvExtFlatPlateVert(Ra_L, Pr)
% Case: Free convection, external, vertical plate
%   Ra_L >= 10^9  ---- Works for all Ra_L, but another better for laminar
% Nu_L = freeConvExtFlatPlateVert(Ra_L, Pr)
% Inputs:
%   Ra_L, Pr
% Outputs:
%   Nu_L -- avg Nusselt

Nu_L = (0.825+ 0.387*Ra_L.^(1/6) ./(1+ (0.492/Pr).^(9/16) ).^(8/27)) .^2; % verified
end

function Nu_L = freeConvExtFlatPlateVertLam(Ra_L, Pr)
% Case: Free convection, external, vertical plate, laminar
%   Ra_L <= 10^9  ---- Works better for laminar
% Nu_L = freeConvExtVertPlate(Ra_L, Pr)
% Inputs:
%   Ra_L, Pr
% Outputs:
%   Nu_L

Nu_L = 0.68 + 0.670*Ra_L.^(1/4)./ (1+(0.492/Pr).^(9/16)).^(4/9);
end

% Inclined plate
function Nu_L = freeConvExtPlateInc(Ra_L,Pr,theta)
% Case: Free convection, external, inclined plate, cold up or hot down
%   0 <= theta <= 60 (degrees)
%   Use g*cos(theta) instead of g for Ra_L calc **
%       Just use the laminar vert plate
% Nu_L = freeConvExtPlateInc(Ra_L,Pr)
% Inputs:
%   Ra_L, Pr
% Outputs:
%   Nu_L

if theta >= 0 && theta <= 60
    Nu_L = freeConvExtFlatPlateVertLam(Ra_L, Pr); % use g*cos(theta) for Ra'
else
    disp('theta outside the range for an inclined plate. Consider using flat plate.')
end
end

% Horizontal plate
function Nu_L = freeConvExtPlateHorizHotUpper(Ra_L,Pr)
% Case: Free convection, external, horizontal plate, hot surf upper, looking @
%           top surface
%   10^4 <= Ra_L <= 10^11
%   Pr >= 0.7 (if Ra_L <= 0.7)
% Nu_L = freeConvExtPlateHorizHotUpLamTop(Ra_L,Pr)
% Inputs:
%   Ra_L, Pr
% Outputs:
%   Nu_L

if Ra_L >= 10^4 && Ra_L <= 10^7 && Pr >= 0.7
    Nu_L = 0.54*Ra_L.^(1/4);
elseif Ra_L >= 10^7 && Ra_L <= 10^11
    Nu_L = 0.15*Ra_L.^(1/3);
else
    disp('Invalid combination of Ra_L and Pr');
    Nu_L = [];
end
end

function Nu_L = freeConvExtPlateHorizHotLower(Ra_L,Pr)
% Case: Free convection, external, horizontal plate, hot surf down, looking @
%           bottom surface (or upper surface of cold plate)
%   10^4 <= Ra_L <= 10^9
%   Pr Nu_L = freeConvExtPlateHorizColdUp(Ra_L,Pr)
% Inputs:
%   Ra_L, Pr
% Outputs:
%   Nu_L

if Ra_L >= 10^4 && Ra_L <= 10^9 && Pr >= 0.7
    Nu_L = 0.52*Ra_L.^(1/5);
else
    disp('Invalid combination of Ra_L and Pr');
    Nu_L = [];
end
end

function Nu_D = freeConvExtHorizCyl(Ra_D,Pr)
% Case: Free convection, external, horizontal cylinder
%   Ra_D <= 10^12
%   Nu_L = freeConvExtHorizCyl(Ra_L,Pr)
% Inputs:
%   Ra_D, Pr
% Outputs:
%   Nu_L

if Ra_D <= 10^12
    Nu_D = (0.60 + 0.387*Ra_D.^(1/6) ./(1+(0.559/Pr).^(9/16)).^(8/27)).^2;
else
    disp('Invalid combination of Ra_D and Pr');
    Nu_D = [];
end
end

function Nu_D = freeConvExtSphere(Ra_D,Pr)
% Case: Free convection, external, sphere
%   Ra_D <= 10^11
%   Pr >= 0.7
%   Nu_D = freeConvExtSphere(Ra_D,Pr)
% Inputs:
%   Ra_D, Pr
% Outputs:
%   Nu_L

if Ra_D <= 10^12 && Pr >= 0.7
    Nu_D = 2 + 0.589*Ra_D.^(1/4) ./(1+(0.469/Pr).^(9/16)).^(4/9);
else
    disp('Invalid combination of Ra_D and Pr');
    Nu_D = [];
end
end

%% Interpolation Ftn
function val = getAirProp(T_K, prop)
% getAirProp: interpolates property tables for air
% val = getAirProp(T, prop)
% Inputs:
%   T_K = temperature (K)
%   prop = string property name to interpolate for
%           rho, cp, mu, k, alpha, Pr
% Output:
%   val = interpolated value

if T_K < 100 || T_K > 3000
    disp('Temperature outside acceptable range.');
    val = [];
    return;
end

% Thermophysical Properties of Air at Atmospheric Pressure
% (Table A.4)

T_air = [100 150 200 250 300 350 400 450 500 550 600 650 700 750 800 ...
     850 900 950 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900 ...
     2000 2100 2200 2300 2400 2500 3000];            % Temperature (K)

rho_air = [3.5562 2.3364 1.7458 1.3947 1.1614 0.9950 0.8711 0.7740 0.6964 ...
       0.6329 0.5804 0.5356 0.4975 0.4643 0.4354 0.4097 0.3868 0.3666 ...
       0.3482 0.3166 0.2902 0.2679 0.2488 0.2322 0.2177 0.2049 0.1935 ...
       0.1833 0.1741 0.1658 0.1582 0.1513 0.1448 0.1389 0.1135];        % Density (kg/m^3)

cp_air = [1.032 1.012 1.007 1.006 1.007 1.009 1.014 1.021 1.030 1.040 ...
      1.051 1.063 1.075 1.087 1.099 1.110 1.121 1.131 1.141 1.159 ...
      1.175 1.189 1.207 1.230 1.248 1.267 1.286 1.307 1.337 1.372 ...
      1.417 1.478 1.558 1.665 2.726];                                   % Specific heat (kJ/kg·K)

mu_air = [71.1 103.4 132.5 159.6 184.6 208.2 230.1 250.7 270.1 288.4 305.8 ...
      322.5 338.8 354.6 369.8 384.3 398.1 411.3 424.4 449.0 473.0 496.0 ...
      530.0 557.0 584.0 611.0 637.0 663.0 689.0 715.0 740.0 766.0 792.0 ...
      818.0 955.0] * 1e-7;                                              % Dynamic viscosity (N·s/m^2)

nu_air = [2.00 4.426 7.590 11.44 15.89 20.92 26.41 32.39 38.79 45.57 52.69 ...
      60.21 68.10 76.37 84.93 93.80 102.9 112.2 121.9 141.8 162.9 185.1 ...
      213.0 240.0 268.0 298.0 329.0 362.0 396.0 431.0 468.0 506.0 547.0 ...
      589.0 841.0] * 1e-6;                                              % Kinematic viscosity (m^2/s)

k_air = [9.34 13.8 18.1 22.3 26.3 30.0 33.8 37.3 40.7 43.9 46.9 49.7 52.4 ...
     54.9 57.3 59.6 62.0 64.3 66.7 71.5 76.3 82.0 91.0 100.0 106.0 ...
     113.0 120.0 128.0 137.0 147.0 160.0 175.0 196.0 222.0 486.0] * 1e-3; % Thermal conductivity (W/m·K)

alpha_air = [2.54 5.84 10.3 15.9 22.5 29.9 38.3 47.2 56.7 66.7 76.9 87.3 ...
          98.0 109.0 120.0 130.0 143.0 155.0 168.0 195.0 224.0 257.0 ...
          303.0 350.0 390.0 435.0 482.0 534.0 589.0 646.0 714.0 783.0 ...
          869.0 960.0 1570.0] * 1e-6;                                   % Thermal diffusivity (m^2/s)

Pr_air = [0.786 0.758 0.737 0.720 0.707 0.700 0.690 0.686 0.684 0.683 ...
      0.685 0.690 0.695 0.702 0.709 0.716 0.720 0.723 0.726 0.728 ...
      0.728 0.719 0.703 0.685 0.688 0.685 0.683 0.677 0.672 0.667 ...
      0.655 0.647 0.630 0.613 0.536];                                  % Prandtl number (-)

switch lower(prop)
    case 'rho'
        % pchip interpolation (see Ch 18.5)
        % density (kg/m^3)
        val = interp1(T_air, rho_air, T_K, 'pchip');

    case 'mu'
        % Dynamic viscosity (Pa*s)
        val = interp1(T_air, mu_air, T_K, 'pchip');

    case 'cp'
        % Specific heat (kJ/kg*K)
        val = interp1(T_air, cp_air, T_K, 'pchip');

    case 'k'
        % Thermal conductivity (W/m*K)
        val = interp1(T_air, k_air, T_K, 'pchip');

    case 'nu'
        % Kinematic viscosity (m^2/s)
        val = interp1(T_air, nu_air, T_K, 'pchip');

    case 'alpha'
        % Thermal diffusivity (m^2/s)
        val = interp1(T_air, alpha_air, T_K, 'pchip');

    case 'pr'
        % Prandtl number (-)
        val = interp1(T_air, Pr_air, T_K, 'pchip');

    otherwise
        error('Property not recognized. Use: rho, mu, cp, or k.');
end
end

function val = getWaterProp(T_K, prop)
% getAirProp: interpolates property tables for air
% val = getWaterProp(T, prop)
% Inputs:
%   T_K = temperature (K)
%   prop = string property name to interpolate for
%           rho, cp, mu, k, alpha, Pr
% Output:
%   val = interpolated value

%Cesar
% mirror the getAirProps ftn
end
