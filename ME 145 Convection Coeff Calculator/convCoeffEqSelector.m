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

                Re_L = getRe();
                disp('***REMOVE*** Make sure its T_f for Pr');
                Pr = getPr(); % make sure it's T_f

                [Nu_L, convType] = ExtConvFlatPlate(Re_L, Pr);
            case 2
                clc;
                disp('External conv, cylinder');

                Re_D = getRe();
                disp('***REMOVE*** Make sure its T_f for Pr');
                Pr = getPr(); % make sure it's T_f


                Nu_D = ExtConvCyl(Re_D, Pr);
            case 3
                clc;
                disp('External conv, sphere');

                Re_D = getRe();
                disp('***REMOVE*** Make sure its T_inf for Pr');
                Pr = getPr(); % make sure it's T_inf
                mu = input('Dynamic viscosity (mu): ');
                mu_s = input('Dynamic viscosity @ surface (mu_s):');

                Nu_D = ExtConvSphere(Re_D, Pr, mu, mu_s);
            case 4
                clc;
                disp('External conv, Bank of tubes');

                Re_Dmax = getRe();
                disp('***REMOVE*** Make sure its T_bar for Pr');
                Pr = getPr(); % make sure it's T_bar
                disp('***REMOVE*** Make sure its T_s for Pr_s')
                Pr = getPr(); % make sure it's T_bar

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
        disp('  1. External free convection');
        disp('  2. Internal free convection');
        freeCase = input('Free convection case (1–2): ');

        switch freeCase
            case 1  % External
                clc;
                disp('Select an external free convection:');
                disp('  1. Flat plate');
                disp('  2. Horizontal cylinder');
                disp('  3. Sphere');
                extFree = input('External free conv case (1–3): ');

                switch extFree
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
                                disp('Free conv, external, vert flat plate');

                                % Get Ra_L, Pr
                                convType = 'Length'; % change
                                Ra_L = getRa();
                                Pr = getPr();
                                
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
                                    disp('Free conv, external, flat plate, hot surface up / cold surface down');

                                    Ra_L = getRa();
                                    Pr = getPr();
                                    Nu_L = freeConvExtPlateHorizHotUpper(Ra_L,Pr);
                                elseif hotcold == 2
                                    clc;
                                    disp('Free conv, external, flat plate, cold surface up / hot surface down');

                                    Ra_L = getRa();
                                    Pr = getPr();
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
                                    disp('Free conv, external, inclined plate, old surface up / hot surface down');

                                    theta = input('Inclination angle theta (deg):');
                                    %   Ra calc'd with g*cos(theta) *****
                                    Ra_D = getRa();
                                    Pr = getPr();

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
                        disp('Free conv, external, horizontal cylinder');
                        
                        Ra_D = getRa();
                        Pr = getPr();
                        Nu_D = freeConvExtHorizCyl(Ra_D,Pr);
                    case 3
                        clc;
                        disp('Free conv, external, sphere');

                        Ra_D = getRa();
                        Pr = getPr();
                        
                        Nu_D = freeConvExtSphere(Ra_D,Pr);
                    otherwise
                        clc;
                        disp('Invalid selection.');
                end

            case 2  % Internal free convection
                clc;
                disp('Select an internal free convection geometry:');
                disp('  1. Rectangular cavity');
                disp('  2. Concentric cylinders');
                disp('  3. Concentric spheres');
                intFreeGeom = input('Geometry (1–3): ');

                switch intFreeGeom
                    case 1
                        clc;
                        disp('Select a rectangular cavity orientation:');
                        disp('  1. Vertical');
                        disp('  2. Horizontal');
                        disp('  3. Inclined');
                        cavity = input('Orientation (1–3): ');
                        if cavity == 1
                            clc;
                            disp('Free convect, internal, rectangular, vertical rectangular cavity');
                        elseif cavity == 2
                            clc;
                            disp('Free convect, internal, rectangular, horizontal rectangular cavity');
                        elseif cavity == 3
                            clc;
                            disp('Free convect, internal, rectangular, inclined rectangular cavity');
                        else
                            clc;
                            disp('Invalid selection.');
                        end
                    case 2
                        clc;
                        disp('Free convect, internal, concentric cylinders');
                    case 3
                        clc;
                        disp('Free convect, internal, concentric spheres');
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
disp('   We love you Troy, so please give us the A :)');
disp('----------------------------------');


%% Functions ============================================================

% Get user input
function Re = getRe()
% Gets Reynolds number
Re = input('Enter Reynolds: ');
end

function Pr = getPr()
% Pr = getPr()
% Returns Prandtl number
Pr = input('Enter Prandtl: ');
end

function Ra = getRa()
% Gets Reynolds number
Ra = input('Enter Rayleigh: ');
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
    disp('Re_L or Pr outside acceptable range.');
    disp('Make sure to use T_f for Pr.')
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

if Re_D.*Pr >= 0.2
    Nu_D = 0.3 + (0.62*Re_D.^(1/2).*Pr.^(1/3).* (1+(0.4/Pr).^(2/3)).^(1/4)).* (1+(Re_D/282000).^(5/8)).^(4/5);
else
    disp('Re_D*Pr outside acceptable range.');
    disp('Make sure to use T_f for Pr.')
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
%   T_bar = mean(T_i, T_o)
% Calc dT_lm
%   Need T_s, T_o, T_i
% Calc C1 and C2
%   Switch statement to select C1 and m based on Re_Dmax & tubeType,
%       S_T/S_L
%   Calc C2 based on N_L

if Re_Dmax 
    Nu_D = C1*C2*Re_Dmax.^(m).*Pr.^0.36.*(Pr./Pr_s).^(1/4);
else
    disp('Re_D, Pr, or mu/mu_s outside acceptable range.');
    disp('Make sure to use T_f for Pr.')
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