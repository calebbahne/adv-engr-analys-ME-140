%% Convection Coefficient Calculator
% Purpose: Calculate the avg convection coefficient h for a given condition
%   Determines the appropriate case, then does calcs from there.

% To Do:
%   - Characteristic length calculator for Ra
%   - Print Ra or Re depending on case in show final results
%   - Tube bank verify 
%   - Water properties
%   - Other missing ones
%   - Test and verify, put 'verified' on what's tested

clc; clear;
Re_L = NaN; Re_D = NaN; Ra_L = NaN; Ra_D = NaN; L = NaN; D = NaN;

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

                T_inf = input('Enter the Upstream Fluid Temperature (C): ')+273.15;
                T_s = input('Enter the Surface Temperature (C): ')+273.15;
                T_f = (T_s+T_inf)/2; % Tfilm
                clc;

                [Re_L L] = getRe(T_f, fluid, 'Re_L');
                Pr = getPr(fluid, T_f); % make sure it's T_f

                [Nu_L, convType] = ExtConvFlatPlate(Re_L, Pr);
            case 2
                clc;
                disp('External conv, cylinder');

                T_inf = input('Enter the Upstream Fluid Temperature (C): ')+273.15;
                T_s = input('Enter the Surface Temperature (C): ')+273.15;
                T_f = (T_s+T_inf)/2; % Tfilm
                clc;

                [Re_D D]  = getRe(T_f, fluid, 'Re_D');
                Pr = getPr(fluid, T_f); % make sure it's T_f


                Nu_D = ExtConvCyl(Re_D, Pr);
            case 3
                clc;
                disp('External conv, sphere');

                T_inf = input('Enter the Upstream Fluid Temperature (C): ')+273.15;
                T_s = input('Enter the Surface Temperature (C): ')+273.15;
                T_f = (T_s+T_inf)/2; % Tfilm
                clc;

                [Re_D D]  = getRe(T_f, fluid, 'Re_D');
                Pr = getPr(fluid, T_f); % make sure it's T_inf
                mu = getFluidProp(fluid, T_f, 'mu');
                mu_s = getFluidProp(fluid, T_s, 'mu');

                Nu_D = ExtConvSphere(Re_D, Pr, mu, mu_s);
            case 4
                clc;
                disp('External conv, Bank of tubes');

                [Nu_D, Re_Dmax, V_max, T_i, T_o, T_s, D, V, S_T, S_L, ...
                    N_L, N_T, C1, C2, m, Pr_s, tubeType, q_p, DT_lm, iter, tol, converged, T_calc_hist, T_update_hist] = handleTubeBank(fluid);
                if converged
                    fprintf('\nTube Bank calcs converged in %d iterations within %.2f%% tolerance.\n', iter, tol);
                else
                    warndlg('Tube Bank calcs did not converge.');
                end

                Re_D = Re_Dmax; % for display at end??
                T_f = (T_i+T_o)/2; % use T_i and T_o for Nu calcs
            otherwise
                disp('Invalid selection.');
        end

    case 2  % Internal Flow
        clc;
        disp('Select an internal flow:');
        disp('  1. Circular tube');
        disp('  2. Noncircular tube');
        intConvCase = input('Internal flow case (1–2): ');

        switch intConvCase
            case 1 % Circular tube
                clc;
                disp('Int conv, circular tube');

                [T_s, T_mi, T_mo, T_m, T_f, SurfCondit, assumedVar, tol] = getTubeVals();
                [Re_D,D] = getRe(T_m, fluid, 'Re_D');
                L = input('Enter the tube length, L (m): '); 
                [Nu_D, convType] = IntTube(fluid, Re_D, T_s, T_m, L, D, SurfCondit);
               
            case 2 % Noncircular tube
                clc;

                [T_s, T_mi, T_mo, T_m, T_f, SurfCondit, assumedVar, tol] = getTubeVals();

                [Re_D,D_h, geom,a, b] = getReHydraul(T_m, fluid); % handles shape, gets hydraulic diam
                L = input('Enter the tube length, L (m): '); 
                [Nu_D, convType] = IntTubeNonCirc(fluid, Re_D, T_s, T_m, L, D_h, SurfCondit, geom, a, b);
                D = D_h; % saves time later
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
                                [Ra_L, T_s, T_inf, T_f, L] = getRa(fluid, convType);
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

                                    [Ra_L, T_s, T_inf, T_f, L] = getRa(fluid, convType);
                                    Pr = getPr(fluid, T_f);
                                    Nu_L = freeConvExtPlateHorizHotUpper(Ra_L,Pr);
                                elseif hotcold == 2
                                    clc;
                                    disp('Free conv, immersed, flat plate, cold surface up / hot surface down');
                                    convType = 'Ra_L';

                                    [Ra_L, T_s, T_inf, T_f, L] = getRa(fluid, convType);
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
                                    disp('Free conv, immersed, inclined plate, cold surface up / hot surface down');

                                    theta = input('Inclination angle theta (deg):');
                                    convType = 'Ra_L';
                                    %   Ra calc'd with g*cosd(theta)
                                    clc;
                                    [Ra_L, T_s, T_inf, T_f, L] = getRa(fluid, convType);
                                    Ra_L = Ra_L*cosd(theta);
                                    Pr = getPr(fluid, T_f);

                                    Nu_L = freeConvExtPlateInc(Ra_L,Pr,theta);
                                else
                                    clc;
                                    warndlg('Not currently supported (not in textbook)');
                                    % this one's not in table 9.3
                                end
                            otherwise
                                disp('Invalid selection.');
                        end

                    case 2
                        clc;
                        disp('Free conv, immersed, horizontal cylinder');
                        convType = 'Ra_D';
                        
                        [Ra_D, T_s, T_inf, T_f, D] = getRa(fluid, convType);
                        T_f = mean([T_s T_inf]);
                        Pr = getPr(fluid, T_f);
                        Nu_D = freeConvExtHorizCyl(Ra_D,Pr);
                    case 3
                        clc;
                        disp('Free conv, immersed, sphere');
                        convType = 'Ra_D';
                        
                        [Ra_D, T_s, T_inf, T_f, D] = getRa(fluid, convType);
                        T_f = mean([T_s T_inf]);
                     
                        Pr = getPr(fluid, T_f);
                        
                        Nu_D = freeConvExtSphere(Ra_D,Pr);
                    otherwise
                        clc;
                        disp('Invalid selection.');
                end

            case 2  % Enclosed free convection
                clc;
                disp('Select an enclosed free convection geometry:');
                disp('  1. Rectangular cavity');
                disp('  2. Other (not yet supported)');
                encFreeGeom = input('Geometry (1–2): ');

                switch encFreeGeom
                    case 1
                        clc;
                        disp('Select a rectangular cavity orientation:');
                        disp('  1. Vertical');
                        disp('  2. Horizontal');
                        cavity = input('Orientation (1–2): ');
                        if cavity == 1
                            clc;
                            disp('Free convect, enclosed, rectangular, vertical rectangular cavity');

                            H = input('Enter the height of the enclosure (H): ');
                            L = input('Enter the length of the enclosure (L): ');
                            T1 = input('Enter the hotter wall temperature, T1 (C): ')+273.15; % K
                            T2 = input('Enter the colder wall temperature, T2 (C): ')+273.15; % K
                            clc;
                            [Ra_L, T_s, T_inf, T_f, L] = getRa(fluid, 'Ra_L', T1, T2, L); 
                            Pr = getFluidProp(fluid, T_f, 'Pr'); 
                            Nu_L = freeConvEncRectVertCav(Ra_L,Pr,H,L);
                        elseif cavity == 2
                            clc;
                            disp('Free convect, enclosed, rectangular, horizontal rectangular cavity');

                            H = input('Enter the height of the enclosure (H): ');
                            L = input('Enter the length of the enclosure (L): ');
                            T1 = input('Enter the hotter wall temperature, T1 (C): ')+273.15; % K
                            T2 = input('Enter the colder wall temperature, T2 (C): ')+273.15; % K
                            clc;
                            [Ra_L, T_s, T_inf, T_f, L] = getRa(fluid, 'Ra_L', T1, T2, L); 
                            Pr = getFluidProp(fluid, T_f, 'Pr'); 
                            Nu_L = freeConvEncRectCavHorizHeatBelow(Ra_L,Pr);
                        else
                            clc;
                            warndlg('Invalid selection.');
                        end
                    case 2
                        warndlg('Nonrectangular cavities not yet supported.');
                end
            otherwise
                clc;
                warndlg('Invalid selection.');
        end

    otherwise
        warndlg('Invalid main selection.');
end

%% Display results
disp('----------------------------------');
disp('Convection coefficient calculation complete');
disp('----------------------------------');

% Shows the final values based on what exists
% Lets the user say if they want all values or just a few

if exist('Nu_L', 'var') && exist('L', 'var')
    Nu = Nu_L;
    charDim = L;
    Re_end = Re_L;
    Ra_end = Ra_L
elseif exist('Nu_D', 'var') && exist('D', 'var')
    Nu = Nu_D;
    charDim = D;
    Re_end = Re_D;
    Ra_end = Ra_D;
else
    warndlg('Error: Neither (Nu_L, L) nor (Nu_D, D) found.');
end

if ~exist('T_f', 'var')
    T_f = T_m; % happens if we have T_K and not T_f
end

k = getFluidProp(fluid, T_f, 'k');

% Calculate convection coefficient (this is the whole point)
h = Nu * k / charDim;

clc;

rho   = getFluidProp(fluid, T_f, 'rho');
cp    = getFluidProp(fluid, T_f, 'cp');
mu    = getFluidProp(fluid, T_f, 'mu');
nu    = getFluidProp(fluid, T_f, 'nu');
Pr    = getFluidProp(fluid, T_f, 'Pr');

% Display property table
fprintf('\n\n--- %s Properties at %.2f K ---\n', upper(fluid), T_f);
fprintf('Density (rho):             %.4f kg/m³\n', rho);
fprintf('Specific Heat (cp):        %.4e J/kg·K\n', cp);
fprintf('Dynamic Viscosity (mu):    %.4e Pa·s\n', mu);
fprintf('Kinematic Viscosity (nu):  %.4e Pa·s\n', nu);
fprintf('Thermal Conductivity (k):  %.4f W/m·K\n', k);
fprintf('Prandtl Number (Pr):       %.4f\n', Pr);

fprintf('\n----------------------------------------\n');

% Display main result
fprintf('\n=== Convection Coefficient Results ===\n');
fprintf('Fluid:                             %s\n', fluid);
fprintf('Film Temperature (T_f):            %.2f C\n', T_f-273.15);
fprintf('Characteristic Dimension (L or D): %.4f m\n', charDim);
if ~isnan(Re_end)
    fprintf('Reynolds Number (Re):              %d \n',Re_end);
end
if ~isnan(Ra_end)
    fprintf('Rayleigh Number (Ra):              %d \n',Ra_end);
end
fprintf('Nusselt Number (Nu):               %.4f\n', Nu);
fprintf('Convection Coefficient (h):        %.4f W/m²·K\n', h);

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
        warndlg('Invalid entry. Air is autoselected');
end
end

function [Re,charDim] = getRe(T_K, fluid, re_type)
% getRe:  Calculate Reynolds number using nu
%   Re = getRe(T_K, fluid, re_type)
%
%   T_K = temp in Kelvin
%   fluid = 'air' or 'water'
%   re_type = 'Re_D'/'D' (use diameter) or 'Re_L'/'L' (use length)

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
disp('   2. Mass flow rate (mdot, kg/s)');
disp('   3. Direct input (Re)');
ReChoice = input('Re input method: ');
clc;

if ReChoice == 3
    Re = input('Enter Re: ');
    charDim = input(sprintf('Enter characteristic %s (m): ', charName));
    return
elseif ReChoice == 2
    % Compute from mass flow rate mdot (kg/s)
    mdot = input('Enter mass flow rate mdot (kg/s): ');
    mu = getFluidProp(fluid, T_K, 'mu');
    if contains(lower(charName),'diameter')
        charDim = input(sprintf('Enter characteristic %s (m): ', charName));
        D = charDim;
        Re = (4 * mdot) / (pi * D * mu);
        return;
    else
        % For non-diameter , ask for cross-sectional area A (m^2)
        A = input('Enter cross-sectional area A (m^2): ');
        Re = mdot * charDim / (A * mu);
        return;
    end
end
clc;

% Otherwise compute from V and characteristic length
V = input('Enter velocity V (m/s): ');
charDim = input(sprintf('Enter characteristic %s (m): ', charName));

clc;

nu = getFluidProp(fluid, T_K, 'nu');

% Compute Reynolds number
Re = V * charDim / nu;
end

function [Re_D, D_h, geom, a, b] = getReHydraul(T_m, fluid)
% getReHydraul: gets the Reynolds number for hydraulic diameter
% [Re_D, D_h, geom] = getReHydraul(T_m, fluid)
% Input:
%   T_m = mean temp (input and output means)
%   fluid 
% Output
%   Re_D
%   D_h = hydraulic diameter
%   geom = 1 if rect/square, 2 if triangle
    clc;
disp('Select a noncircular cross-section:');
disp('  1. Rectangle / square');
disp('  2. Triangle');
geom = input('Select geometry (1-2): ');
clc;

switch geom
    case 1 % rectangle/square
        fprintf('Rectangle duct selected.\n');
        a = input('Enter shorter side length, a (m): ');
        b = input('Enter longer side length, b (m): ');
        % Ensure positive
        assert(a > 0 && b > 0, 'Sides must be positive.');
        % Cross-sectional area and perimeter
        A_c = a * b;
        P = 2*(a + b);
        b_over_a = b / a;
    case 2 % triangle
        fprintf('Triangular duct selected.\n');
        base = input('Enter base (m): ');
        height = input('Enter height (m): ');
        use_perim = input('Do you know the wetted perimeter (y/n)? ', 's');
        if lower(use_perim) == 'y'
            P = input('Enter wetted perimeter P (m): ');
            A_c = 0.5 * base * height;
        else
            % assume an isosceles with two equal sides: compute sides from base & height
            side = sqrt((base/2)^2 + height^2);
            P = base + 2*side;
            A_c = 0.5 * base * height;
        end
        a = NaN; b = NaN;
    otherwise
        error('Invalid geometry selection.');
end

% Hydraulic diameter
D_h = 4*A_c/P;
clc;

disp('Select Re input method for this hydraulic diameter:');
disp('  1. Velocity (V)');
disp('  2. Mass flow rate (mdot)');
disp('  3. Direct input (Re)');
ReChoice = input('Re input method (1-3): ');
clc;

% Get fluid properties needed
nu = getFluidProp(fluid, T_m, 'nu');   % kinematic viscosity m^2/s
mu  = getFluidProp(fluid, T_m, 'mu');  % dynamic viscosity Pa*s
rho = getFluidProp(fluid, T_m, 'rho'); % density kg/m^3

if ReChoice == 3
    Re_D = input('Enter Re based on D_h: ');
else
    if ReChoice == 1
        V = input('Enter velocity V (m/s): ');
        Re_D = V * D_h / nu;
    elseif ReChoice == 2
        mdot = input('Enter mass flow rate mdot (kg/s): ');
        um = mdot / (rho * A_c);
        Re_D = um * D_h / nu;
    else
        error('Invalid Re input selection.');
    end
end
end


function Pr = getPr(fluid, T_K)
% Pr = getPr(fluid, T_K)
% Returns Prandtl number
% fluid = 'air' or 'water'
% T_K = temp of fluid

Pr = getFluidProp(fluid, T_K, 'pr');
end

function [Ra, T_s, T_inf, T_f, L] = getRa(fluid, ra_type, T1, T2, L)
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

if nargin < 4 || isempty(T1) || isempty(T2)
    % Prompt user to select input method
    disp('Select Ra input method:');
    fprintf('   1. Compute from characteristic (%s)\n', charName);
    disp('   2. Direct input (Ra)');
    RaChoice = input('Ra input method: ');
    clc;

    if RaChoice == 2
        Ra = input('Enter Rayleigh number (Ra): ');
        return
    end
    clc;
    
    % Prompt only if temps not passed in
    T_s  = input('Enter surface temperature, T_s (C): ')+273.15;
    T_inf = input('Enter ambient temperature, T_inf (C): ')+273.15;
else
    % Use the provided function inputs
    T_s  = T1;
    T_inf = T2;
end

if nargin < 5 || isempty(L)
    L = input(sprintf('Enter characteristic %s (m): ', charName));
end
clc;

T_f = mean([T_s T_inf]);
% alpha = getFluidProp(fluid, T_f, 'alpha';
nu = getFluidProp(fluid, T_f, 'nu');
Pr = getFluidProp(fluid, T_f, 'pr');

% Constants
g = 9.81;                  % grav accel (m/s^2)

dT = 0.1;

if strcmpi(fluid,'water') % can't use ideal gas approx
        rho_m = getFluidProp(fluid, T_f, 'rho');         
        rho_plus = getFluidProp(fluid, T_f + dT, 'rho'); 
        rho_minus = getFluidProp(fluid, T_f - dT, 'rho');
        % I found out we needed to do this from AI, not in the textbook

    % derivative 
    drho_dT = (rho_plus - rho_minus) / (2*dT);
    beta = - (1 ./ rho_m) .* drho_dT; % vol expansion coeff

    % sanity check: if beta is negative or zero, warn and clamp to reasonable value
    if beta <= 0
        warndlg('beta not positive, falling back to lit val');
        beta = 4.6e-4; % at 50 C
    end
else
    % Ideal gas approx, doesn't work for water
    beta = 1 / T_f;
end

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
    Nu_L = NaN;
    convType = 'RE_OUTSIDE';
    warndlg(sprintf('Re_L or Pr outside acceptable range.\nRe_L = %.2e\nPr = %.3f', Re_L, Pr), 'Warning'); 
    % I found warndlg online: it pops up a GUI so we can find after fails
    warndlg('Make sure to use T_f for Pr.');
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
    warndlg('Re_D*Pr outside acceptable range.');
    warndlg('Make sure to use T_f for Pr.')
    Nu_D = NaN;
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

if Re_D <= 7.6e4 && Re_D >= 3.5 && Pr >= 0.71 && Pr <= 380 && mu/mu_s >= 1 && mu/mu_s <= 3.2 % expanded range a little
    Nu_D = 2 + (0.4*Re_D.^(1/2) + 0.06*Re_D.^(2/3)).*Pr.^(0.4).*(mu./mu_s).^(1/4);
else
    warndlg('Re_D, Pr, or mu/mu_s outside acceptable range.');
    Nu_D = 2 + (0.4*Re_D.^(1/2) + 0.06*Re_D.^(2/3)).*Pr.^(0.4).*(mu./mu_s).^(1/4);
end
end

%% Forced external, tube 
function [Nu_D, Re_Dmax, V_max, T_i, T_o, T_s, D, V, S_T, S_L, ...
          N_L, N_T, C1, C2, m, Pr_s, tubeType, q_p, DT_lm, iter, tol, converged, T_calc_hist, T_update_hist] = handleTubeBank(fluid)
% handleTubeBank - does all the calcs for tube bank, including iteration,
% finds log mean temp and q per unit length too

[T_i, T_o, T_s, D, V, S_T, S_L, N_L, N_T, tubeType, assumedVar, tol] = getTubeBankInputs();
T_f = mean([T_i, T_o]);
rho = getFluidProp(fluid, T_f, 'rho');
cp  = getFluidProp(fluid, T_f, 'cp');

maxIter = 50; iter = 0; converged = false;
T_calc_hist = [];     % stores T_o_calc or T_i_calc
T_update_hist = [];   % stores updated T_o or T_i

while ~converged && iter < maxIter
    iter = iter + 1;

    % Re and Nu
    [Re_Dmax, V_max] = getReTubeBank(fluid, D, V, T_i, T_o, S_T, S_L, tubeType);
    [Nu_D, C1, C2, m, Pr_s] = ExtConvTubeBank(fluid, Re_Dmax, T_i, T_o, T_s, S_T, S_L, N_L, tubeType);
    
    T_f = mean([T_i T_o]);

    % h
    k = getFluidProp(fluid, T_f, 'k');
    h = Nu_D * k / D;
    N = N_T*N_L;

    % Props
    rho = getFluidProp(fluid, T_f, 'rho');
    cp  = getFluidProp(fluid, T_f, 'cp');

    % Calculate the temperature ratio
    expo = exp(-pi*D*N*h / (rho*V*N_T*S_T*cp)); % exponent stuff

    switch assumedVar % T_i, T_o
        % case 1 = none assumed
        case 2  % assumed T_o
            T_o_calc = T_s - (T_s - T_i) * expo

            % Save iteration data
            T_calc_hist(end+1)   = T_o_calc-273.15;
            T_update_hist(end+1) = T_o-273.15;

            err = abs((T_o_calc - T_o) / T_o) * 100;
            if err <= tol
                converged = true;
            else
                T_o = T_o + 0.5 * (T_o_calc - T_o)
            end

        case 3  % assumed T_i
            T_i_calc = T_s - (T_s - T_o) / expo;

            % Save iteration data
            T_calc_hist(end+1)   = T_i_calc-273.15;
            T_update_hist(end+1) = T_i-273.15;

            err = abs((T_i_calc - T_i) / T_i) * 100;
            if err <= tol
                converged = true;
            else
                T_i = T_i + 0.5 * (T_i_calc - T_i)
            end

        otherwise
            converged = true;  % no assumed variable
    end
end

% After convergence
DT_lm = ((T_s-T_i)-(T_s-T_o))/log((T_s-T_i)/(T_s-T_o)); % log mean temp change
q_p = N_T*N_L*(h*pi*D*DT_lm); % heat rate per unit length

end

function [T_i, T_o, T_s, D, V, S_T, S_L, N_L, N_T, tubeType, assumedVar, tol] = getTubeBankInputs()
% getTubeBankInputs: get stuff for tube bank calcs
% [T_i, T_o, T_s, D, S_T, S_L, tubeType] = getTubeBankInputs()
% Outputs:
%   T_i, T_o, T_s = inlet, outlet, and surface temperatures (K)
%   D = tube diameter (m)
%   S_T, S_L = transverse and longitudinal pitches (m)
%   tubeType = 1 = Aligned, 2 = Staggered

% Temperature inputs
disp('Temperature:');
T_i = input('  Enter fluid inlet temperature T_i (C): ');
T_o = input('  Enter fluid outlet temperature T_o (C): ');
T_s = input('  Enter tube surface temperature T_s (C): ');
clc;

T_i = T_i + 273.15; % convert to kelvin
T_o = T_o + 273.15;
T_s = T_s + 273.15;

disp('Was any temperature assumed?');
disp('  1. None');
disp('  2. Outlet temperature (T_o)');
disp('  3. Inlet temperature (T_i)');
assumedVar = input('Assumption: ');
if assumedVar ~= 1
    tol = input('Enter tolerance for iteration (percent): ');
else
    tol = NaN;
end
clc;

% Tube arrangement type
disp('Select tube arrangement type:');
disp('  1. Aligned');
disp('  2. Staggered');
tubeType = input('Arrangement: ');
clc;

% Geometry
disp('Geometry:');
D = input('  Enter tube outer diameter D (m): ');
S_T = input('  Enter transverse pitch S_T (m): ');
S_L = input('  Enter longitudinal pitch S_L (m): ');
N_L = input('  Enter the rows of tubes (N_L): ');
N_T = input('  Enter the number of tubes per row: (N_T): ');
clc;
V = input('Enter free stream velocity (m/s): ');
clc;
end

function [Re_Dmax, V_max] = getReTubeBank(fluid, D, V, T_i, T_o, S_T, S_L, tubeType)
% getReTubeBank: calculate Reynolds for tube bank
%   [Re_Dmax, V_max] = getReTubeBank(fluid, D, V, T_i, T_o, T_s, S_T, S_L, tubeType
% Inputs:
%   fluid = 'air' or 'water'
%   D = diam (m)
%   V = free stream velo (m/s)
%   T_i = inlet temp (K)
%   T_o = outlet temp (K) (may be assumed?)
%   T_s = surface temp (K)
%   S_T = transverse pitch (m)
%   S_L = longitudinal pitch (m)
%   tubeType = 1 (aligned) or 2 (staggered)
% Outputs:
%   Re_Dmax = Reynolds number (based on max velo thru tube bank)
%   V_max = max velo between tubes (m/s)
%
% Notes:
%   Uses fluid properties at film temperature (T_f = (T_i + T_o) / 2)

% just like getRe
fluid = lower(string(fluid));

% We've gotta use film temp (not log mean)
T_f = 0.5 * (T_i + T_o);  

nu = getFluidProp(fluid, T_f, 'nu');

S_D = (S_L^2+(S_T/2)^2)^0.5;

% Get max velo
if tubeType == 1 % aligned
    V_max = V*S_T/(S_T-D);
elseif tubeType == 2 && S_D < (S_T+D)/2% staggered, A1
    V_max = V*S_T/(2*(S_T-D));
elseif tubeType == 2 && 2*(S_D-D) >= (S_T-D)/2% staggered, A2
    V_max = V*S_T/(S_T-D);
else
    warndlg('Unknown tube type (1 = aligned, 2 = staggered).');
end

Re_Dmax = V_max * D / nu;
end

function [Nu_D, C1, C2, m, Pr_s] = ExtConvTubeBank(fluid, Re_Dmax, T_i, T_o, T_s, S_T, S_L, N_L, tubeType)
% Case: External convection, tube bank
% Nu_D = ExtConvTubeBank(fluid, Re_Dmax, T_i, T_o, T_s, S_T, S_L, N_L, tubeType)
%   10 <= Re_D <= 2e6
%   0.7 <= Pr <= 500
%   Use T_bar (need a calc)
% Inputs:
%   Re_D, Pr, tubeType
%   tubeType = 1 (Aligned) or 2 (Staggered)
%   N_L = longitudinal pitch
%   fluid = 'air' or 'water'
%   T_i = temp coming into the tube bank (K)
%   T_o = temp coming out of the tube bank (K)
%   S_T = transverse pitch (m)
%   S_L = longitudinal pitch (m)
%   N_L = number of rows
%   V = velocity of free stream (m/s)
% Outputs:
%   Nu_D = avg Nusselt
%   C1 & C2 = from table, for double checking

% Get film temps
T_f = mean([T_i, T_o]);

Pr = getFluidProp(fluid, T_f, 'pr');
Pr_s = getFluidProp(fluid, T_s, 'pr');

% Find constants C1 and m from Table 7.5
ST_SL = S_T / S_L;

if tubeType == 1  % Aligned
    if Re_Dmax < 1e3
        C1 = 0.8; m = 0.4;
    elseif Re_Dmax < 2e5
        C1 = 0.27; m = 0.63;
    else
        C1 = 0.021; m = 0.84;
    end
elseif tubeType == 2  % Staggered
    if Re_Dmax < 1e3
        C1 = 0.9; m = 0.4;
    elseif Re_Dmax < 2e5
        if ST_SL < 2
            C1 = 0.35*(ST_SL)^(1/5); m = 0.6;
        else
            C1 = 0.4; m = 0.6;
        end
    else
        C1 = 0.022; m = 0.84;
    end
else
    warndlg('Invalid tube type (1 = aligned, 2 = staggered).');
end

% Interpolate correction factor C2 from Table 7.6 (use ch 18 interp1)
NL_table = [1 2 3 4 5 7 10 13 16];
if tubeType == 1
    C2_table = [0.70 0.80 0.86 0.90 0.92 0.95 0.97 0.98 0.99];
else
    C2_table = [0.64 0.76 0.84 0.89 0.92 0.95 0.97 0.98 0.99];
end

if N_L < 1
    warndlg('Error: longitudinal tube count < 1');
    C2 = NaN;
elseif N_L < 20
    C2 = interp1(NL_table, C2_table, N_L, 'linear');
    % we're using interp1 so we don't have to call polyint or another ftn
else
    C2 = 1;
end

% calc Nu_D
if Re_Dmax >= 10 && Re_Dmax <= 2e6 && Pr >= 0.68 && Pr <= 500 % really Pr >= .70, but some miss ok
    Nu_D = C1*C2*Re_Dmax^m*Pr^0.36*(Pr/Pr_s)^(1/4);
else
    warndlg(sprintf('Re_Dmax or Pr not ok.\nRe_L = %.2e\nPr = %.3f', Re_Dmax, Pr), 'Warning');
    Nu_D = C1*C2*Re_Dmax^m*Pr^0.36*(Pr/Pr_s)^(1/4);
end
end

%% Forced Internal
%Flow in a circular tube [Chapter 8]
% Case: 
function [Nu_D, convType] = IntCircularTubeV0(Re_D, Pr, Ts, Tm, L, D)

if Re_D <= 2300 && Pr >= 0.5 % laminar, thermal entry, uniform Ts
   convType = 'laminar'; 
       Nu_D = 3.66+(0.0668*G_zd)./(1+0.04*G_zd.^(2/3)); 
elseif Re_L >= 2300 && Pr >= 0.1 % laminar, combined entry, uniform Ts
   convType = 'laminar';
       %Nu_D=((3.66/tanh(2.264*G_zd.^(-1/3)+1.7*G_zd.^(-2/3))+0.0499*G_zd.*tanh(G_zd.^(-1)))./(tanh(2.432*Pr.^(1/6).*G_zd.^(-1/6));
       %has combined entry
elseif Re_D >= 10000 && 0.6 <= Pr && Pr <= 160 && L/D>10 && Ts > Tm % turbulent, fully developed (one of the most common)
   convType = 'turbulent';
       Nu_D=(0.023*Re_D.^(4/5).*Pr.^(0.3));

elseif Re_D >= 10000 && 0.6 <= Pr && Pr <= 160 && Ts < Tm % turbulent, fully developed (but Ts > Tm)
   convType = 'turbulent';
       Nu_D=(0.023*Re_D.^(4/5).*Pr.^(0.4));

elseif Re_D >= 10000 && 0.6 <= Pr && Pr <= 160 && Ts > Tm
   convType = 'turbulent';
       Nu_D=(0.023*Re_D.^(4/5).*Pr.^(0.3));
 
elseif 3600 <= Re_D && Re_D <= 905000 && 0.003 <= Pr && Pr <= 0.05 && 100 <= Re_D*Pr && Re_D*Pr <= 10000 % liquid metals, we dont' need
   convType = 'turbulent'; 
       Nu_D = 4.82+0.0185*(Re_D.*Pr).^(0.827);

elseif Re_D*Pr >= 100 % liquid metal
   convType = 'turbulent';
       Nu_D=5.0+0.025*(Re_D.*Pr).^(0.8);

else
    warning('Re_L or Pr outside acceptable range.');
    convType = 'RE_OUTSIDE';
end
end

function [Nu_D convType] = IntTube(fluid, Re_D, T_s, T_m, L, D, SurfCondit)
% Case: Internal convection, tube
% Nu_D = IntTube(fluid, Re_D, T_s, T_m, L, D, SurfCondit, x)
%   ___ <= Re_D <= ___
%   Assuming L/D > 10
% Inputs:
%   fluid, Re_D
%   T_s = surface temp (K)
%   T_m = mean temp in the tube (K)
%   LD = tube length L / diam D
%   SurfCondit = 0 -> none. 1 -> T_s uniform. 2 -> q" uniform
% Outputs:
%   Nu_D = avg Nusselt

LD = L/D; x = L; % if we wanted non x = L, we could prompt for it

if LD < 10
    warndlg('IntTube: L/D < 10. Results may be inaccurate.');
end

% Get fluid properties
Pr = getFluidProp(fluid, T_m, 'Pr');
mu = getFluidProp(fluid, T_m, 'mu');   
mu_s = getFluidProp(fluid, T_s, 'mu');

% Graetz number
G_zd = (D/x)*Re_D*Pr;

% Determine lam/turb
if Re_D < 2300
    lamTurb = 'laminar';
else
    lamTurb = 'turbulent';
end

% Handle lam/turb
switch lamTurb
    case 'laminar'
        G_zd_thresh = 100;  % we just picked this, kinda arbitrary
        
        if SurfCondit ~= 1 && SurfCondit ~= 2
            warndlg('SurfCondit missing or invalid. Assuming T_s uniform (SurfCondit=1).');
            SurfCondit = 1; % default to const temp
        end
        
        if SurfCondit == 2 && G_zd < G_zd_thresh
            Nu_D = 4.36;
            convType = 'Fully developed laminar, uniform q"';
            return;        
        else % T_s uniform
            if G_zd < G_zd_thresh
                % Fully developed laminar
                Nu_D = 3.66;
                convType = 'Fully developed laminar, uniform T_s';
                return;
            else % we just assume const surf temp
                if Pr >= 0.1
                    % Nu_D = [ 3.66 / tanh(2.264 G_zd^{-1/3} + 1.7 G_zd^{-2/3}) + 0.0499 G_zd tanh(G_zd^{-1}) ] ...
                    %        / tanh(2.432 Pr^{1/6} G_zd^{-1/6})
                    % guard against zero/negative arguments in power operations
                    A = 2.264 * (G_zd^(-1/3)) + 1.7 * (G_zd^(-2/3));
                    term1 = 3.66 / tanh(A);
                    term2 = 0.0499 * G_zd * tanh(G_zd^(-1));
                    denom = tanh(2.432 * Pr^(1/6) * G_zd^(-1/6));
                    Nu_D = (term1 + term2) / denom;
                    convType = 'Laminar, combined entry, uniform T_s';
                    return;
                else
                    Nu_D = 3.66 + (0.0668 * G_zd) / (1 + 0.04 * G_zd^(2/3));
                    convType = 'Laminar, thermal entry, uniform T_s';
                    return;
                end
            end
        end
        
    case 'turbulent'
        % Must be fully developed hydraulically & thermally
        if LD < 10
            warndlg('IntTube: L/D < 10. We assume L/D >= 10.');
        end
     
        if (Pr >= 0.6) && (Pr <= 160) && Re_D >= 10000 
            n = 0.3;
            if T_s > T_m
                n = 0.4;
            end
            Nu_D = 0.023 * Re_D^(4/5) * Pr^(n);
            convType = 'Turbulent, fully developed, eq 8.60';
            return;
        else
            if (Pr >= 0.7) && (Re_D >= 10000) && (LD >= 10)
                Nu_D = 0.027 * Re_D^(4/5) * Pr^(1/3) * (mu / mu_s)^0.14;
                convType = 'Turbulent, fully developed, eq 8.61';
                return;
            else
                warndlg('Pr or Re outside the acceptable range');
                Nu_D = 0.027 * Re_D^(4/5) * Pr^(1/3) * (mu / mu_s)^0.14;
                convType = 'Turbulent, fully developed, eq 8.61';
                return;
            end
        end
        
    otherwise
        error('IntTube: Unknown lamTurb detection error.');
end
end

function [T_s, T_mi, T_mo, T_m, T_f, SurfCondit, assumedVar, tol] = getTubeVals()
    clc;
    disp('Select a surface condition:');
    disp('  1. Uniform surface temperature T_s');
    disp('  2. Uniform surface heat flux q"');
    disp('  3. Other (enter T_inf)');
    SurfCondit = input('Surface condition (1–3): ');
    clc;

    if SurfCondit == 1 
        T_s = input('Enter the surface temperature, T_s (C): ')+273.15; % K
    elseif SurfCondit == 2 
        T_s = 273.15; % assume, it shouldn't matter
    elseif SurfCondit == 3
        T_inf = input('Enter the surrounding temperature, T_inf (C): ')+273.15; % K
        T_s = T_inf;
        % assume a thin walled tube
    end
    clc;

    T_mi = input('Enter the tube mean inlet temperature, T_mi (C): ')+273.15; % K
    T_mo = input('Enter the tube mean outlet temperature, T_mo (C): ')+273.15; % K
    clc;
    
    disp('Was any temperature assumed?');
    disp('  1. None');
    disp('  2. Outlet temperature (T_mo)');
    disp('  3. Inlet temperature (T_mi)');
    assumedVar = input('Assumption: ');
    if assumedVar ~= 1
        tol = input('Enter tolerance for iteration (percent): ');
    else
        tol = NaN;
    end
    clc;

    T_m = (T_mi+T_mo)/2;
    T_f = T_m; % helps us out in the end
end

function [Nu_D convType] = IntTubeNonCirc(fluid, Re_D, T_s, T_m, L, D_h, SurfCondit, geom, a,b)
% Case: Internal convection, tube, noncircular
% Nu_D = IntTube(fluid, Re_D, T_s, T_m, L, D, SurfCondit, x)
%   ___ <= Re_D <= ___
%   Assuming L/D > 10
% Inputs:
%   fluid, Re_D
%   T_s = surface temp (K)
%   T_m = mean temp in the tube (K)
%   LD = tube length L / diam D
%   SurfCondit = 0 -> none. 1 -> T_s uniform. 2 -> q" uniform
% Outputs:
%   Nu_D = avg Nusselt

LD = L / D_h;
if LD < 10
    warndlg('NonCircTube: L/D < 10. Fully-developed assumptions may be invalid.');
end

if Re_D < 2300 && LD >= 10
    if geom == 1 % rect
        b_over_a = b / a;
        b_a_vals = [1.0, 1.43, 2.0, 3.0, 4.0, 8.0, 1e6]; % b/a
        Nu_q_vals = [3.61, 3.73, 4.12, 4.79, 5.33, 6.49, 8.23]; % Nu for uniform q" 
        Nu_T_vals = [2.98, 3.08, 3.39, 3.96, 4.44, 5.60, 7.54]; % Nu for uniform T_s
        if SurfCondit == 1
            Nu_D = interp1(b_a_vals, Nu_T_vals, b_over_a, 'linear');
            convType = 'laminar, rect cross sect, fully developed, uniform T_s';
        else

            Nu_D = interp1(b_a_vals, Nu_q_vals, b_over_a, 'linear'); % chapter 18, linear interpolation ftn
            % this one's built in, so no need to include polyint, etc
            convType = 'laminar, rect cross sect, fully developed, uniform q"';
        end

        return;
        
    elseif geom == 2
        % Triangle (fully developed laminar)
        Nu_q_tri = 3.11;
        Nu_T_tri = 2.49;
        if SurfCondit == 2 % q" uniform
            Nu_D = Nu_q_tri;
            convType = 'laminar, triangle, fully developed, uniform q"';
        else % T_s uniform
            Nu_D = Nu_T_tri;
            convType = 'laminar, triangle, fully developed, uniform T_s';
        end
        return;
    end
else
    [Nu_D convType] = IntTube(fluid, Re_D, T_s, T_m, L, D_h, SurfCondit);
    return;
end
end

%% Free Convection

% Vertical plate - verified
function Nu_L = freeConvExtFlatPlateVert(Ra_L, Pr)
% Case: Free convection, external, vertical plate
%   Ra_L >= 10^4  ---- Works for all Ra_L, but another better for laminar
% Nu_L = freeConvExtFlatPlateVert(Ra_L, Pr)
% Inputs:
%   Ra_L, Pr
% Outputs:
%   Nu_L -- avg Nusselt
if Ra_L <= 10^9 && Ra_L >= 10^4
    Nu_L = 0.68 + 0.670*Ra_L.^(1/4)./ (1+(0.492/Pr).^(9/16)).^(4/9); % laminar, verified
elseif Ra_L >= 10^9
    Nu_L = (0.825+ 0.387*Ra_L.^(1/6) ./(1+ (0.492/Pr).^(9/16) ).^(8/27)) .^2; % verified
else
    Nu_L = NaN;
    warndlg(sprintf('Ra_L outside acceptable range.\nRe_L = %.2e\nPr = %.3f', Ra_L), 'Warning');
    % found online for troubleshooting
end
end

% Inclined plate
function Nu_L = freeConvExtPlateInc(Ra_L,Pr,theta)
% Case: Free convection, external, inclined plate, cold up or hot down
%   0 <= theta <= 60 (degrees)
%   Use g*cos(theta) instead of g for Ra_L calc (we do that outside getRa)
%       Just use the laminar vert plate
% Nu_L = freeConvExtPlateInc(Ra_L,Pr)
% Inputs:
%   Ra_L, Pr
% Outputs:
%   Nu_L

if theta >= 0 && theta <= 60
    Nu_L = freeConvExtFlatPlateVert(Ra_L, Pr); % use g*cos(theta) for Ra' % should default to vert lam
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
    warndlg('HotUpper: Invalid combination of Ra_L and Pr');
    Nu_L = 0.54*Ra_L.^(1/4);
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
    warndlg('HotLower: Invalid combination of Ra_L and Pr');
    Nu_L = 0.52*Ra_L.^(1/5);
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
    warndlg('ExtCyl: Invalid combination of Ra_D and Pr');
    Nu_D = (0.60 + 0.387*Ra_D.^(1/6) ./(1+(0.559/Pr).^(9/16)).^(8/27)).^2;
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
    warndlg('ExtSph: Invalid combination of Ra_D and Pr');
    Nu_D = 2 + 0.589*Ra_D.^(1/4) ./(1+(0.469/Pr).^(9/16)).^(4/9);
end
end

% Enclosures
function Nu_L = freeConvEncRectCavHorizHeatBelow(Ra_L,Pr)
% Case: Free convection, enclosed, rectangular cavity. Just a first approx
%   Ra_L <= 7e9, Ra_L >= 3e5
%   Nu_L = freeConvEncRectCav(Ra_L,Pr)
% Inputs:
%   Ra_L, Pr
% Outputs:
%   Nu_L

if Ra_L <= 7e9 && Ra_L >= 3e5
    Nu_L = 0.069*Ra_L^(1/3)*Pr^0.074;
else
    Nu_L = 0.069*Ra_L^(1/3)*Pr^0.074;
    warndlg('HeatBel: Invalid combination of Ra_L and Pr');
end
end

function Nu_L = freeConvEncRectVertCav(Ra_L,Pr,H,L)
% Case: Free convection, enclosed, rectangular vertical cavity.
%   Range of Ra_L
%   Nu_L = 
% Inputs:
%   Ra_L, Pr
%   H = height (along surfaces w/ dif temp)
%   L = length (dist btwn T1 and T2 faces)
% Outputs:
%   Nu_L

if Ra_L*Pr/(0.2+Pr) >= 10^3   && H/L >= 1 && H/L <= 2   && Pr >= 10^-3 && Pr <= 10^5
    Nu_L = 0.18*(Pr/(0.2+Pr)*Ra_L)^0.29;
elseif Ra_L >= 10^3 && Ra_L <= 10^10   && H/L >= 1 && H/L <= 10   && Pr <= 10^5
    Nu_L = 0.22*(Pr/(0.2+Pr)*Ra_L)^0.28*(H/L)^(-1/4);
elseif Ra_L >= 10^6 && Ra_L <= 10^9   && H/L >= 1 && H/L <= 40   && Pr >= 0.68 && Pr <= 20
    Nu_L = 0.046*Ra_L^(1/3);
    if Pr < 1
        warndlg('freeConvVertCav: Pr < 1');
    end
elseif Ra_L >= 10^4 && Ra_L <= 10^7   && H/L >= 10 && H/L <= 40   && Pr >= 0.68 && Pr <= 2e4
    Nu_L = 0.42*Ra_L^0.25*Pr^0.012*(H/L)^-0.3;
    if Pr < 1
        warndlg('freeConvVertCav: Pr < 1');
    end
else
    warndlg('VertCav: Invalid combination of Ra_L and Pr');
    Nu_L = NaN;
end
end

%% Interpolation Ftns
function val = getFluidProp(fluid, T_K, prop)
% getFluidProp: picks the right fluid interpolation ftn
% val = getFluidProp(T_K, fluid, prop)
%
% Inputs:
%   T_K = temperature (K)
%   fluid = 'air' or 'water'
%   prop = property name ('rho', 'cp', 'mu', 'k', 'Pr', etc.)
%
% Output:
%   val = interpolated property value

    switch lower(fluid)
        case 'air'
            val = getAirProp(T_K, prop);
        case 'water'
            val = getWaterProp(T_K, prop);
        otherwise
            error('Fluid not recognized. Use ''air'' or ''water''.');
    end
end

function val = getAirProp(T_K, prop)
% getAirProp: interpolates property tables for air
% val = getAirProp(T, prop)
% Inputs:
%   T_K = temperature (K)
%   prop = string property name to interpolate for
%           rho, cp, mu, k, Pr
% Output:
%   val = interpolated value

if T_K < 100 || T_K > 3000
    warndlg('Temperature outside acceptable range.');
    val = NaN;
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
      1.417 1.478 1.558 1.665 2.726]*1000;                              % Specific heat (J/kg·K)

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

    case 'pr'
        % Prandtl number (-)
        val = interp1(T_air, Pr_air, T_K, 'pchip');

    otherwise
        warndlg('Property not recognized. Use: rho, mu, cp, or k.');
end
end

function val = getWaterProp(T_K, prop)
% getAirProp: interpolates property tables for air
% val = getWaterProp(T, prop)
% Inputs:
%   T_K = temperature (K)
%   prop = string property name to interpolate for
%           rho, cp, mu, k, Pr
% Output:
%   val = interpolated value

if T_K < 260 || T_K > 650
    disp('Temperature outside acceptable range.');
    val = [];
    return;
end

% Temperature (K)
T_water = [ ...
 273.15 275 280 285 290 295 300 305 310 315 320 325 330 335 340 ...
 345 350 355 360 365 370 373.15 380 385 390 400 410 420 430 ...
 440 450 460 470 480 490 500 510 520 530 540 550 560 570 580 ...
 590 600 610 620 625 630 635 640 645 647.3];

% Specific volume (vf * 1e3) → will invert to rho
vf_water_1e3 = [ ...
 1.000 1.000 1.000 0.9994 0.9990 1.002 1.003 1.005 1.007 1.009 1.011 ...
 1.013 1.016 1.018 1.021 1.024 1.027 1.030 1.034 1.038 1.041 1.044 ...
 1.049 1.053 1.058 1.067 1.077 1.088 1.099 1.110 1.123 1.137 1.152 ...
 1.167 1.184 1.203 1.222 1.244 1.268 1.294 1.323 1.355 1.392 1.433 ...
 1.482 1.541 1.612 1.705 1.778 1.856 1.935 2.075 2.351 3.170];

% Specific heat (kJ/kg·K) — kept in kJ/kg·K here and multiplied to J/kg·K below
cp_water = [ ...
 4.217 4.211 4.198 4.189 4.184 4.181 4.179 4.178 4.178 4.179 4.180 ...
 4.182 4.184 4.186 4.188 4.191 4.195 4.199 4.203 4.209 4.214 4.217 ...
 4.220 4.232 4.239 4.256 4.273 4.302 4.331 4.360 4.400 4.440 4.480 ...
 4.530 4.590 4.660 4.740 4.840 4.950 5.080 5.240 5.430 5.680 6.000 ...
 6.410 7.000 7.850 9.350 10.600 12.600 16.400 26.000 90.000 Inf] * 1000; 

% Dynamic viscosity (μ_f * 1e6 Pa·s)
mu_water_1e6 = [ ...
 1750 1652 1422 1225 1080 959 855 769 695 631 577 528 489 453 420 ...
 389 365 343 324 306 289 279 260 248 237 217 200 185 173 162 152 ...
 143 136 129 124 118 113 108 104 101 97 94 91 88 84 81 77 72 ...
 70 67 64 59 54 45];

% Thermal conductivity (k_f * 1e3 W/m·K) 
k_water_1e3 = [ ...
 569 574 582 590 598 606 613 620 628 634 640 645 650 656 660 ...
 664 668 671 674 677 679 681 683 685 686 688 688 688 685 682 ...
 678 673 667 660 651 642 631 621 608 594 580 563 548 528 513 ...
 497 467 444 430 412 392 367 331 238];

% Prandtl number (-)
Pr_water = [ ...
 12.99 12.22 10.26 8.81 7.56 6.62 5.83 5.20 4.62 4.16 3.77 3.42 3.15 ...
 2.88 2.66 2.45 2.29 2.14 2.02 1.91 1.80 1.76 1.70 1.61 1.53 1.47 1.34 ...
 1.24 1.09 1.04 1.12 0.99 1.14 0.95 0.92 0.86 0.85 0.84 0.85 0.86 0.87 ...
 0.90 0.94 0.99 1.05 1.14 1.30 1.52 1.65 2.00 2.70 4.20 12.00 Inf];

% --- Convert to SI
vf_water   = vf_water_1e3 * 1e-3;        % m^3/kg
rho_water  = 1 ./ vf_water;              % kg/m^3
mu_water   = mu_water_1e6 * 1e-6;        % Pa·s
k_water    = k_water_1e3  * 1e-3;        % W/m·K
nu_water   = mu_water ./ rho_water;      % m^2/s

% Interpolate w/ interp1
switch lower(prop)
    case 'rho'
        val = interp1(T_water, rho_water, T_K, 'pchip');
    case 'mu'
        val = interp1(T_water, mu_water, T_K, 'pchip');
    case 'cp'
        val = interp1(T_water, cp_water, T_K, 'pchip');   % J/kg·K
    case 'k'
        val = interp1(T_water, k_water, T_K, 'pchip');
    case 'nu'
        val = interp1(T_water, nu_water, T_K, 'pchip');
    case 'pr'
        val = interp1(T_water, Pr_water, T_K, 'pchip');
    otherwise
        error('Property not recognized. Use: rho, mu, cp, k, nu, pr.');
end
end