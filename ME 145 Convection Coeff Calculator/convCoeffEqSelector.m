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
disp('  1. External Flow');
disp('  2. Internal Flow');
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
            case 2
                clc;
                disp('External conv, cylinder');
            case 3
                clc;
                disp('External conv, sphere');
            case 4
                clc;
                disp('External conv, Bank of tubes');
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
                            case 2
                                clc;
                                disp('Select a horizontal plate:');
                                disp('  1. Hot surface up / cold surface down');
                                disp('  2. Cold surface up / hot surface down');
                                hotcold = input('Orientation (1–2): ');
                                if hotcold == 1
                                    clc;
                                    disp('Free conv, external, flat plate, hot surface up / cold surface down');
                                elseif hotcold == 2
                                    clc;
                                    disp('Free conv, external, flat plate, cold surface up / hot surface down');
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
                    case 3
                        clc;
                        disp('Free conv, external, sphere');
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
disp(' ');
disp('Output the following:');
disp('  Selected case type');
disp('  Re, Pr, Ra, geometry stuff, Nu');
disp('  h')

disp('----------------------------------');
disp('Convection coefficient calculation complete');
disp('   We love you Troy, so please give us the A :)');
disp('----------------------------------');
