%Reynolds Number Equation Solver
clear
clc
%Coefficients Used
%V = Velocity (m/s)
%D = Diameter (m)
%nu = Kinematic Viscosity
%mu = Dynamic Viscosity (Pa*s)
%rho = Fluid Density
%Write Information Here!

Options = input(['Enter 1 if you have Density and Dynamic Viscosity or ' ...
    '2 if you have Kinematic viscosity:\n']);

if Options == 1
    rho = input('Enter your fluid density:\n');
    mu = input('Enter your Dynamic Viscosity:\n');
    D = input('Enter your Diameter:\n');
    V = input('Enter your Velocity:\n');
    v = mu / rho; %Just to have it already be here
    Re = (rho * V * D) / mu;

elseif Options == 2
    v = input('Enter your Kinematic Viscosity:\n');
    D = input('Enter your Diameter:\n');
    V = input('Enter your Velocity:\n');
    Re = (V * D) / v;
else
    error('Please either choose 1 or 2')

end

if Re < 2300
    flow = 'Laminar';
elseif Re <= 4000
    flow = 'Transitional';
else
    flow = 'Turbulent';
end

Pr = input('Enter your Prandtl Number:\n');
k = input('Please input k (Thermal Conducitiy):\n');
%---------------------------------------------
%For Laminar
if flow == "Laminar"
Option2 = input(['Enter 1 if you have uniform qs" Enter 2 if you have' ...
    ' uniform T_s \n or Enter 3 if both options dont apply:\n']);
if Option2 == 1

    Nu_D = 4.36;
    h = (Nu_D*k)/D;
    fprintf('Your convection coefficient is:\n')
    disp(h)

elseif Option2 == 2
    Nu_D = 3.66;
    h = (Nu_D*k)/D;
    fprintf('Your convection coefficient is:\n')
    disp(h)

elseif Option2 == 3
    x = input(['Please input the distance from the tube inlet in meters' ...
        ':\n']);
    Gzd = (D/x)*Re*Pr;
    %This is the second condition which consists of Pr
    if Pr >= 5
    Nu_D = 3.66 + (0.0668 * Gzd) / (1 + 0.04 * Gzd^(2/3));
    h = (Nu_D*k)/D;
    fprintf('Your convection coefficient is:\n')
    disp(h)
    elseif Pr > 0.1
    Nu_D = ((3.66/(tanh(2.264*Gzd^(-1/3) + 1.7 * Gzd^(-2/3))) + 0.0499*Gzd*tanh(Gzd^(-1)))/ (tanh(2.432*Pr^(1/6) * Gzd^(-1/6))));
    h = (Nu_D*k)/D;
    fprintf('Your convection coefficient is:\n')
    disp(h)
    else
        error('You made a mistake somewhere, please recheck your values')
    end
else
    error('Please input a number that is available')
end
%----------------------------------------------
%For Turbulent
elseif flow == "Turbulent"
f = (0.790*log(Re) - 1.64)^-2;
Option3 = input(['Please choose 1 if you are dealing with Liquid Metals\n' ...
    'or Choose 2 if you are NOT dealing with Liquid Metals:\n']);
if Option3 == 1 %These options are going for 8.6 and over
    if Pr*Re > 100
        Nu_D = 5 + 0.025*(Re*Pr)^0.8;
        h = (Nu_D*k)/D;
        fprintf('Your convection coefficient is:\n')
        disp(h)
    else
        Nu_D = 4.82 + 0.0185*(Re*Pr)^0.827;
        h = (Nu_D*k)/D;
        fprintf('Your convection coefficient is:\n')
        disp(h)
    end
elseif Option3 == 2
    L = input('Please enter the length:\n');
    LD = L/D;
    if Re > 3000 && Re < 5e6 && Pr > 0.5 && Pr < 2000 && LD > 10
       Nu_D = ((f/8)*(Re - 1000)*Pr) / (1 + 12.7*(f/8)^0.5*(Pr^(2/3) - 1));
       h = (Nu_D*k)/D;
       fprintf('Your convection coefficient is:\n')
       disp(h)

    elseif Re >= 1e4 && Pr > 0.7 && Pr < 16700 && LD > 10
        mus = input('Please input your Dynamic Viscosity at the wall:\n');
        Nu_D = 0.027*Re^(4/5) * Pr^(1/3) * (mu/mus)^0.14;
        h = (Nu_D*k)/D;
        fprintf('Your convection coefficient is:\n')
        disp(h)

    elseif Re >= 1e4 && Pr > 0.6 && Pr < 160 && LD > 10
       T_s = input('Please input your Surface Temperture (C):\n');
       T_f = input('Please input your fluid Temperture (C):\n');
       if T_s > T_f
       n = 0.4;
       else
       n = 0.3;
       end
       Nu_D = 0.023 * Re^0.8 * Pr^n;
       h = (Nu_D*k)/D;
       fprintf('Your convection coefficient is:\n')
       disp(h)

    else
        error('Please recheck your inputs')
    end
else
    error('Choose either 1 or 2')
end

else 
    %There is nothing for Transitional Flow
    error('Use another method to solve for convection')
end