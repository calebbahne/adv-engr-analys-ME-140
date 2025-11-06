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
