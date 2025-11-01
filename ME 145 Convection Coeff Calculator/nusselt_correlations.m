% nusselt_correlations.m
%
% Collection of average Nusselt number correlations. Each function implements
% one correlation from the provided table and includes short comments about
% the case and prerequisites. All functions accept scalars or arrays (elementwise ops).
%
% NOTE: The user requested only the function definitions. No example calls or
% additional logic are included here.

%% 1) Laminar, fully-developed, uniform heat flux q'' (constant Nu)

function Nu = Nu_D_laminar_fullydev_qflux()
% Case: Laminar, fully-developed internal flow, uniform heat flux (q'')
% Returns: constant Nusselt number for this case.
% Prereqs: fully developed thermally, laminar (Re_D small), constant q''.
Nu = 4.36;
end

%% 2) Laminar, fully-developed, uniform surface temperature Ts (constant Nu)

function Nu = Nu_D_laminar_fullydev_Ts()
% Case: Laminar, fully-developed internal flow, uniform surface temperature Ts
% Returns: constant Nusselt number for this case.
% Prereqs: fully developed thermally, laminar, constant Ts.
Nu = 3.66;
end

%% 3) Laminar, thermal entry (or combined entry with Pr >= ~5), uniform Ts

function Nu = Nu_D_laminar_thermal_entry(GzD)
% Case: Laminar, thermal entry (or combined entry with Pr >= ~5), uniform Ts
% Input:
%   GzD = Graetz-like parameter Gz_D = (D/x) * Re_D * Pr  (or appropriate definition)
% Prereqs: GzD > 0. Use when entrance effects important.
Nu = 3.66 + 0.0668 .* GzD ./ (1 + 0.04 .* (GzD).^(2/3));
end


%% 4) Laminar, combined entry, Pr >= ~0.1, uniform Ts

function Nu = Nu_D_laminar_combined_entry(GzD, Pr)
% Case: Laminar, combined thermal/velocity entry (Pr >= ~0.1), uniform Ts
% Inputs:
%   GzD = Graetz-like parameter (D/x * Re_D * Pr or as used in your reference)
%   Pr   = Prandtl number
% Prereqs: laminar flow, combined entry conditions; valid for Pr >= ~0.1.
num = 3.66 .* tanh(2.264 .* (GzD).^(-1/3) + 1.7 .* (GzD).^- (2/3)) + 0.0499 .* GzD .* tanh((GzD).^-1);
den = tanh(2.432 .* Pr.^(1/6) .* (GzD).^- (1/6));
Nu = num ./ den;
end

%% 7) Turbulent fully-developed (Dittus-Boelter)

function Nu = Nu_D_turbulent_dittus_boelter(ReD, Pr, n)
% Case: Turbulent, fully-developed internal flow, Dittus-Boelter correlation.
% Inputs:
%   ReD = Reynolds number
%   Pr  = Prandtl number
%   n   = exponent (typically 0.4 for heating (Ts > Tbulk), 0.3 for cooling)
% Prereqs: turbulent, Re_D > ~10,000 typical; valid Pr range often 0.6 < Pr < 160
Nu = 0.023 .* (ReD).^0.8 .* (Pr).^n;
end


%% 8) Turbulent fully-developed (Sieder–Tate style viscosity correction)

function Nu = Nu_D_turbulent_sieder_tate(ReD, Pr, mu_bulk, mu_surface)
% Case: Turbulent, fully-developed with viscosity correction (Sieder-Tate style)
% Inputs:
%   ReD        = Reynolds number
%   Pr         = Prandtl number
%   mu_bulk    = dynamic viscosity at bulk/film temperature (or mu)
%   mu_surface = dynamic viscosity at surface temperature (mu_s)
% Prereqs: turbulent, ReD >= 10,000 typically; valid where viscosity variation matters.
Nu = 0.027 .* (ReD).^0.8 .* Pr.^(1/3) .* (mu_bulk ./ mu_surface).^0.14;
end


%% 9) Gnielinski correlation (turbulent, fully-developed)

function Nu = Nu_D_gnielinski(ReD, Pr, f)
% Case: Turbulent, fully-developed internal flow — Gnielinski correlation.
% Inputs:
%   ReD = Reynolds number
%   Pr  = Prandtl number
%   f   = Darcy friction factor (use friction_factor_smooth or colebrook_friction)
% Prereqs: 0.5 <= Pr <= 2000 typically; 3000 <= ReD <= 5e6 and (L/D) >= ~10 often suggested
Nu = (f./8) .* (ReD - 1000) .* Pr ./ (1 + 12.7 .* sqrt(f./8) .* (Pr.^(2/3) - 1));
end


%% 10) Liquid-metal / very low Pr case (correlation variant)

function Nu = Nu_D_liquid_metal_1(ReD, Pr)
% Case: Liquid metals, turbulent, fully developed, uniform q'' (from table)
% Inputs:
%   ReD, Pr
% Prereqs: low-Pr fluids (liquid metals); check table limits (e.g., Pr ~ 3.6e-3 ...)
Nu = 4.82 + 0.0185 .* (ReD .* Pr).^0.827;
end


%% 11) Liquid-metal / alternate (from table)

function Nu = Nu_D_liquid_metal_2(ReD, Pr)
% Case: Liquid metals, turbulent, fully developed, uniform Ts
% Inputs:
%   ReD, Pr
% Prereqs: low-Pr fluids; check the reference ranges before using.
Nu = 5.0 + 0.025 .* (ReD .* Pr).^0.8;
end
