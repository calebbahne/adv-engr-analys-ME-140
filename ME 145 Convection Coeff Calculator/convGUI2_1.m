function convGUI2_1
% Convection Coefficient Calculator - GUI wrapper (single-page style)
% - Hidden tabs (use full-page panels and a Back button).
% - Centered main inputs (Property, Film Temperature) with choices beneath.
% - Larger Results text. Responsive layout via relayout().

%% ===== Figure & five "pages" as panels =====
fig = uifigure('Name','Convection Coefficient Calculator','Position',[80 60 1040 700]);
fig.AutoResizeChildren = 'off';         % <-- allow SizeChangedFcn to run
% (prevents: "SizeChangedFcn callback will not execute while AutoResizeChildren is on")

% One bottom-right Back button (hidden on first page)
btnBack = uibutton(fig,'Text','← Back','FontSize',14, ...
    'Position',[fig.Position(3)-140, 16, 120, 36], ...
    'ButtonPushedFcn',@(~,~) goBack(), ...
    'Visible','off');

% Full-page panels (simulate tabs but hidden)
pMain   = uipanel(fig,'Title','General Input','Position',[10 60 1020 630]);      % Step 1
pForced = uipanel(fig,'Title','Forced External','Position',[10 60 1020 630], 'Visible','off');
pInt    = uipanel(fig,'Title','Forced Internal','Position',[10 60 1020 630], 'Visible','off');
pFree   = uipanel(fig,'Title','Free Convection','Position',[10 60 1020 630], 'Visible','off');
pRes    = uipanel(fig,'Title','Results','Position',[10 60 1020 630], 'Visible','off');

% Track current page (main|forced|internal|free|results)
currentPage = 'main';

%% ===== Main (General Input) =====
% Centered inputs (label above control)
lblProp = uilabel(pMain,'Text','Property:','FontSize',16, ...
    'HorizontalAlignment','center');
ddFluid = uidropdown(pMain,'Items',{'Air','Water'},'Value','Air','FontSize',16);

lblTf   = uilabel(pMain,'Text','Film Temperature T_f (°C):','FontSize',16, ...
    'HorizontalAlignment','center');
efTfC   = uieditfield(pMain,'numeric','Value',25,'FontSize',16);

% Big centered "choose next" title + buttons
lblNextTitle = uilabel(pMain, ...
    'Text','Choose where to go next', ...
    'FontSize',20, 'FontWeight','bold', ...
    'HorizontalAlignment','center');

btnGoForced   = uibutton(pMain,'Text','Click here for Forced External →', ...
    'FontSize',14, 'ButtonPushedFcn',@(~,~) showPage('forced'));
btnGoInternal = uibutton(pMain,'Text','Click here for Forced Internal →', ...
    'FontSize',14, 'ButtonPushedFcn',@(~,~) showPage('internal'));
btnGoFree     = uibutton(pMain,'Text','Click here for Free Convection →', ...
    'FontSize',14, 'ButtonPushedFcn',@(~,~) showPage('free'));

%% ===== Forced External =====
uilabel(pForced,'Text','Choose subcase:','Position',[20 580 120 24]);
ddExt = uidropdown(pForced,'Items',{'Flat plate','Cylinder','Sphere'}, ...
    'Value','Flat plate','Position',[140 580 170 24],'ValueChangedFcn',@(~,~)switchExt());

% Flat plate panel (ACTIVE by default; keep border/title)
pFlat = uipanel(pForced,'Title','Flat Plate in Parallel Flow','Position',[20 260 980 300], ...
    'BorderType','line');
uilabel(pFlat,'Text','Velocity V (m/s)','Position',[20 200 160 24]);
efV_flat = uieditfield(pFlat,'numeric','Value',1.0,'Position',[190 200 120 24]);
uilabel(pFlat,'Text','Plate Length L (m)','Position',[20 160 160 24]);
efL_flat = uieditfield(pFlat,'numeric','Value',1.0,'Position',[190 160 120 24]);
uibutton(pFlat,'Text','Compute heat coefficient (Flat Plate)','Position',[20 110 260 34],...
    'ButtonPushedFcn',@(~,~)onFlatPlate());

% Cylinder panel (INACTIVE initially: hide border/title to avoid ghost slivers)
pCyl = uipanel(pForced,'Title','', 'Position',[20 260 980 300], ...
    'Visible','off','BorderType','none');
uilabel(pCyl,'Text','Velocity V (m/s)','Position',[20 200 160 24]);
efV_cyl = uieditfield(pCyl,'numeric','Value',1.0,'Position',[190 200 120 24]);
uilabel(pCyl,'Text','Diameter D (m)','Position',[20 160 160 24]);
efD_cyl = uieditfield(pCyl,'numeric','Value',0.05,'Position',[190 160 120 24]);
uibutton(pCyl,'Text','Compute heat coefficient (Cylinder)','Position',[20 110 260 34],...
    'ButtonPushedFcn',@(~,~)onCylinder());

% Sphere panel (INACTIVE initially: hide border/title)
pSph = uipanel(pForced,'Title','', 'Position',[20 260 980 300], ...
    'Visible','off','BorderType','none');
uilabel(pSph,'Text','Velocity V (m/s)','Position',[20 200 160 24]);
efV_sph = uieditfield(pSph,'numeric','Value',1.0,'Position',[190 200 120 24]);
uilabel(pSph,'Text','Diameter D (m)','Position',[20 160 160 24]);
efD_sph = uieditfield(pSph,'numeric','Value',0.05,'Position',[190 160 120 24]);

bgMuS = uibuttongroup(pSph,'Title','Surface viscosity μ_s source','Position',[360 160 430 90]);
rbTs = uiradiobutton(bgMuS,'Text','Use surface temperature T_s (recommended)','Position',[10 50 360 20],'Value',true);
rbManual = uiradiobutton(bgMuS,'Text','Enter μ_s manually','Position',[10 25 300 20]);
bgMuS.SelectionChangedFcn = @(~,evt) toggleMuS(evt);

uilabel(pSph,'Text','Surface temperature T_s (°C)','Position',[20 120 180 24]);
efTsC = uieditfield(pSph,'numeric','Value',25,'Position',[210 120 120 24]);
uilabel(pSph,'Text','Surface viscosity μ_s (Pa·s)','Position',[20 90 180 24]);
efMuS_sph = uieditfield(pSph,'numeric','Value',1.5e-5,'Position',[210 90 120 24]); efMuS_sph.Enable='off';
uibutton(pSph,'Text','Compute heat coefficient (Sphere)','Position',[20 50 260 34],...
    'ButtonPushedFcn',@(~,~)onSphere());

    function toggleMuS(evt)
        if evt.NewValue == rbManual
            efMuS_sph.Enable = 'on';  efTsC.Enable = 'off';
        else
            efMuS_sph.Enable = 'off'; efTsC.Enable = 'on';
        end
    end

    function switchExt()
        % Anti-ghosting: only the active subpanel shows border/title
        switch ddExt.Value
            case 'Flat plate'
                pFlat.Visible = 'on';  pFlat.BorderType = 'line'; pFlat.Title = 'Flat Plate in Parallel Flow';
                pCyl.Visible  = 'off'; pCyl.BorderType  = 'none'; pCyl.Title  = '';
                pSph.Visible  = 'off'; pSph.BorderType  = 'none'; pSph.Title  = '';
                uistack(pFlat,'top');

            case 'Cylinder'
                pCyl.Visible  = 'on';  pCyl.BorderType  = 'line'; pCyl.Title  = 'Cylinder in Cross Flow';
                pFlat.Visible = 'off'; pFlat.BorderType = 'none'; pFlat.Title = '';
                pSph.Visible  = 'off'; pSph.BorderType  = 'none'; pSph.Title  = '';
                uistack(pCyl,'top');

            case 'Sphere'
                pSph.Visible  = 'on';  pSph.BorderType  = 'line'; pSph.Title  = 'Sphere';
                pFlat.Visible = 'off'; pFlat.BorderType = 'none'; pFlat.Title = '';
                pCyl.Visible  = 'off'; pCyl.BorderType  = 'none'; pCyl.Title  = '';
                uistack(pSph,'top');
        end
    end

%% ===== Forced Internal =====
uilabel(pInt,'Text','Circular Tube Flow:','Position',[20 580 140 24]);
ddInt = uidropdown(pInt,'Items',{'Laminar (fully developed, constant wall T)','Turbulent (Dittus–Boelter)'},...
    'Value','Laminar (fully developed, constant wall T)','Position',[170 580 330 24], ...
    'ValueChangedFcn',@(~,~)switchInt());

% Laminar panel (active default)
pIntLam = uipanel(pInt,'Title','Circular Tube — Laminar','Position',[20 260 980 300], 'BorderType','line');
uilabel(pIntLam,'Text','Velocity V (m/s)','Position',[20 200 160 24]);
efV_lam = uieditfield(pIntLam,'numeric','Value',0.2,'Position',[190 200 120 24]);
uilabel(pIntLam,'Text','Diameter D (m)','Position',[20 160 160 24]);
efD_lam = uieditfield(pIntLam,'numeric','Value',0.02,'Position',[190 160 120 24]);
uibutton(pIntLam,'Text','Compute heat coefficient (Laminar)','Position',[20 110 260 34],...
    'ButtonPushedFcn',@(~,~)onIntLam());

% Turbulent panel (inactive default: no border/title)
pIntTurb = uipanel(pInt,'Title','', 'Position',[20 260 980 300], 'Visible','off','BorderType','none');
uilabel(pIntTurb,'Text','Velocity V (m/s)','Position',[20 200 160 24]);
efV_turb = uieditfield(pIntTurb,'numeric','Value',2.0,'Position',[190 200 120 24]);
uilabel(pIntTurb,'Text','Diameter D (m)','Position',[20 160 160 24]);
efD_turb = uieditfield(pIntTurb,'numeric','Value',0.02,'Position',[190 160 120 24]);
uilabel(pIntTurb,'Text','Wall heat direction','Position',[360 200 160 24]);
ddDir = uidropdown(pIntTurb,'Items',{'Heating fluid (n=0.4)','Cooling fluid (n=0.3)'},...
    'Value','Heating fluid (n=0.4)','Position',[520 200 200 24]);
uibutton(pIntTurb,'Text','Compute heat coefficient (Turbulent)','Position',[20 110 260 34],...
    'ButtonPushedFcn',@(~,~)onIntTurb());

    function switchInt()
        switch ddInt.Value
            case 'Laminar (fully developed, constant wall T)'
                pIntLam.Visible  = 'on';  pIntLam.BorderType  = 'line'; pIntLam.Title  = 'Circular Tube — Laminar';
                pIntTurb.Visible = 'off'; pIntTurb.BorderType = 'none'; pIntTurb.Title = '';
                uistack(pIntLam,'top');

            case 'Turbulent (Dittus–Boelter)'
                pIntTurb.Visible = 'on';  pIntTurb.BorderType = 'line'; pIntTurb.Title = 'Circular Tube — Turbulent (Dittus–Boelter)';
                pIntLam.Visible  = 'off'; pIntLam.BorderType  = 'none'; pIntLam.Title  = '';
                uistack(pIntTurb,'top');
        end
    end

%% ===== Free Convection (Immersed) =====
uilabel(pFree,'Text','Geometry:','Position',[20 580 100 24]);
ddFree = uidropdown(pFree,'Items',{'Flat plate','Horizontal cylinder','Sphere'},...
    'Value','Flat plate','Position',[100 580 160 24],'ValueChangedFcn',@(~,~)switchFree());

% ---- Sub-panels (same size so we can reuse positions) ----
pFreePlate = uipanel(pFree,'Title','Flat Plate','Position',[20 210 980 360]);
pFreeCyl   = uipanel(pFree,'Title','Horizontal Cylinder','Position',[20 210 980 360],'Visible','off');
pFreeSph   = uipanel(pFree,'Title','Sphere','Position',[20 210 980 360],'Visible','off');

% Common positions for temp inputs INSIDE sub-panels
tsLblPos   = [20 300 200 24];
tsEditPos  = [200 300 120 24];
tinfLblPos = [425 300 260 24];
tinfEditPos= [640 300 90  24];

% ---- Shared temperature controls (initially parented to Flat Plate panel) ----
lblTs   = uilabel(pFreePlate,'Text','Surface temperature T_s (°C)','Position',tsLblPos);
efTs_free   = uieditfield(pFreePlate,'numeric','Value',50,'Position',tsEditPos);

lblTinf = uilabel(pFreePlate,'Text','Free Stream Temperature T_inf (°C)','Position',tinfLblPos);
efTinf_free = uieditfield(pFreePlate,'numeric','Value',25,'Position',tinfEditPos);

% ---- Flat plate sub-options (same as your original) ----
uilabel(pFreePlate,'Text','Orientation','Position',[20 260 120 24]);
ddOrient = uidropdown(pFreePlate,'Items',{'Vertical','Horizontal: hot surface up','Horizontal: hot surface down','Inclined (cold up/hot down)'},...
    'Value','Vertical','Position',[120 260 260 24],'ValueChangedFcn',@(~,~)toggleIncl());
uilabel(pFreePlate,'Text','Length L (m)','Position',[20 220 140 24]);
efL_free = uieditfield(pFreePlate,'numeric','Value',0.5,'Position',[120 220 120 24]);
uilabel(pFreePlate,'Text','Inclination θ (deg)','Position',[360 220 140 24]);
efTheta = uieditfield(pFreePlate,'numeric','Value',30,'Position',[510 220 100 24]); 
efTheta.Enable='off';
uibutton(pFreePlate,'Text','Compute heat coefficient (Flat Plate)','Position',[20 170 260 34],...
    'ButtonPushedFcn',@(~,~)onFreePlate());

% ---- Cylinder sub-options ----
uilabel(pFreeCyl,'Text','Diameter D (m)','Position',[20 260 140 24]);
efD_free_cyl = uieditfield(pFreeCyl,'numeric','Value',0.05,'Position',[120 260 120 24]);
uibutton(pFreeCyl,'Text','Compute heat coefficient (Cylinder)','Position',[20 210 260 34],...
    'ButtonPushedFcn',@(~,~)onFreeCyl());

% ---- Sphere sub-options ----
uilabel(pFreeSph,'Text','Diameter D (m)','Position',[20 260 140 24]);
efD_free_sph = uieditfield(pFreeSph,'numeric','Value',0.05,'Position',[120 260 120 24]);
uibutton(pFreeSph,'Text','Compute heat coefficient (Sphere)','Position',[20 210 260 34],...
    'ButtonPushedFcn',@(~,~)onFreeSph());

% ---- Switcher + helpers ----
    function switchFree()
        % toggle panels
        plate = strcmp(ddFree.Value,'Flat plate');
        cyl   = strcmp(ddFree.Value,'Horizontal cylinder');
        sph   = strcmp(ddFree.Value,'Sphere');
        pFreePlate.Visible = plate;  pFreeCyl.Visible = cyl;  pFreeSph.Visible = sph;

        % reparent the shared temperature inputs to whichever panel is visible
        if plate
            lblTs.Parent = pFreePlate;     efTs_free.Parent = pFreePlate;
            lblTinf.Parent = pFreePlate;   efTinf_free.Parent = pFreePlate;
        elseif cyl
            lblTs.Parent = pFreeCyl;       efTs_free.Parent = pFreeCyl;
            lblTinf.Parent = pFreeCyl;     efTinf_free.Parent = pFreeCyl;
        else % sphere
            lblTs.Parent = pFreeSph;       efTs_free.Parent = pFreeSph;
            lblTinf.Parent = pFreeSph;     efTinf_free.Parent = pFreeSph;
        end
        % keep positions consistent inside the chosen panel
        lblTs.Position    = tsLblPos;    efTs_free.Position    = tsEditPos;
        lblTinf.Position  = tinfLblPos;  efTinf_free.Position  = tinfEditPos;
    end

    function toggleIncl()
        if strcmp(ddOrient.Value,'Inclined (cold up/hot down)'); efTheta.Enable='on';
        else; efTheta.Enable='off'; end
    end

%% ===== Results (bigger text, fills panel with margins) =====
ta = uitextarea(pRes,'Position',[24 24 972 582],'Editable','off','FontName','Consolas','FontSize',16);
ta.Value = "Results will appear here.";

%% ===== Navigation helpers =====
% Show a page and manage Back visibility
function showPage(which)
    currentPage = which;
    pMain.Visible   = strcmp(which,'main');
    pForced.Visible = strcmp(which,'forced');
    pInt.Visible    = strcmp(which,'internal');
    pFree.Visible   = strcmp(which,'free');
    pRes.Visible    = strcmp(which,'results');

    if strcmp(which,'main')
        btnBack.Visible = 'off';
    else
        btnBack.Visible = 'on';
    end
    relayout();  % keep things centered
end

% Back: all pages -> main; results -> main
function goBack()
    switch currentPage
        case {'forced','internal','free','results'}
            showPage('main');
        otherwise
            showPage('main');
    end
end

% Responsive relayout (center the main inputs and choices; keep Back anchored)
fig.SizeChangedFcn = @(~,~) relayout();
function relayout()
    pad = 20;
    fw  = fig.InnerPosition(3);

    % Bottom-right Back button
    btnBack.Position = [fw-140-pad, 16, 120, 36];

    % Use panel width for centering inside pMain
    pmPos = pMain.InnerPosition; pw = pmPos(3);
    centerX = pmPos(1) + pw/2;

    % Sizes
    h = 28; gap = 10;                      % control height, label→control gap
    lblW = 300; fldW = 260;                % widths for label and control

    % --- Block 1: Property (label on its own line, dropdown below) ---
    yLabel1 = 540;                         % label Y
    yCtrl1  = yLabel1 - (h + gap);         % control Y just below
    lblProp.Position = [centerX - lblW/2, yLabel1, lblW, h];
    ddFluid.Position = [centerX - fldW/2, yCtrl1,  fldW, h];

    % --- Block 2: Film Temperature (label on its own line, edit below) ---
    yLabel2 = yCtrl1 - 2*(h + gap);        % add some vertical spacing
    yCtrl2  = yLabel2 - (h + gap);
    lblTf.Position   = [centerX - lblW/2, yLabel2, lblW, h];
    efTfC.Position   = [centerX - fldW/2, yCtrl2,  fldW, h];

    % --- Title + buttons block beneath ---
    titleY = yCtrl2 - 50;                  % space below temp input
    lblWbig = 520; lblH = 34;
    lblNextTitle.Position = [centerX - lblWbig/2, titleY, lblWbig, lblH];

    btnW = 420; btnH = 48; btnGap = 16;
    btnGoForced.Position   = [centerX - btnW/2, titleY - 60, btnW, btnH];
    btnGoInternal.Position = [centerX - btnW/2, titleY - 60 - (btnH+btnGap), btnW, btnH];
    btnGoFree.Position     = [centerX - btnW/2, titleY - 60 - 2*(btnH+btnGap), btnW, btnH];
end


% Start at main
showPage('main');

%% ===== Shared helpers =====
    function [fluidStr, T_K] = getFluidAndTempK()
        fluidStr = lower(string(ddFluid.Value));   % 'air' or 'water'
        T_K = efTfC.Value + 273.15;
    end
    function [nu, Pr, k, mu, rho, cp_kJ, alpha] = getProps(fluidStr, T_K)
        switch fluidStr
            case "air"
                nu = getAirProp(T_K,'nu');  Pr = getAirProp(T_K,'pr');  k  = getAirProp(T_K,'k');
                mu = getAirProp(T_K,'mu');  rho= getAirProp(T_K,'rho'); cp_kJ = getAirProp(T_K,'cp');
                alpha = getAirProp(T_K,'alpha');
            case "water"
                nu = getWaterProp(T_K,'nu'); Pr = getWaterProp(T_K,'pr'); k  = getWaterProp(T_K,'k');
                mu = getWaterProp(T_K,'mu'); rho= getWaterProp(T_K,'rho'); cp_kJ = getWaterProp(T_K,'cp');
                alpha = getWaterProp(T_K,'alpha');
            otherwise
                error('Unknown fluid');
        end
    end
    function out = propBlock(fluidStr,T_K)
        switch fluidStr
            case "air"
                rho   = getAirProp(T_K,'rho');  cp_kJ = getAirProp(T_K,'cp');
                mu    = getAirProp(T_K,'mu');   alpha = getAirProp(T_K,'alpha');
                nu    = getAirProp(T_K,'nu');
            case "water"
                rho   = getWaterProp(T_K,'rho'); cp_kJ = getWaterProp(T_K,'cp');
                mu    = getWaterProp(T_K,'mu');  alpha = getWaterProp(T_K,'alpha');
                nu    = getWaterProp(T_K,'nu');
        end
        fmtProps = [ ...
            '\n--- %s Properties at %.2f K ---\n' ...
            'rho:    %.6g kg/m³\n' ...
            'cp:     %.6g J/kg·K\n' ...
            'mu:     %.6g Pa·s\n' ...
            'nu:     %.6g m²/s\n' ...
            'alpha:  %.6g m²/s\n' ];
        out = sprintf(fmtProps, char(upper(string(fluidStr))), T_K, rho, cp_kJ*1000, mu, nu, alpha);
    end
    function lines = composeResults(fluidStr, caseTitle, T_K, dimName, dimVal, V, Re, Pr, Nu, k, h)
        fmtHeader = [ ...
            '=== Convection Coefficient Results ===\n' ...
            'Fluid:                             %s\n' ...
            'Case:                              %s\n' ...
            'Film Temperature (T_f):            %.2f K\n' ...
            'Velocity (V):                      %.6g m/s\n' ...
            '%s:                                %.6g m\n' ...
            'Re:                                %.6g\n' ...
            'Pr:                                %.6g\n' ...
            'Nu:                                %.6g\n' ...
            'k:                                 %.6g W/m·K\n' ...
            'h:                                 %.6g W/m²·K\n' ];
        header = sprintf(fmtHeader, char(upper(string(fluidStr))), char(caseTitle), ...
            T_K, V, char(dimName), dimVal, Re, Pr, Nu, k, h);
        lines = string(header) + newline + string(propBlock(fluidStr,T_K));
    end
    function lines = composeResultsFree(fluidStr, caseTitle, T_K, dimName, dimVal, Ra, Pr, Nu, k, h)
        fmtHeader = [ ...
            '=== Convection Coefficient Results ===\n' ...
            'Fluid:                             %s\n' ...
            'Case:                              %s\n' ...
            'Film Temperature (T_f):            %.2f K\n' ...
            '%s:                                %.6g m\n' ...
            'Ra:                                %.6g\n' ...
            'Pr:                                %.6g\n' ...
            'Nu:                                %.6g\n' ...
            'k:                                 %.6g W/m·K\n' ...
            'h:                                 %.6g W/m²·K\n' ];
        header = sprintf(fmtHeader, char(upper(string(fluidStr))), char(caseTitle), ...
            T_K, char(dimName), dimVal, Ra, Pr, Nu, k, h);
        lines = string(header) + newline + string(propBlock(fluidStr,T_K));
    end

%% ===== Forced External callbacks =====
    function onFlatPlate()
        [fluidStr, T_K] = getFluidAndTempK();
        [nu, Pr, k, ~, ~, ~, ~] = getProps(fluidStr, T_K);
        V = efV_flat.Value; L = efL_flat.Value; Re_L = V*L/nu;
        [Nu_L, convType] = ExtConvFlatPlate(Re_L, Pr);
        if isempty(Nu_L), uialert(fig,'Inputs outside correlation range for flat plate.','Warning'); return; end
        h = Nu_L * k / L;
        out = composeResults(fluidStr, ['Forced External → Flat Plate (' char(convType) ')'], ...
            T_K, 'L', L, V, Re_L, Pr, Nu_L, k, h);
        ta.Value = out; showPage('results');
    end

    function onCylinder()
        [fluidStr, T_K] = getFluidAndTempK();
        [nu, Pr, k, ~, ~, ~, ~] = getProps(fluidStr, T_K);
        V = efV_cyl.Value; D = efD_cyl.Value; Re_D = V*D/nu;
        Nu_D = ExtConvCyl(Re_D, Pr);
        if isempty(Nu_D), uialert(fig,'Inputs outside correlation range for cylinder.','Warning'); return; end
        h = Nu_D * k / D;
        out = composeResults(fluidStr, 'Forced External → Cylinder', ...
            T_K, 'D', D, V, Re_D, Pr, Nu_D, k, h);
        ta.Value = out; showPage('results');
    end

    function onSphere()
        [fluidStr, T_K] = getFluidAndTempK();
        [nu, Pr, k, mu, ~, ~, ~] = getProps(fluidStr, T_K);
        V = efV_sph.Value; D = efD_sph.Value; Re_D = V*D/nu;

        % μ_s source
        if rbTs.Value
            Ts_K = efTsC.Value + 273.15;
            if fluidStr == "air"
                mu_s = getAirProp(Ts_K,'mu');
            elseif fluidStr == "water"
                mu_s = getWaterProp(Ts_K,'mu');
            else
                error('Unknown fluid: %s', char(fluidStr));
            end
        else
            mu_s = efMuS_sph.Value;
        end

        ratio = mu/mu_s;
        if ~(Re_D>=3.5 && Re_D<=7.6e4) || ~(Pr>=0.71 && Pr<=380) || ~(ratio>=1 && ratio<=3.2)
            uialert(fig, sprintf( ...
                ['Inputs outside sphere correlation range:\n' ...
                 'Re_D=%.3g (need 3.5–7.6e4)\nPr=%.3g (need 0.71–380)\nμ/μ_s=%.3g (need 1.0–3.2)\n' ...
                 'Tip: set T_s ≈ T_f so μ≈μ_s; adjust V or D.'], Re_D, Pr, ratio),'Out of range');
        end
        Nu_D = ExtConvSphere(Re_D, Pr, mu, mu_s);
        if isempty(Nu_D), return; end
        h = Nu_D*k/D;
        out = composeResults(fluidStr, 'Forced External → Sphere', T_K, 'D', D, V, Re_D, Pr, Nu_D, k, h);
        out = [out; sprintf('\nμ/μ_s used: %.3g   (μ=%.3g Pa·s, μ_s=%.3g Pa·s)\n', ratio, mu, mu_s)];
        ta.Value = out; showPage('results');
    end

%% ===== Forced Internal callbacks =====
    function onIntLam()
        % Fully developed laminar, constant wall T: Nu = 3.66 (valid Re<2300)
        [fluidStr, T_K] = getFluidAndTempK();
        [nu, Pr, k, ~, ~, ~, ~] = getProps(fluidStr, T_K);
        V = efV_lam.Value; D = efD_lam.Value; Re = V*D/nu;
        if Re >= 2300
            uialert(fig,sprintf('Warning: Re=%.3g ≥ 2300, not laminar. Result may be invalid.',Re),'Warning');
        end
        Nu = 3.66;
        h  = Nu * k / D;
        out = composeResults(fluidStr, 'Forced Internal → Circular tube (Laminar FD, Tw const)', ...
            T_K, 'D', D, V, Re, Pr, Nu, k, h);
        ta.Value = out; showPage('results');
    end

    function onIntTurb()
        % Dittus–Boelter: Nu = 0.023 Re^0.8 Pr^n ; n=0.4 heating, n=0.3 cooling
        [fluidStr, T_K] = getFluidAndTempK();
        [nu, Pr, k, ~, ~, ~, ~] = getProps(fluidStr, T_K);
        V = efV_turb.Value; D = efD_turb.Value; Re = V*D/nu;
        n = 0.4; if contains(ddDir.Value,'Cooling'), n=0.3; end
        if Re < 1e4
            uialert(fig,sprintf('Warning: Re=%.3g < 10,000, Dittus–Boelter outside standard range.',Re),'Warning');
        end
        Nu = 0.023 * Re^0.8 * Pr^n;
        h  = Nu * k / D;
        out = composeResults(fluidStr, ['Forced Internal → Circular tube (Dittus–Boelter, n=' num2str(n) ')'], ...
            T_K, 'D', D, V, Re, Pr, Nu, k, h);
        ta.Value = out; showPage('results');
    end

%% ===== Free Convection (Immersed) callbacks =====
    function [Ra, Pr, k, TfK] = freeRa(fluidStr, TsC, TinfC, Lc)
        TsK = TsC + 273.15; TinfK = TinfC + 273.15; TfK = (TsK+TinfK)/2;
        [nu, Pr, k, ~, ~, ~, ~] = getProps(fluidStr, TfK);
        g = 9.81; beta = 1/TfK; dT = TsK - TinfK;
        Ra = abs(g*beta*dT*Lc^3 / (nu^2) * Pr);
    end

    function onFreePlate()
        [fluidStr, ~] = getFluidAndTempK();
        TsC = efTs_free.Value; TinfC = efTinf_free.Value; L = efL_free.Value;
        [Ra, Pr, k, TfK] = freeRa(fluidStr, TsC, TinfC, L);
        switch ddOrient.Value
            case 'Vertical'
                Nu = freeConvExtFlatPlateVert(Ra,Pr);
                caseTitle = 'Free Convection → Flat Plate (Vertical)';
            case 'Horizontal: hot surface up'
                Nu = freeConvExtPlateHorizHotUpper(Ra,Pr);
                caseTitle = 'Free Convection → Flat Plate (Horizontal, hot up)';
            case 'Horizontal: hot surface down'
                Nu = freeConvExtPlateHorizHotLower(Ra,Pr);
                caseTitle = 'Free Convection → Flat Plate (Horizontal, hot down)';
            case 'Inclined (cold up/hot down)'
                theta = efTheta.Value;
                Nu = freeConvExtPlateInc(Ra,Pr,theta); % user note: Ra should use g*cosθ ideally
                caseTitle = sprintf('Free Convection → Flat Plate (Inclined, θ=%.1f°)',theta);
        end
        if isempty(Nu), uialert(fig,'Inputs outside correlation range for flat plate.','Warning'); return; end
        h = Nu * k / L;
        out = composeResultsFree(fluidStr, caseTitle, TfK, 'L', L, Ra, Pr, Nu, k, h);
        ta.Value = out; showPage('results');
    end

    function onFreeCyl()
        [fluidStr, ~] = getFluidAndTempK();
        TsC = efTs_free.Value; TinfC = efTinf_free.Value; D = efD_free_cyl.Value;
        [Ra, Pr, k, TfK] = freeRa(fluidStr, TsC, TinfC, D);
        Nu = freeConvExtHorizCyl(Ra,Pr);
        if isempty(Nu), uialert(fig,'Inputs outside correlation range for horizontal cylinder.','Warning'); return; end
        h = Nu * k / D;
        out = composeResultsFree(fluidStr, 'Free Convection → Horizontal Cylinder', TfK, 'D', D, Ra, Pr, Nu, k, h);
        ta.Value = out; showPage('results');
    end

    function onFreeSph()
        [fluidStr, ~] = getFluidAndTempK();
        TsC = efTs_free.Value; TinfC = efTinf_free.Value; D = efD_free_sph.Value;
        [Ra, Pr, k, TfK] = freeRa(fluidStr, TsC, TinfC, D);
        Nu = freeConvExtSphere(Ra,Pr);
        if isempty(Nu), uialert(fig,'Inputs outside correlation range for sphere.','Warning'); return; end
        h = Nu * k / D;
        out = composeResultsFree(fluidStr, 'Free Convection → Sphere', TfK, 'D', D, Ra, Pr, Nu, k, h);
        ta.Value = out; showPage('results');
    end

end % ===== end convGUI2 (GUI layer) =====

%% ===== Your correlation functions (unchanged names) =====
function [Nu_L, convType] = ExtConvFlatPlate(Re_L, Pr)
if Re_L <= 5e5 && Pr >= 0.6
    convType = 'laminar'; Nu_L = 0.664*Re_L.^0.5 .* Pr.^(1/3);
elseif Re_L >= 5e5 && Re_L <= 1e8 && Pr >= 0.6 && Pr <= 60
    convType = 'mixed';   Nu_L = (0.037*Re_L.^0.8 - 871) .* Pr.^(1/3);
else
    warning('Re_L or Pr outside acceptable range.'); Nu_L = []; convType = 'RE_OUTSIDE';
end
end

function Nu_D = ExtConvCyl(Re_D, Pr)
if Re_D*Pr >= 0.2
    Nu_D = 0.3 + (0.62*Re_D.^0.5 .* Pr.^(1/3) .* (1+(0.4./Pr).^(2/3)).^(-1/4)) ...
              .* (1 + (Re_D/282000).^(5/8)).^(4/5);
else, disp('Re_D*Pr outside acceptable range.'); Nu_D = [];
end
end

function Nu_D = ExtConvSphere(Re_D, Pr, mu, mu_s)
if (Re_D >= 3.5 && Re_D <= 7.6e4) && (Pr >= 0.71 && Pr <= 380) && (mu./mu_s >= 1 && mu./mu_s <= 3.2)
    Nu_D = 2 + (0.4*Re_D.^0.5 + 0.06*Re_D.^(2/3)) .* Pr.^0.4 .* (mu./mu_s).^(1/4);
else, disp('Re_D, Pr, or μ/μ_s outside acceptable range.'); Nu_D = [];
end
end

function Nu_L = freeConvExtFlatPlateVert(Ra_L, Pr)
if Ra_L <= 1e9 && Ra_L >= 1e4
    Nu_L = 0.68 + 0.670*Ra_L.^(1/4)./(1+(0.492/Pr).^(9/16)).^(4/9);
elseif Ra_L >= 1e9
    Nu_L = (0.825+0.387*Ra_L.^(1/6)./(1+(0.492/Pr).^(9/16)).^(8/27)).^2;
else, Nu_L = []; disp('Error in freeConvExtFlatPlateVert, Ra_L outside accepted range');
end
end

function Nu_L = freeConvExtPlateInc(Ra_L,Pr,theta)
if theta>=0 && theta<=60
    Nu_L = freeConvExtFlatPlateVert(Ra_L, Pr);
else, disp('theta outside range for inclined plate.'); Nu_L=[];
end
end

function Nu_L = freeConvExtPlateHorizHotUpper(Ra_L,Pr)
if Ra_L >= 1e4 && Ra_L <= 1e7 && Pr >= 0.7
    Nu_L = 0.54*Ra_L.^(1/4);
elseif Ra_L >= 1e7 && Ra_L <= 1e11
    Nu_L = 0.15*Ra_L.^(1/3);
else, disp('Invalid combination Ra_L/Pr'); Nu_L=[];
end
end

function Nu_L = freeConvExtPlateHorizHotLower(Ra_L,Pr)
if Ra_L >= 1e4 && Ra_L <= 1e9 && Pr >= 0.7
    Nu_L = 0.52*Ra_L.^(1/5);
else, disp('Invalid combination Ra_L/Pr'); Nu_L=[];
end
end

function Nu_D = freeConvExtHorizCyl(Ra_D,Pr)
if Ra_D <= 1e12
    Nu_D = (0.60 + 0.387*Ra_D.^(1/6) ./(1+(0.559/Pr).^(9/16)).^(8/27)).^2;
else, disp('Invalid Ra_D'); Nu_D=[];
end
end

function Nu_D = freeConvExtSphere(Ra_D,Pr)
if Ra_D <= 1e12 && Pr >= 0.7
    Nu_D = 2 + 0.589*Ra_D.^(1/4) ./(1+(0.469/Pr).^(9/16)).^(4/9);
else, disp('Invalid Ra_D/Pr'); Nu_D=[];
end
end


%% ===== Property tables (as you provided) =====
function val = getAirProp(T_K, prop)
if T_K < 100 || T_K > 3000, disp('Temperature outside acceptable range.'); val = []; return; end
T_air = [100 150 200 250 300 350 400 450 500 550 600 650 700 750 800 ...
     850 900 950 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900 ...
     2000 2100 2200 2300 2400 2500 3000];
rho_air = [3.5562 2.3364 1.7458 1.3947 1.1614 0.9950 0.8711 0.7740 0.6964 ...
       0.6329 0.5804 0.5356 0.4975 0.4643 0.4354 0.4097 0.3868 0.3666 ...
       0.3482 0.3166 0.2902 0.2679 0.2488 0.2322 0.2177 0.2049 0.1935 ...
       0.1833 0.1741 0.1658 0.1582 0.1513 0.1448 0.1389 0.1135];
cp_air = [1.032 1.012 1.007 1.006 1.007 1.009 1.014 1.021 1.030 1.040 ...
      1.051 1.063 1.075 1.087 1.099 1.110 1.121 1.131 1.141 1.159 ...
      1.175 1.189 1.207 1.230 1.248 1.267 1.286 1.307 1.337 1.372 ...
      1.417 1.478 1.558 1.665 2.726];
mu_air = [71.1 103.4 132.5 159.6 184.6 208.2 230.1 250.7 270.1 288.4 305.8 ...
      322.5 338.8 354.6 369.8 384.3 398.1 411.3 424.4 449.0 473.0 496.0 ...
      530.0 557.0 584.0 611.0 637.0 663.0 689.0 715.0 740.0 766.0 792.0 ...
      818.0 955.0] * 1e-7;
nu_air = [2.00 4.426 7.590 11.44 15.89 20.92 26.41 32.39 38.79 45.57 52.69 ...
      60.21 68.10 76.37 84.93 93.80 102.9 112.2 121.9 141.8 162.9 185.1 ...
      213.0 240.0 268.0 298.0 329.0 362.0 396.0 431.0 468.0 506.0 547.0 ...
      589.0 841.0] * 1e-6;
k_air = [9.34 13.8 18.1 22.3 26.3 30.0 33.8 37.3 40.7 43.9 46.9 49.7 52.4 ...
     54.9 57.3 59.6 62.0 64.3 66.7 71.5 76.3 82.0 91.0 100.0 106.0 ...
     113.0 120.0 128.0 137.0 147.0 160.0 175.0 196.0 222.0 486.0] * 1e-3;
alpha_air = [2.54 5.84 10.3 15.9 22.5 29.9 38.3 47.2 56.7 66.7 76.9 87.3 ...
          98.0 109.0 120.0 130.0 143.0 155.0 168.0 195.0 224.0 257.0 ...
          303.0 350.0 390.0 435.0 482.0 534.0 589.0 646.0 714.0 783.0 ...
          869.0 960.0 1570.0] * 1e-6;
Pr_air = [0.786 0.758 0.737 0.720 0.707 0.700 0.690 0.686 0.684 0.683 ...
      0.685 0.690 0.695 0.702 0.709 0.716 0.720 0.723 0.726 0.728 ...
      0.728 0.719 0.703 0.685 0.688 0.685 0.683 0.677 0.672 0.667 ...
      0.655 0.647 0.630 0.613 0.536];
switch lower(prop)
    case 'rho',   val = interp1(T_air, rho_air,  T_K, 'pchip');
    case 'mu',    val = interp1(T_air, mu_air,   T_K, 'pchip');
    case 'cp',    val = interp1(T_air, cp_air,   T_K, 'pchip');
    case 'k',     val = interp1(T_air, k_air,    T_K, 'pchip');
    case 'nu',    val = interp1(T_air, nu_air,   T_K, 'pchip');
    case 'alpha', val = interp1(T_air, alpha_air,T_K, 'pchip');
    case 'pr',    val = interp1(T_air, Pr_air,   T_K, 'pchip');
    otherwise, error('Property not recognized. Use: rho, mu, cp, k, nu, alpha, pr.');
end
end

function val = getWaterProp(T_K, prop)
if T_K < 273.15 || T_K > 647.3, disp('Temperature outside acceptable range.'); val = []; return; end
T_water = [273.15 275 280 285 290 295 300 305 310 315 320 325 330 335 340 ...
 345 350 355 360 365 370 373.15 380 400 420 440 460 480 500 520 ...
 540 560 580 600 620 640 647.3];
vf_water_1e3 = [1.000 1.000 1.000 1.000 1.001 1.001 1.003 1.003 1.007 1.012 1.011 ...
 1.013 1.016 1.017 1.013 1.024 1.032 1.030 1.030 1.038 1.041 1.043 ...
 1.049 1.058 1.088 1.110 1.137 1.152 1.184 1.244 1.294 1.392 1.433 ...
 1.541 1.770 2.197 2.110];
cp_water = [4.217 4.215 4.198 4.191 4.184 4.179 4.179 4.178 4.178 4.178 4.180 ...
 4.182 4.186 4.192 4.188 4.191 4.195 4.203 4.209 4.216 4.224 4.217 ...
 4.226 4.256 4.232 4.236 4.198 4.138 4.059 3.984 3.725 2.589 2.263 ...
 2.007 0.941 0.563 0.000];
mu_water_1e6 = [1750 1654 1422 1225 1080 959 855 765 695 639 577 532 489 453 420 ...
 389 364 343 326 306 289 279 264 237 215 162 138 123 112 104 ...
 95 91 90 81 72 64 45];
k_water_1e3 = [561 574 582 590 598 606 613 619 626 632 640 645 650 656 660 ...
 664 671 679 686 677 688 682 688 688 688 682 673 667 651 628 ...
 594 548 515 497 444 330 238];
Pr_water = [12.8 12.3 11.0 10.2 9.41 8.73 7.57 6.83 6.24 5.73 5.10 4.65 4.26 ...
 3.96 3.71 3.41 3.20 3.02 2.81 2.61 2.44 1.76 1.61 1.35 1.16 ...
 1.04 0.99 0.92 0.87 0.84 0.86 0.94 0.99 1.14 1.52 2.70 Inf];

vf_water   = vf_water_1e3 * 1e-3;
rho_water  = 1 ./ vf_water;
mu_water   = mu_water_1e6 * 1e-6;
k_water    = k_water_1e3  * 1e-3;
nu_water   = mu_water ./ rho_water;
alpha_water= k_water ./ (rho_water .* (cp_water*1000));

switch lower(prop)
    case 'rho',   val = interp1(T_water, rho_water,   T_K, 'pchip');
    case 'mu',    val = interp1(T_water, mu_water,    T_K, 'pchip');
    case 'cp',    val = interp1(T_water, cp_water,    T_K, 'pchip');
    case 'k',     val = interp1(T_water, k_water,     T_K, 'pchip');
    case 'nu',    val = interp1(T_water, nu_water,    T_K, 'pchip');
    case 'alpha', val = interp1(T_water, alpha_water, T_K, 'pchip');
    case 'pr',    val = interp1(T_water, Pr_water,    T_K, 'pchip');
    otherwise,    error('Property not recognized. Use: rho, mu, cp, k, nu, alpha, pr.');
end
end