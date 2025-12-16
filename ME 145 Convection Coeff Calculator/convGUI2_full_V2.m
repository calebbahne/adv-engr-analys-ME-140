function convGUI2_full
%% Convection Coefficient Calculator
% NOTE: The GUI functionality of this code was generated using AI. 
%       The GUI was added on top of the convCoeffEqSelector code we wrote.

% Purpose: Calculate the avg convection coefficient h for a given condition
%   Determines the appropriate case, then does calcs from there.
%   -  14 cases are currently supported. 
%   -  All cases have been verified for at least one example problem
%       - 24 total test trials recorded, with mean error 1.9%
%       - 63% of trials had less than 3% error, 96% had less than 5% error
%       - Maximum error was 8% (free convection vertical plate)
%   - Most problems can be input to the script within 45 sec (max 90 sec)
%       - Most problems take ~ 12+ minutes to compute by hand
%       - Tube bank problems can be entered within 90 sec 
%           - Tube problems take ~ 30+ min to hand calc (per iteration)
%% ===== Figure & five "pages" as panels =====
fig = uifigure('Name','Convection Coefficient Calculator','Position',[80 60 1120 740]);
fig.AutoResizeChildren = 'off';
btnBack = uibutton(fig,'Text','← Back','FontSize',14, ...
    'Position',[fig.Position(3)-160, 16, 140, 36], ...
    'ButtonPushedFcn',@(~,~) goBack(), 'Visible','off');

pMain   = uipanel(fig,'Title','General Input','Position',[10 60 1100 670]);
pForced = uipanel(fig,'Title','Forced External','Position',[10 60 1100 670], 'Visible','off');
pInt    = uipanel(fig,'Title','Forced Internal','Position',[10 60 1100 670], 'Visible','off');
pFree   = uipanel(fig,'Title','Free Convection (Immersed)','Position',[10 60 1100 670], 'Visible','off');
pEnc    = uipanel(fig,'Title','Free Convection (Enclosed)','Position',[10 60 1100 670], 'Visible','off');
pRes    = uipanel(fig,'Title','Results','Position',[10 60 1100 670], 'Visible','off');

lastNonResultsPage = 'main';   % remembers where user was before navigating to Results
currentPage        = 'main';   % track current page for global back

%% ===== Main (Property and T_f with label-above-control layout) =====
lblProp = uilabel(pMain,'Text','Property (Fluid):','FontSize',16,'HorizontalAlignment','center');
ddFluid = uidropdown(pMain,'Items',{'Air','Water'},'Value','Air','FontSize',16);

lblTf   = uilabel(pMain,'Text','Film Temperature T_f (°C):','FontSize',16,'HorizontalAlignment','center');
efTfC   = uieditfield(pMain,'numeric','Value',25,'FontSize',16);

lblNextTitle = uilabel(pMain,'Text','Choose where to go next','FontSize',20,'FontWeight','bold','HorizontalAlignment','center');

btnGoForced   = uibutton(pMain,'Text','Forced External →','FontSize',14,'ButtonPushedFcn',@(~,~) showPage('forced'));
btnGoInternal = uibutton(pMain,'Text','Forced Internal →','FontSize',14,'ButtonPushedFcn',@(~,~) showPage('internal'));
btnGoFree     = uibutton(pMain,'Text','Free Convection (Immersed) →','FontSize',14,'ButtonPushedFcn',@(~,~) showPage('free'));
btnGoEnclosed = uibutton(pMain,'Text','Free Convection (Enclosed) →','FontSize',14,'ButtonPushedFcn',@(~,~) showPage('enclosed'));

%% ===== Forced External =====
uilabel(pForced,'Text','Subcase:','Position',[20 620 80 24]);
ddExt = uidropdown(pForced,'Items',{'Flat plate','Cylinder','Sphere','Bank of tubes'}, ...
    'Value','Flat plate','Position',[100 620 200 24],'ValueChangedFcn',@(~,~)switchExt());

% ---- Flat plate panel
pFlat = uipanel(pForced,'Title','Flat Plate in Parallel Flow','Position',[20 240 1060 360],'BorderType','line');
cbReFlat = uicheckbox(pFlat,'Text','Enter Re_L directly','Position',[20 300 200 22],'Value',false);
uilabel(pFlat,'Text','Velocity V (m/s)','Position',[20 260 160 24]);
efV_flat = uieditfield(pFlat,'numeric','Value',1.0,'Position',[190 260 120 24]);
uilabel(pFlat,'Text','Plate Length L (m)','Position',[20 220 160 24]);
efL_flat = uieditfield(pFlat,'numeric','Value',1.0,'Position',[190 220 120 24]);
uilabel(pFlat,'Text','Re_L (direct)','Position',[360 260 120 24]);
efReL_flat = uieditfield(pFlat,'numeric','Value',1e5,'Position',[480 260 120 24],'Enable','off');
uibutton(pFlat,'Text','Compute h (Flat Plate)','Position',[20 170 220 34],'ButtonPushedFcn',@(~,~)onFlatPlate());
cbReFlat.ValueChangedFcn = @(~,~) setReDirect(cbReFlat,efReL_flat,[efV_flat,efL_flat]);

% ---- Cylinder panel
pCyl = uipanel(pForced,'Title','Cylinder in Cross Flow','Position',[20 240 1060 360], ...
    'BorderType','none','Visible','off');
cbReCyl = uicheckbox(pCyl,'Text','Enter Re_D directly','Position',[20 300 200 22],'Value',false);
uilabel(pCyl,'Text','Velocity V (m/s)','Position',[20 260 160 24]);
efV_cyl = uieditfield(pCyl,'numeric','Value',1.0,'Position',[190 260 120 24]);
uilabel(pCyl,'Text','Diameter D (m)','Position',[20 220 160 24]);
efD_cyl = uieditfield(pCyl,'numeric','Value',0.05,'Position',[190 220 120 24]);
uilabel(pCyl,'Text','Re_D (direct)','Position',[360 260 120 24]);
efReD_cyl = uieditfield(pCyl,'numeric','Value',1e4,'Position',[480 260 120 24],'Enable','off');
uibutton(pCyl,'Text','Compute h (Cylinder)','Position',[20 170 220 34],'ButtonPushedFcn',@(~,~)onCylinder());
cbReCyl.ValueChangedFcn = @(~,~) setReDirect(cbReCyl,efReD_cyl,[efV_cyl,efD_cyl]);

% ---- Sphere panel
pSph = uipanel(pForced,'Title','Sphere','Position',[20 240 1060 360], ...
    'BorderType','none','Visible','off');
cbReSph = uicheckbox(pSph,'Text','Enter Re_D directly','Position',[20 300 200 22],'Value',false);
uilabel(pSph,'Text','Velocity V (m/s)','Position',[20 260 160 24]);
efV_sph = uieditfield(pSph,'numeric','Value',1.0,'Position',[190 260 120 24]);
uilabel(pSph,'Text','Diameter D (m)','Position',[20 220 160 24]);
efD_sph = uieditfield(pSph,'numeric','Value',0.05,'Position',[190 220 120 24]);
uilabel(pSph,'Text','Re_D (direct)','Position',[360 260 120 24]);
efReD_sph = uieditfield(pSph,'numeric','Value',1e4,'Position',[480 260 120 24],'Enable','off');
bgMuS = uibuttongroup(pSph,'Title','Surface viscosity μ_s source','Position',[640 230 380 110]);
rbTs = uiradiobutton(bgMuS,'Text','Use surface temperature T_s (recommended)','Position',[10 70 360 20],'Value',true);
rbManual = uiradiobutton(bgMuS,'Text','Enter μ_s manually','Position',[10 45 300 20]);
uilabel(pSph,'Text','T_s (°C)','Position',[640 190 80 24]);
efTsC = uieditfield(pSph,'numeric','Value',25,'Position',[720 190 100 24]);
uilabel(pSph,'Text','μ_s (Pa·s)','Position',[840 190 80 24]);
efMuS_sph = uieditfield(pSph,'numeric','Value',1.5e-5,'Position',[920 190 100 24]); efMuS_sph.Enable='off';
uibutton(pSph,'Text','Compute h (Sphere)','Position',[20 170 220 34],'ButtonPushedFcn',@(~,~)onSphere());
cbReSph.ValueChangedFcn = @(~,~) setReDirect(cbReSph,efReD_sph,[efV_sph,efD_sph]);
bgMuS.SelectionChangedFcn = @(~,evt) toggleMuS(evt,efTsC,efMuS_sph);

% ---- Tube bank panel
pBank = uipanel(pForced,'Title','Bank of Tubes','Position',[20 70 1060 530],'BorderType','none','Visible','off');
% Temperatures
uilabel(pBank,'Text','T_i (°C)','Position',[20 470 80 24]);   efTi = uieditfield(pBank,'numeric','Value',20,'Position',[90 470 80 24]);
uilabel(pBank,'Text','T_o (°C)','Position',[180 470 80 24]);  efTo = uieditfield(pBank,'numeric','Value',40,'Position',[250 470 80 24]);
uilabel(pBank,'Text','T_s (°C)','Position',[340 470 80 24]);  efTs = uieditfield(pBank,'numeric','Value',60,'Position',[410 470 80 24]);
% Assumed variable + tolerance
uilabel(pBank,'Text','Assume:','Position',[520 470 60 24]);
ddAssume = uidropdown(pBank,'Items',{'None','Outlet T_o','Inlet T_i'},'Value','None','Position',[580 470 140 24]);
uilabel(pBank,'Text','Tolerance (%)','Position',[740 470 100 24]); efTol = uieditfield(pBank,'numeric','Value',1,'Position',[840 470 80 24]);
% Arrangement
uilabel(pBank,'Text','Arrangement','Position',[20 430 90 24]);
ddArrange = uidropdown(pBank,'Items',{'Aligned','Staggered'},'Value','Aligned','Position',[110 430 120 24]);
% Geometry
uilabel(pBank,'Text','D (m)','Position',[250 430 60 24]);    efD_b = uieditfield(pBank,'numeric','Value',0.02,'Position',[300 430 80 24]);
uilabel(pBank,'Text','S_T (m)','Position',[390 430 60 24]);  efST  = uieditfield(pBank,'numeric','Value',0.04,'Position',[450 430 80 24]);
uilabel(pBank,'Text','S_L (m)','Position',[540 430 60 24]);  efSL  = uieditfield(pBank,'numeric','Value',0.04,'Position',[600 430 80 24]);
uilabel(pBank,'Text','Rows N_L','Position',[690 430 80 24]); efNL  = uieditfield(pBank,'numeric','Value',5,'Position',[770 430 80 24]);
uilabel(pBank,'Text','Cols N_T','Position',[860 430 80 24]); efNT  = uieditfield(pBank,'numeric','Value',8,'Position',[940 430 80 24]);
% Flow
uilabel(pBank,'Text','Free-stream V (m/s)','Position',[20 390 140 24]); efVb = uieditfield(pBank,'numeric','Value',5,'Position',[160 390 90 24]);

uibutton(pBank,'Text','Compute Tube Bank','Position',[20 340 200 34],'ButtonPushedFcn',@(~,~)onTubeBank());

tbOut = uitextarea(pBank,'Position',[20 20 1020 300],'Editable','off','FontName','Consolas');

%% ===== Forced Internal =====
uilabel(pInt,'Text','Subcase:','Position',[20 620 80 24]);
ddInt = uidropdown(pInt,'Items',{'Circular — Laminar (fully developed, T_w const)', ...
                                 'Circular — Turbulent (Dittus–Boelter)', ...
                                 'Noncircular ducts'}, ...
                    'Value','Circular — Laminar (fully developed, T_w const)', ...
                    'Position',[100 620 420 24], 'ValueChangedFcn',@(~,~)switchInt());

% ---- Circular — Laminar
pIntLam = uipanel(pInt,'Title','Circular Tube — Laminar FD, T_w const','Position',[20 240 1060 360],'BorderType','line');
cbReIntL = uicheckbox(pIntLam,'Text','Enter Re directly','Position',[20 300 160 22],'Value',false);
uilabel(pIntLam,'Text','V (m/s)','Position',[20 260 80 24]);  efV_lam = uieditfield(pIntLam,'numeric','Value',0.2,'Position',[90 260 100 24]);
uilabel(pIntLam,'Text','D (m)','Position',[210 260 80 24]);    efD_lam = uieditfield(pIntLam,'numeric','Value',0.02,'Position',[260 260 100 24]);
uilabel(pIntLam,'Text','Re (direct)','Position',[380 260 100 24]); efRe_intL = uieditfield(pIntLam,'numeric','Value',1000,'Position',[480 260 120 24],'Enable','off');
uibutton(pIntLam,'Text','Compute h (Laminar)','Position',[20 210 220 34],'ButtonPushedFcn',@(~,~)onIntLam());
cbReIntL.ValueChangedFcn = @(~,~) setReDirect(cbReIntL,efRe_intL,[efV_lam,efD_lam]);

% ---- Circular — Turbulent
pIntTurb = uipanel(pInt,'Title','Circular Tube — Turbulent (Dittus–Boelter)','Position',[20 240 1060 360],'BorderType','none','Visible','off');
cbReIntT = uicheckbox(pIntTurb,'Text','Enter Re directly','Position',[20 300 160 22],'Value',false);
uilabel(pIntTurb,'Text','V (m/s)','Position',[20 260 80 24]);  efV_turb = uieditfield(pIntTurb,'numeric','Value',2.0,'Position',[90 260 100 24]);
uilabel(pIntTurb,'Text','D (m)','Position',[210 260 80 24]);    efD_turb = uieditfield(pIntTurb,'numeric','Value',0.02,'Position',[260 260 100 24]);
uilabel(pIntTurb,'Text','Re (direct)','Position',[380 260 100 24]); efRe_intT = uieditfield(pIntTurb,'numeric','Value',2e4,'Position',[480 260 120 24],'Enable','off');
uilabel(pIntTurb,'Text','Wall heat direction','Position',[620 260 140 24]);
ddDir = uidropdown(pIntTurb,'Items',{'Heating fluid (n=0.4)','Cooling fluid (n=0.3)'},'Value','Heating fluid (n=0.4)','Position',[760 260 220 24]);
uibutton(pIntTurb,'Text','Compute h (Turbulent)','Position',[20 210 220 34],'ButtonPushedFcn',@(~,~)onIntTurb());
cbReIntT.ValueChangedFcn = @(~,~) setReDirect(cbReIntT,efRe_intT,[efV_turb,efD_turb]);

% ---- Noncircular ducts (GUI, hydraulic-diameter method)
pIntNon = uipanel(pInt,'Title','Noncircular ducts (hydraulic diameter method)', ...
                  'Position',[20 240 1060 360],'BorderType','none','Visible','off');

% Reynolds entry mode
cbReNon = uicheckbox(pIntNon,'Text','Enter Re_D directly', ...
    'Position',[20 300 180 22],'Value',false);

% Temperatures and length
uilabel(pIntNon,'Text','Mean fluid temperature T_m (°C)', ...
    'Position',[20 260 220 24]);
efTm_non = uieditfield(pIntNon,'numeric','Value',25, ...
    'Position',[250 260 100 24]);

uilabel(pIntNon,'Text','Surface temperature T_s (°C)', ...
    'Position',[380 260 200 24]);
efTs_non = uieditfield(pIntNon,'numeric','Value',35, ...
    'Position',[580 260 100 24]);

uilabel(pIntNon,'Text','Duct length L (m)', ...
    'Position',[700 260 120 24]);
efL_non = uieditfield(pIntNon,'numeric','Value',1.0, ...
    'Position',[820 260 100 24]);

% Geometry (rectangular a x b)
uilabel(pIntNon,'Text','Cross-section type', ...
    'Position',[20 220 150 24]);
ddGeom_non = uidropdown(pIntNon,'Items',{'Rectangular duct a×b'}, ...
    'Value','Rectangular duct a×b','Position',[170 220 220 24]);

uilabel(pIntNon,'Text','Width a (m)', ...
    'Position',[20 180 120 24]);
efA_non = uieditfield(pIntNon,'numeric','Value',0.05, ...
    'Position',[140 180 100 24]);

uilabel(pIntNon,'Text','Height b (m)', ...
    'Position',[260 180 120 24]);
efB_non = uieditfield(pIntNon,'numeric','Value',0.02, ...
    'Position',[380 180 100 24]);

% Flow
uilabel(pIntNon,'Text','Bulk velocity V (m/s)', ...
    'Position',[20 140 150 24]);
efV_non = uieditfield(pIntNon,'numeric','Value',1.0, ...
    'Position',[170 140 100 24]);

uilabel(pIntNon,'Text','Re_D (direct)', ...
    'Position',[300 140 120 24]);
efRe_non = uieditfield(pIntNon,'numeric','Value',4000, ...
    'Position',[420 140 120 24],'Enable','off');

cbReNon.ValueChangedFcn = @(~,~) setReDirect(cbReNon,efRe_non,[efV_non]);

% Compute button + note
uibutton(pIntNon,'Text','Compute h (Noncircular duct)', ...
    'Position',[20 90 260 34], ...
    'ButtonPushedFcn',@(~,~)onIntNon());

uilabel(pIntNon,'Text','Uses hydraulic diameter D_h = 4A_c/P and pipe correlations (laminar: Nu=3.66; turbulent: Dittus–Boelter with D_h).', ...
    'Position',[20 40 1020 24],'WordWrap','on');

%% ===== Free Convection (Immersed) =====
uilabel(pFree,'Text','Geometry:','Position',[20 620 100 24]);
ddFree = uidropdown(pFree,'Items',{'Flat plate','Horizontal cylinder','Sphere'},'Value','Flat plate','Position',[110 620 200 24],'ValueChangedFcn',@(~,~)switchFree());

% subpanels same size
pFreePlate = uipanel(pFree,'Title','Flat Plate','Position',[20 210 1060 400]);
pFreeCyl   = uipanel(pFree,'Title','Horizontal Cylinder','Position',[20 210 1060 400],'Visible','off');
pFreeSph   = uipanel(pFree,'Title','Sphere','Position',[20 210 1060 400],'Visible','off');

% Common temp inputs (shared controls re-parented)
tsLblPos   = [20 340 200 24];  tsEditPos  = [200 340 120 24];
tinfLblPos = [425 340 260 24]; tinfEditPos= [690 340 100 24];

lblTs   = uilabel(pFreePlate,'Text','Surface temperature T_s (°C)','Position',tsLblPos);
efTs_free   = uieditfield(pFreePlate,'numeric','Value',50,'Position',tsEditPos);
lblTinf = uilabel(pFreePlate,'Text','Free Stream Temperature T_\infty (°C)','Position',tinfLblPos);
efTinf_free = uieditfield(pFreePlate,'numeric','Value',25,'Position',tinfEditPos);

% Direct Ra toggle
cbRaDirect = uicheckbox(pFree,'Text','Enter Rayleigh number directly','Position',[340 620 230 24],'Value',false);
uilabel(pFree,'Text','Ra (direct)','Position',[580 620 90 24]);
efRaDirect = uieditfield(pFree,'numeric','Value',1e7,'Position',[670 620 140 24],'Enable','off');
cbRaDirect.ValueChangedFcn = @(~,~) setEnable(efRaDirect, cbRaDirect.Value);

% Flat-plate options
uilabel(pFreePlate,'Text','Orientation','Position',[20 300 120 24]);
ddOrient = uidropdown(pFreePlate,'Items',{'Vertical','Horizontal: hot surface up','Horizontal: hot surface down','Inclined (cold up/hot down)'},...
    'Value','Vertical','Position',[120 300 260 24],'ValueChangedFcn',@(~,~)toggleIncl());
uilabel(pFreePlate,'Text','Length L (m)','Position',[20 260 140 24]);  efL_free = uieditfield(pFreePlate,'numeric','Value',0.5,'Position',[120 260 120 24]);
uilabel(pFreePlate,'Text','Inclination \theta (deg)','Position',[360 260 140 24]); efTheta = uieditfield(pFreePlate,'numeric','Value',30,'Position',[510 260 100 24]); efTheta.Enable='off';
uibutton(pFreePlate,'Text','Compute h (Flat Plate)','Position',[20 210 260 34],'ButtonPushedFcn',@(~,~)onFreePlate());

% Cylinder options
uilabel(pFreeCyl,'Text','Diameter D (m)','Position',[20 300 140 24]); efD_free_cyl = uieditfield(pFreeCyl,'numeric','Value',0.05,'Position',[120 300 120 24]);
uibutton(pFreeCyl,'Text','Compute h (Cylinder)','Position',[20 250 260 34],'ButtonPushedFcn',@(~,~)onFreeCyl());

% Sphere options
uilabel(pFreeSph,'Text','Diameter D (m)','Position',[20 300 140 24]); efD_free_sph = uieditfield(pFreeSph,'numeric','Value',0.05,'Position',[120 300 120 24]);
uibutton(pFreeSph,'Text','Compute h (Sphere)','Position',[20 250 260 34],'ButtonPushedFcn',@(~,~)onFreeSph());

%% ===== Free Convection (Enclosed) — Rectangular cavity only =====
uilabel(pEnc,'Text','Geometry:','Position',[20 620 100 24]);
ddEnc = uidropdown(pEnc,'Items',{'Rectangular cavity'},'Value','Rectangular cavity', ...
                   'Position',[110 620 220 24]);

pEncRect = uipanel(pEnc,'Title','Rectangular cavity','Position',[20 210 1060 400],'BorderType','line');

% Configuration selection
uilabel(pEncRect,'Text','Configuration:','Position',[20 340 120 24]);
ddEncMode = uidropdown(pEncRect,'Items',{'Horizontal: hot below','Vertical cavity'}, ...
                       'Value','Horizontal: hot below', 'Position',[140 340 220 24]);

% Temperatures and dimensions
uilabel(pEncRect,'Text','Hot surface temperature T_{hot} (°C)','Position',[20 300 260 24]);
efTh_enc = uieditfield(pEncRect,'numeric','Value',60,'Position',[290 300 100 24]);

uilabel(pEncRect,'Text','Cold surface temperature T_{cold} (°C)','Position',[20 260 260 24]);
efTc_enc = uieditfield(pEncRect,'numeric','Value',20,'Position',[290 260 100 24]);

uilabel(pEncRect,'Text','Gap between surfaces L (m)','Position',[20 220 220 24]);
efL_enc = uieditfield(pEncRect,'numeric','Value',0.05,'Position',[290 220 100 24]);

uilabel(pEncRect,'Text','Cavity height H (m)  (used for vertical only)','Position',[20 180 320 24]);
efH_enc = uieditfield(pEncRect,'numeric','Value',0.20,'Position',[350 180 100 24]);

uibutton(pEncRect,'Text','Compute h (Rectangular cavity)','Position',[20 130 260 34], ...
         'ButtonPushedFcn',@(~,~)onEncRect());

%% ===== Results =====
ta = uitextarea(pRes,'Position',[24 70 1052 576],'Editable','off','FontName','Consolas','FontSize',16);
ta.Value = "Results will appear here.";

btnResBack   = uibutton(pRes,'Text','← Back','FontSize',14, ...
    'Position',[24, 20, 120, 36], 'ButtonPushedFcn',@(~,~) resBack());
btnRestart   = uibutton(pRes,'Text','⟲ Restart','FontSize',14, ...
    'Position',[154, 20, 120, 36], 'ButtonPushedFcn',@(~,~) restartToMain());

%% ===== Navigation helpers & layout =====
function showPage(which)
    if ~strcmp(which,'results')
        lastNonResultsPage = which;
    end
    currentPage = which;

    pMain.Visible   = strcmp(which,'main');
    pForced.Visible = strcmp(which,'forced');
    pInt.Visible    = strcmp(which,'internal');
    pFree.Visible   = strcmp(which,'free');
    pEnc.Visible    = strcmp(which,'enclosed');
    pRes.Visible    = strcmp(which,'results');

    btnBack.Visible = iff(strcmp(which,'main') || strcmp(which,'results'),'off','on');

    relayout();
end

function goBack()
    if strcmp(currentPage,'results')
        showPage(lastNonResultsPage);
    else
        showPage('main');
    end
end

fig.SizeChangedFcn = @(~,~) relayout();
function relayout()
    pad = 20; fw  = fig.InnerPosition(3);
    btnBack.Position = [fw-160-pad, 16, 140, 36];

    % ----- Main page layout -----
    pmPos = pMain.InnerPosition; pw = pmPos(3); centerX = pmPos(1) + pw/2;
    h = 28; gap = 10; lblW = 300; fldW = 260;

    yLabel1 = 560; yCtrl1 = yLabel1 - (h+gap);
    lblProp.Position = [centerX - lblW/2, yLabel1, lblW, h];
    ddFluid.Position = [centerX - fldW/2, yCtrl1,  fldW, h];

    yLabel2 = yCtrl1 - 2*(h+gap); yCtrl2 = yLabel2 - (h+gap);
    lblTf.Position = [centerX - lblW/2, yLabel2, lblW, h];
    efTfC.Position = [centerX - fldW/2, yCtrl2,  fldW, h];

    titleY = yCtrl2 - 50; lblWbig = 540; lblH = 34;
    lblNextTitle.Position = [centerX - lblWbig/2, titleY, lblWbig, lblH];

    btnW = 420; btnH = 48; btnGap = 16;
    btnGoForced.Position   = [centerX - btnW/2, titleY - 60, btnW, btnH];
    btnGoInternal.Position = [centerX - btnW/2, titleY - 60 - (btnH+btnGap), btnW, btnH];
    btnGoFree.Position     = [centerX - btnW/2, titleY - 60 - 2*(btnH+btnGap), btnW, btnH];
    btnGoEnclosed.Position = [centerX - btnW/2, titleY - 60 - 3*(btnH+btnGap), btnW, btnH];

    % ----- Results page layout -----
    pr = pRes.InnerPosition;
    padR = 24;
    btnWr = 120; btnHr = 36; gapR = 10;

    btnResBack.Position = [padR, padR, btnWr, btnHr];
    btnRestart.Position = [pr(3) - padR - btnWr, padR, btnWr, btnHr];
    ta.Position = [padR, padR + btnHr + gapR, pr(3) - 2*padR, pr(4) - (padR + btnHr + gapR) - padR];
end

function resBack()
    showPage(lastNonResultsPage);
end

function restartToMain()
    ta.Value = "Results will appear here.";
    showPage('main');
end

showPage('main');

%% ===== Small UI helpers =====
function setEnable(ctrl, tf)
    if tf, ctrl.Enable='on'; else, ctrl.Enable='off'; end
end
function setReDirect(cb, efRe, others)
    setEnable(efRe, cb.Value);
    for c = others, c.Enable = iff(cb.Value,'off','on'); end
end
function out = iff(tf,a,b)
    if tf, out=a; else, out=b; end
end
function toggleMuS(evt,efTs,efMuS)
    if strcmp(evt.NewValue.Text,'Enter μ_s manually'), efMuS.Enable='on'; efTs.Enable='off';
    else, efMuS.Enable='off'; efTs.Enable='on'; end
end

function switchExt()
    switch ddExt.Value
        case 'Flat plate'
            showOnly(pFlat, {pCyl,pSph,pBank});
        case 'Cylinder'
            showOnly(pCyl, {pFlat,pSph,pBank});
        case 'Sphere'
            showOnly(pSph, {pFlat,pCyl,pBank});
        case 'Bank of tubes'
            showOnly(pBank,{pFlat,pCyl,pSph});
    end
end
function showOnly(onPanel, offPanels)
    onPanel.Visible='on'; onPanel.BorderType='line';
    for p = offPanels, p{1}.Visible='off'; p{1}.BorderType='none'; end
    uistack(onPanel,'top');
end

function switchInt()
    switch ddInt.Value
        case 'Circular — Laminar (fully developed, T_w const)'
            showOnly(pIntLam,{pIntTurb,pIntNon});
        case 'Circular — Turbulent (Dittus–Boelter)'
            showOnly(pIntTurb,{pIntLam,pIntNon});
        case 'Noncircular ducts'
            showOnly(pIntNon,{pIntLam,pIntTurb});
    end
end

function switchFree()
    plate = strcmp(ddFree.Value,'Flat plate'); cyl = strcmp(ddFree.Value,'Horizontal cylinder'); sph = strcmp(ddFree.Value,'Sphere');
    pFreePlate.Visible = plate;  pFreeCyl.Visible = cyl;  pFreeSph.Visible = sph;
    if plate
        [lblTs.Parent, efTs_free.Parent, lblTinf.Parent, efTinf_free.Parent] = deal(pFreePlate,pFreePlate,pFreePlate,pFreePlate);
    elseif cyl
        [lblTs.Parent, efTs_free.Parent, lblTinf.Parent, efTinf_free.Parent] = deal(pFreeCyl,pFreeCyl,pFreeCyl,pFreeCyl);
    else
        [lblTs.Parent, efTs_free.Parent, lblTinf.Parent, efTinf_free.Parent] = deal(pFreeSph,pFreeSph,pFreeSph,pFreeSph);
    end
    lblTs.Position = tsLblPos; efTs_free.Position = tsEditPos; lblTinf.Position = tinfLblPos; efTinf_free.Position = tinfEditPos;
end

function toggleIncl()
    if strcmp(ddOrient.Value,'Inclined (cold up/hot down)'), efTheta.Enable='on'; else, efTheta.Enable='off'; end
end

%% ===== Shared helpers=====
    function [fluidStr, T_K] = getFluidAndTempK()
        fluidStr = lower(string(ddFluid.Value));        % "air" or "water"
        T_K = efTfC.Value + 273.15;
    end

    % ---- now calls getFluidProp instead of getAirProp/getWaterProp ----
    function [nu, Pr, k, mu, rho, cp, alpha] = getProps(fluidStr, T_K)
        nu    = getFluidProp(fluidStr,T_K,'nu');
        Pr    = getFluidProp(fluidStr,T_K,'pr');
        k     = getFluidProp(fluidStr,T_K,'k');
        mu    = getFluidProp(fluidStr,T_K,'mu');
        rho   = getFluidProp(fluidStr,T_K,'rho');
        cp    = getFluidProp(fluidStr,T_K,'cp');  
        alpha = getFluidProp(fluidStr,T_K,'alpha');
    end

    % ---- same formatting, but uses getFluidProp ----
    function out = propBlock(fluidStr,T_K)
        rho   = getFluidProp(fluidStr,T_K,'rho');
        cp    = getFluidProp(fluidStr,T_K,'cp');
        mu    = getFluidProp(fluidStr,T_K,'mu');
        nu    = getFluidProp(fluidStr,T_K,'nu');
        alpha = getFluidProp(fluidStr,T_K,'alpha');

        fmtProps = [ ...
            '\n--- %s Properties at %.2f K ---\n' ...
            'rho:    %.6g kg/m³\n' ...
            'cp:     %.6g J/kg·K\n' ...
            'mu:     %.6g Pa·s\n' ...
            'nu:     %.6g m²/s\n' ...
            'alpha:  %.6g m²/s\n' ];
        out = sprintf(fmtProps, upper(char(fluidStr)), T_K, rho, cp, mu, nu, alpha);
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
        header = sprintf(fmtHeader, upper(char(fluidStr)), char(caseTitle), ...
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
        header = sprintf(fmtHeader, upper(char(fluidStr)), char(caseTitle), ...
            T_K, char(dimName), dimVal, Ra, Pr, Nu, k, h);
        lines = string(header) + newline + string(propBlock(fluidStr,T_K));
    end

%% ===== Forced External callbacks =====
    function onFlatPlate()
        [fluidStr, T_K] = getFluidAndTempK(); [nu, Pr, k, ~, ~, ~, ~] = getProps(fluidStr, T_K);
        if cbReFlat.Value
            Re_L = efReL_flat.Value; L = efL_flat.Value; V = Re_L*nu/max(L,eps);
        else
            V = efV_flat.Value; L = efL_flat.Value; Re_L = V*L/nu;
        end
        [Nu_L, convType] = ExtConvFlatPlate(Re_L, Pr);
        if isempty(Nu_L), uialert(fig,'Inputs outside correlation range for flat plate.','Warning'); return; end
        h = Nu_L * k / L;
        ta.Value = composeResults(fluidStr, ['Forced External → Flat Plate (' char(convType) ')'], T_K, 'L', L, V, Re_L, Pr, Nu_L, k, h);
        showPage('results');
    end

    function onCylinder()
        [fluidStr, T_K] = getFluidAndTempK(); [nu, Pr, k, ~, ~, ~, ~] = getProps(fluidStr, T_K);
        if cbReCyl.Value
            Re_D = efReD_cyl.Value; D = efD_cyl.Value; V = Re_D*nu/max(D,eps);
        else
            V = efV_cyl.Value; D = efD_cyl.Value; Re_D = V*D/nu;
        end
        Nu_D = ExtConvCyl(Re_D, Pr); if isempty(Nu_D), uialert(fig,'Inputs outside correlation range for cylinder.','Warning'); return; end
        h = Nu_D * k / D;
        ta.Value = composeResults(fluidStr, 'Forced External → Cylinder', T_K, 'D', D, V, Re_D, Pr, Nu_D, k, h);
        showPage('results');
    end

    function onSphere()
        [fluidStr, T_K] = getFluidAndTempK(); [nu, Pr, k, mu, ~, ~, ~] = getProps(fluidStr, T_K);
        if cbReSph.Value
            Re_D = efReD_sph.Value; D = efD_sph.Value; V = Re_D*nu/max(D,eps);
        else
            V = efV_sph.Value; D = efD_sph.Value; Re_D = V*D/nu;
        end
        if rbTs.Value
            Ts_K = efTsC.Value + 273.15;
            mu_s = getFluidProp(fluidStr,Ts_K,'mu');
        else
            mu_s = efMuS_sph.Value;
        end
        ratio = mu/mu_s;
        if ~(Re_D>=3.5 && Re_D<=7.6e4) || ~(Pr>=0.71 && Pr<=380) || ~(ratio>=1 && ratio<=3.2)
            uialert(fig, sprintf('Sphere correlation limits:\nRe_D=%.3g (3.5–7.6e4)\nPr=%.3g (0.71–380)\nμ/μ_s=%.3g (1.0–3.2)', Re_D, Pr, ratio),'Out of range');
        end
        Nu_D = ExtConvSphere(Re_D, Pr, mu, mu_s); if isempty(Nu_D), return; end
        h = Nu_D*k/D;
        lines = composeResults(fluidStr, 'Forced External → Sphere', T_K, 'D', D, V, Re_D, Pr, Nu_D, k, h);
        ta.Value = [lines; sprintf('\nμ/μ_s used: %.3g   (μ=%.3g Pa·s, μ_s=%.3g Pa·s)\n', ratio, mu, mu_s)];
        showPage('results');
    end

    function onTubeBank()
        [fluidStr, ~] = getFluidAndTempK();   % 'air' or 'water'

        % Temperatures from GUI (convert °C -> K)
        Ti = efTi.Value + 273.15;
        To = efTo.Value + 273.15;
        Ts = efTs.Value + 273.15;

        % Geometry & flow from GUI
        D   = efD_b.Value;
        ST  = efST.Value;
        SL  = efSL.Value;
        NL  = efNL.Value;
        NT  = efNT.Value;
        V   = efVb.Value;
        tubeType = iff(strcmp(ddArrange.Value,'Aligned'), 1, 2);  % 1=aligned, 2=staggered

        % Assume / tolerance from GUI
        assumedVar = 1;      % 1=None, 2=Outlet T_o, 3=Inlet T_i
        tol        = NaN;
        if strcmp(ddAssume.Value,'Outlet T_o')
            assumedVar = 2;
            tol = efTol.Value;
        elseif strcmp(ddAssume.Value,'Inlet T_i')
            assumedVar = 3;
            tol = efTol.Value;
        end

        % Call GUI version of tube-bank solver
        [Nu_D, Re_Dmax, V_max, Ti_out, To_out, C1, C2, m, Pr_s, ...
            q_p, DT_lm, iter, tolUsed, converged] = ...
            handleTubeBank_GUI(fluidStr, Ti, To, Ts, D, V, ST, SL, NL, NT, ...
                               tubeType, assumedVar, tol);

        % Final film temperature and properties
        Tf = 0.5 * (Ti_out + To_out);
        k  = getFluidProp(fluidStr, Tf, 'k');
        h  = Nu_D * k / D;

        % ---- Detailed text for the tube-bank tab (tbOut) ----
        if isnan(tolUsed)
            tolStr = 'n/a';
        else
            tolStr = num2str(tolUsed);
        end

        tbOut.Value = sprintf(['Tube Bank Details\n' ...
            '-----------------------------\n' ...
            'Arrangement:           %s\n' ...
            'D:                     %.4g m\n' ...
            'S_T:                   %.4g m\n' ...
            'S_L:                   %.4g m\n' ...
            'N_L (rows):            %d\n' ...
            'N_T (per row):         %d\n' ...
            'V_free:                %.4g m/s\n' ...
            'V_max:                 %.4g m/s\n' ...
            'Re_D,max:              %.4g\n' ...
            'Pr_s:                  %.4g\n' ...
            'C1, C2, m:             %.4g, %.4g, %.4g\n' ...
            'Nu_D:                  %.4g\n' ...
            'h:                     %.4g W/m^2-K\n' ...
            'T_i (inlet):           %.2f K\n' ...
            'T_o (outlet):          %.2f K\n' ...
            'T_s (surface):         %.2f K\n' ...
            'T_f (film):            %.2f K\n' ...
            'DeltaT_lm:             %.4g K\n' ...
            'q'' (per unit length): %.4g W/m\n' ...
            'Iterations:            %d\n' ...
            'Converged:             %s\n' ...
            'Final tolerance:       %s %%\n'], ...
            ddArrange.Value, D, ST, SL, NL, NT, V, V_max, Re_Dmax, Pr_s, ...
            C1, C2, m, Nu_D, h, Ti_out, To_out, Ts, Tf, DT_lm, q_p, ...
            iter, string(converged), tolStr);

        % ---- Main Results page text (ta) in same style as the sphere ----
        header = sprintf(['=== Convection Coefficient Results ===\n' ...
            'Fluid:                             %s\n' ...
            'Case:                              Forced External -> Tube Bank\n' ...
            'Film Temperature (T_f):            %.2f K\n' ...
            'D:                                 %.6g m\n' ...
            'Re_D,max:                          %.6g\n' ...
            'Nu_D:                              %.6g\n' ...
            'k (at T_f):                        %.6g W/m-K\n' ...
            'h:                                 %.6g W/m^2-K\n' ...
            'DeltaT_lm:                         %.6g K\n' ...
            'q'' (per unit length):             %.6g W/m\n'], ...
            upper(char(fluidStr)), Tf, D, Re_Dmax, Nu_D, k, h, DT_lm, q_p);

        ta.Value = string(header) + newline + string(propBlock(fluidStr, Tf));
        showPage('results');
    end
%% ===== Forced Internal callbacks =====
    function onIntLam()
        [fluidStr, T_K] = getFluidAndTempK(); [nu, Pr, k, ~, ~, ~, ~] = getProps(fluidStr, T_K);
        if cbReIntL.Value
            Re = efRe_intL.Value; D = efD_lam.Value; V = Re*nu/max(D,eps);
        else
            V = efV_lam.Value; D = efD_lam.Value; Re = V*D/nu;
        end
        if Re >= 2300, uialert(fig,sprintf('Warning: Re=%.3g ≥ 2300, not laminar (Nu=3.66 is FD laminar).',Re),'Warning'); end
        Nu = 3.66; h  = Nu * k / D;
        ta.Value = composeResults(fluidStr, 'Forced Internal → Circular tube (Laminar FD, T_w const)', T_K, 'D', D, V, Re, Pr, Nu, k, h);
        showPage('results');
    end

    function onIntTurb()
        [fluidStr, T_K] = getFluidAndTempK(); [nu, Pr, k, ~, ~, ~, ~] = getProps(fluidStr, T_K);
        if cbReIntT.Value
            Re = efRe_intT.Value; D = efD_turb.Value; V = Re*nu/max(D,eps);
        else
            V = efV_turb.Value; D = efD_turb.Value; Re = V*D/nu;
        end
        n = 0.4; if contains(ddDir.Value,'Cooling'), n=0.3; end
        if Re < 1e4, uialert(fig,sprintf('Warning: Re=%.3g < 10,000 (Dittus–Boelter range).',Re),'Warning'); end
        Nu = 0.023 * Re^0.8 * Pr^n; h  = Nu * k / D;
        ta.Value = composeResults(fluidStr, ['Forced Internal → Circular tube (Dittus–Boelter, n=' num2str(n) ')'], T_K, 'D', D, V, Re, Pr, Nu, k, h);
        showPage('results');
    end

       % ---- Noncircular ducts callback (GUI; hydraulic diameter method) ----
    function onIntNon()
        [fluidStr, ~] = getFluidAndTempK();

        % Temperatures and length
        TmK = efTm_non.Value + 273.15;
        TsK = efTs_non.Value + 273.15; 
        L   = efL_non.Value;
        if ~(L > 0)
            uialert(fig,'Duct length L must be positive.','Input error');
            return;
        end

        % Geometry 
        a = efA_non.Value;
        b = efB_non.Value;
        if ~(a > 0) || ~(b > 0)
            uialert(fig,'Dimensions a and b must be positive.','Input error');
            return;
        end
        A_c = a*b;
        P   = 2*(a + b);
        D_h = 4*A_c / max(P,eps);

        % Fluid properties at T_m
        [nu, Pr, k, ~, ~, ~, ~] = getProps(fluidStr, TmK);

        % Reynolds number / velocity
        if cbReNon.Value
            Re_D = efRe_non.Value;
            V    = Re_D * nu / max(D_h,eps);
        else
            V    = efV_non.Value;
            Re_D = V * D_h / nu;
        end

        % Correlations using hydraulic diameter
        if Re_D < 2300
            Nu_D = 3.66;   % laminar, fully developed, const T_w (approx.)
            regime = 'Laminar FD, constant T_w (Nu=3.66 with D_h)';
        else
            n = 0.4;       % assume heating
            Nu_D = 0.023 * Re_D^0.8 * Pr^n;
            regime = sprintf('Turbulent, Dittus–Boelter (n=%.1f) with D_h',n);
        end

        h = Nu_D * k / D_h;

        caseTitle = sprintf('Forced Internal → Noncircular duct (rectangular, a=%.4g m, b=%.4g m)',a,b);
        lines = composeResults(fluidStr, caseTitle, TmK, 'D_h', D_h, V, Re_D, Pr, Nu_D, k, h);

        extra = sprintf('\nInputs:\n   T_m = %.2f K   T_s = %.2f K   L = %.4g m\n   a = %.4g m   b = %.4g m   D_h = %.4g m\nModel: %s\n', ...
                        TmK, TsK, L, a, b, D_h, regime);

        ta.Value = [lines; string(extra)];
        showPage('results');
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
        if cbRaDirect.Value
            Ra = efRaDirect.Value; TfK = (efTs_free.Value + efTinf_free.Value)/2 + 273.15;
            [~, Pr, k, ~, ~, ~, ~] = getProps(fluidStr, TfK); L = max(efL_free.Value,eps);
        else
            TsC = efTs_free.Value; TinfC = efTinf_free.Value; L = efL_free.Value;
            [Ra, Pr, k, TfK] = freeRa(fluidStr, TsC, TinfC, L);
        end
        switch ddOrient.Value
            case 'Vertical', Nu = freeConvExtFlatPlateVert(Ra,Pr); caseTitle = 'Free Convection → Flat Plate (Vertical)';
            case 'Horizontal: hot surface up',   Nu = freeConvExtPlateHorizHotUpper(Ra,Pr); caseTitle = 'Free Convection → Flat Plate (Horizontal, hot up)';
            case 'Horizontal: hot surface down', Nu = freeConvExtPlateHorizHotLower(Ra,Pr); caseTitle = 'Free Convection → Flat Plate (Horizontal, hot down)';
            case 'Inclined (cold up/hot down)'
                theta = efTheta.Value; Nu = freeConvExtPlateInc(Ra,Pr,theta); caseTitle = sprintf('Free Convection → Flat Plate (Inclined, \\theta=%.1f°)',theta);
        end
        if isempty(Nu), uialert(fig,'Inputs outside correlation range for flat plate.','Warning'); return; end
        h = Nu * k / L;
        ta.Value = composeResultsFree(fluidStr, caseTitle, TfK, 'L', L, Ra, Pr, Nu, k, h);
        showPage('results');
    end

    function onFreeCyl()
        [fluidStr, ~] = getFluidAndTempK();
        if cbRaDirect.Value
            Ra = efRaDirect.Value; TfK = (efTs_free.Value + efTinf_free.Value)/2 + 273.15; [~, Pr, k, ~, ~, ~, ~] = getProps(fluidStr,TfK); D = max(efD_free_cyl.Value,eps);
        else
            TsC = efTs_free.Value; TinfC = efTinf_free.Value; D = efD_free_cyl.Value;
            [Ra, Pr, k, TfK] = freeRa(fluidStr, TsC, TinfC, D);
        end
        Nu = freeConvExtHorizCyl(Ra,Pr); if isempty(Nu), uialert(fig,'Inputs outside correlation range (cylinder).','Warning'); return; end
        h = Nu * k / D;
        ta.Value = composeResultsFree(fluidStr, 'Free Convection → Horizontal Cylinder', TfK, 'D', D, Ra, Pr, Nu, k, h);
        showPage('results');
    end

    function onFreeSph()
        [fluidStr, ~] = getFluidAndTempK();
        if cbRaDirect.Value
            Ra = efRaDirect.Value; TfK = (efTs_free.Value + efTinf_free.Value)/2 + 273.15; [~, Pr, k, ~, ~, ~, ~] = getProps(fluidStr,TfK); D = max(efD_free_sph.Value,eps);
        else
            TsC = efTs_free.Value; TinfC = efTinf_free.Value; D = efD_free_sph.Value;
            [Ra, Pr, k, TfK] = freeRa(fluidStr, TsC, TinfC, D);
        end
        Nu = freeConvExtSphere(Ra,Pr); if isempty(Nu), uialert(fig,'Inputs outside correlation range (sphere).','Warning'); return; end
        h = Nu * k / D;
        ta.Value = composeResultsFree(fluidStr, 'Free Convection → Sphere', TfK, 'D', D, Ra, Pr, Nu, k, h);
        showPage('results');
    end

%% ===== Free Convection (Enclosed) callback =====
    function onEncRect()
        fluidStr = lower(string(ddFluid.Value));

        T_hot_C = efTh_enc.Value;
        T_cold_C = efTc_enc.Value;
        L = efL_enc.Value;
        H = efH_enc.Value;

        if ~(L > 0) || ~(H > 0)
            uialert(fig,'L and H must be positive.','Input error');
            return;
        end

        T_hot_K = T_hot_C + 273.15;
        T_cold_K = T_cold_C + 273.15;

        % Use your getRa to get Ra_L and T_f 
        [Ra_L, ~, ~, T_f, ~] = getRa(fluidStr, 'L', T_hot_K, T_cold_K, L);

        Pr = getFluidProp(fluidStr, T_f, 'pr');
        k  = getFluidProp(fluidStr, T_f, 'k');

        if strcmp(ddEncMode.Value,'Horizontal: hot below')
            Nu_L = freeConvEncRectCavHorizHeatBelow(Ra_L,Pr);
            caseTitle = 'Free Convection (Enclosed) → Rectangular cavity, horizontal (hot below)';
        else
            Nu_L = freeConvEncRectVertCav(Ra_L,Pr,H,L);
            caseTitle = sprintf('Free Convection (Enclosed) → Rectangular cavity, vertical (H/L=%.3g)',H/L);
        end

        h = Nu_L * k / L;

        lines = composeResultsFree(fluidStr, caseTitle, T_f, 'L', L, Ra_L, Pr, Nu_L, k, h);
        if strcmp(ddEncMode.Value,'Vertical cavity')
            extra = sprintf('\nAdditional geometry:\n   H = %.4g m   L = %.4g m   H/L = %.3g\n',H,L,H/L);
            lines = [lines; string(extra)];
        end

        ta.Value = lines;
        showPage('results');
    end
    function [Nu_D, Re_Dmax, V_max, Ti_out, To_out, C1, C2, m, Pr_s, ...
          q_p, DT_lm, iter, tolUsed, converged] = ...
          handleTubeBank_GUI(fluid, Ti, To, Ts, D, V_free, S_T, S_L, ...
                             N_L, N_T, tubeType, assumedVar, tol)
    % GUI version of handleTubeBank

    Ti_k = Ti;
    To_k = To;
    Ts_k = Ts;

    maxIter   = 50;
    iter      = 0;
    converged = false;
    tolUsed   = NaN;

    % If no assumed temperature, just do a single pass (no iteration)
    if assumedVar == 1 || isnan(tol)
        [Re_Dmax, V_max] = getReTubeBank(fluid, D, V_free, Ti_k, To_k, S_T, S_L, tubeType);
        [Nu_D, C1, C2, m, Pr_s] = ExtConvTubeBank(fluid, Re_Dmax, Ti_k, To_k, Ts_k, S_T, S_L, N_L, tubeType);

        Tf = 0.5 * (Ti_k + To_k);
        k  = getFluidProp(fluid, Tf, 'k');
        h  = Nu_D * k / D;

        N   = N_L * N_T;
        rho = getFluidProp(fluid, Tf, 'rho');
        cp  = getFluidProp(fluid, Tf, 'cp');

        % exponent, same as original code
        expo = exp(-pi * D * N * h / (rho * V_free * N_T * S_T * cp)); %#ok<NASGU>

        DT_lm = ((Ts_k - Ti_k) - (Ts_k - To_k)) / ...
                log((Ts_k - Ti_k) / (Ts_k - To_k));
        q_p   = N_T * N_L * (h * pi * D * DT_lm);

        Ti_out   = Ti_k;
        To_out   = To_k;
        iter     = 1;
        converged = true;
        return;
    end

    % ----- Iterative solution when T_i or T_o is assumed -----
    while ~converged && iter < maxIter
        iter = iter + 1;

        % Re and Nu based on current Ti_k, To_k
        [Re_Dmax, V_max] = getReTubeBank(fluid, D, V_free, Ti_k, To_k, S_T, S_L, tubeType);
        [Nu_D, C1, C2, m, Pr_s] = ExtConvTubeBank(fluid, Re_Dmax, Ti_k, To_k, Ts_k, S_T, S_L, N_L, tubeType);

        Tf = 0.5 * (Ti_k + To_k);
        k  = getFluidProp(fluid, Tf, 'k');
        h  = Nu_D * k / D;

        N   = N_L * N_T;
        rho = getFluidProp(fluid, Tf, 'rho');
        cp  = getFluidProp(fluid, Tf, 'cp');

        expo = exp(-pi * D * N * h / (rho * V_free * N_T * S_T * cp));

        switch assumedVar
            case 2   % Outlet temperature T_o was assumed
                T_o_calc = Ts_k - (Ts_k - Ti_k) * expo;
                tolUsed  = abs((T_o_calc - To_k) / To_k) * 100;

                if tolUsed <= tol
                    converged = true;
                else
                    To_k = To_k + 0.5 * (T_o_calc - To_k);
                end

            case 3   % Inlet temperature T_i was assumed
                T_i_calc = Ts_k - (Ts_k - To_k) / expo;
                tolUsed  = abs((T_i_calc - Ti_k) / Ti_k) * 100;

                if tolUsed <= tol
                    converged = true;
      
                else
                    Ti_k = Ti_k + 0.5 * (T_i_calc - Ti_k);
                end

            otherwise
                converged = true;
        end
    end

    % Final temperatures / outputs
    Ti_out = Ti_k;
    To_out = To_k;

    Tf = 0.5 * (Ti_k + To_k);
    k  = getFluidProp(fluid, Tf, 'k');
    h  = Nu_D * k / D;

    DT_lm = ((Ts_k - Ti_k) - (Ts_k - To_k)) / ...
            log((Ts_k - Ti_k) / (Ts_k - To_k));
    q_p   = N_T * N_L * (h * pi * D * DT_lm);
    end

end % ===== end main GUI function =====

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
    
    T_f = 0.5*(T_i + T_o);

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
% val = getFluidProp(fluid, T_K, prop)
%
% Inputs:
%   fluid = 'air' or 'water'
%   T_K   = temperature (K)
%   prop  = 'rho','mu','cp','k','nu','pr','alpha'
%
% Output:
%   val   = interpolated property value

    fluid = lower(string(fluid));
    prop  = lower(string(prop));

    % --- Handle alpha here so we don't ask air/water tables for it ---
    if prop == "alpha"
        switch fluid
            case "air"
                rho = getAirProp(T_K, "rho");
                cp  = getAirProp(T_K, "cp");   % J/kg·K
                k   = getAirProp(T_K, "k");    % W/m·K
            case "water"
                rho = getWaterProp(T_K, "rho");
                cp  = getWaterProp(T_K, "cp");
                k   = getWaterProp(T_K, "k");
            otherwise
                error('Fluid not recognized. Use ''air'' or ''water''.');
        end

        % thermal diffusivity α = k / (ρ cp)
        val = k ./ (rho .* cp);
        return;
    end

    % --- All other properties are directly from the tables ---
    switch fluid
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