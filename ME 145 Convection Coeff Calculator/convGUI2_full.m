function convGUI2_full
% Convection Coefficient Calculator — Full GUI
% - Matches original CLI coverage with GUI panels
% - Adds Tube Bank (iterative) + results
% - Provides structure for noncircular/annulus/enclosed (as in original)

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

% ---- Tube bank panel (NEW)
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
                                 'Noncircular ducts (placeholder)', ...
                                 'Concentric tube annulus (placeholder)'}, ...
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

% ---- Noncircular ducts (placeholder panel like CLI)
pIntNon = uipanel(pInt,'Title','Noncircular ducts (placeholder)','Position',[20 240 1060 360],'BorderType','none','Visible','off');
uilabel(pIntNon,'Text','This branch was not implemented in the CLI. Add correlations for square/rectangular/triangular ducts as needed.',...
    'Position',[20 280 1000 40]);

% ---- Concentric tube annulus (placeholder)
pIntAnn = uipanel(pInt,'Title','Concentric tube annulus (placeholder)','Position',[20 240 1060 360],'BorderType','none','Visible','off');
uilabel(pIntAnn,'Text','This branch was not implemented in the CLI. Add annulus correlations and hydraulic diameter handling.',...
    'Position',[20 280 1000 40]);

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

% Direct Ra toggle (like CLI)
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

%% ===== Free Convection (Enclosed) — structure matches CLI (placeholders) =====
uilabel(pEnc,'Text','Geometry:','Position',[20 620 100 24]);
ddEnc = uidropdown(pEnc,'Items',{'Rectangular cavity','Concentric cylinders','Concentric spheres'},'Value','Rectangular cavity','Position',[110 620 220 24], 'ValueChangedFcn',@(~,~)switchEnc());

pEncRect = uipanel(pEnc,'Title','Rectangular cavity (placeholder)','Position',[20 210 1060 400],'BorderType','line');
pEncCyl  = uipanel(pEnc,'Title','Concentric cylinders (placeholder)','Position',[20 210 1060 400],'BorderType','none','Visible','off');
pEncSph  = uipanel(pEnc,'Title','Concentric spheres (placeholder)','Position',[20 210 1060 400],'BorderType','none','Visible','off');

uilabel(pEncRect,'Text','Not implemented in CLI; add correlations if desired.','Position',[20 330 600 24]);
uilabel(pEncCyl,'Text','Not implemented in CLI; add correlations if desired.','Position',[20 330 600 24]);
uilabel(pEncSph,'Text','Not implemented in CLI; add correlations if desired.','Position',[20 330 600 24]);

%% ===== Results =====
ta = uitextarea(pRes,'Position',[24 70 1052 576],'Editable','off','FontName','Consolas','FontSize',16);
ta.Value = "Results will appear here.";

btnResBack   = uibutton(pRes,'Text','← Back','FontSize',14, ...
    'Position',[24, 20, 120, 36], 'ButtonPushedFcn',@(~,~) resBack());
btnRestart   = uibutton(pRes,'Text','⟲ Restart','FontSize',14, ...
    'Position',[154, 20, 120, 36], 'ButtonPushedFcn',@(~,~) restartToMain());

%% ===== Navigation helpers & layout =====
function showPage(which)
    % Remember where we were before Results
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

    % Hide the top global back button on main AND results
    btnBack.Visible = iff(strcmp(which,'main') || strcmp(which,'results'),'off','on');

    relayout();
end

function goBack()
    % Global top-left Back behavior
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
pr = pRes.InnerPosition;  % [x y w h]
padR = 24;
btnWr = 120; btnHr = 36; gapR = 10;

% Back stays bottom-left
btnResBack.Position = [padR, padR, btnWr, btnHr];

% Restart moves to bottom-right
btnRestart.Position = [pr(3) - padR - btnWr, padR, btnWr, btnHr];

% Text area fills the space above buttons
ta.Position = [padR, padR + btnHr + gapR, pr(3) - 2*padR, pr(4) - (padR + btnHr + gapR) - padR];
end

function resBack()
    % Go back to the page we were on right before Results
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
            showOnly(pIntLam,{pIntTurb,pIntNon,pIntAnn});
        case 'Circular — Turbulent (Dittus–Boelter)'
            showOnly(pIntTurb,{pIntLam,pIntNon,pIntAnn});
        case 'Noncircular ducts (placeholder)'
            showOnly(pIntNon,{pIntLam,pIntTurb,pIntAnn});
        case 'Concentric tube annulus (placeholder)'
            showOnly(pIntAnn,{pIntLam,pIntTurb,pIntNon});
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

function switchEnc()
    switch ddEnc.Value
        case 'Rectangular cavity'
            showOnly(pEncRect,{pEncCyl,pEncSph});
        case 'Concentric cylinders'
            showOnly(pEncCyl,{pEncRect,pEncSph});
        case 'Concentric spheres'
            showOnly(pEncSph,{pEncRect,pEncCyl});
    end
end

%% ===== Shared helpers =====
    function [fluidStr, T_K] = getFluidAndTempK()
        fluidStr = lower(string(ddFluid.Value));
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
            otherwise, error('Unknown fluid');
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
            Ts_K = efTsC.Value + 273.15; mu_s = iff(fluidStr=="air", getAirProp(Ts_K,'mu'), getWaterProp(Ts_K,'mu'));
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
        [fluidStr, T_K] = getFluidAndTempK(); %#ok<NASGU>
        % Gather inputs (K)
        Ti = efTi.Value + 273.15; To = efTo.Value + 273.15; Ts = efTs.Value + 273.15;
        D  = efD_b.Value; ST = efST.Value; SL = efSL.Value; NL = efNL.Value; NT = efNT.Value;
        V  = efVb.Value;
        tubeType = iff(strcmp(ddArrange.Value,'Aligned'),1,2);
        assumedVar = 1; tol = NaN;
        if strcmp(ddAssume.Value,'Outlet T_o'), assumedVar = 2; tol = efTol.Value;
        elseif strcmp(ddAssume.Value,'Inlet T_i'), assumedVar = 3; tol = efTol.Value; end

        % Run the original logic
        [Nu_D, Re_Dmax, V_max, Ti_out, To_out, ~, ~, ~, ~, ...
            ~, ~, C1, C2, m, Pr_s, ~, q_p, DT_lm, iter, tolUsed, converged] = ...
            handleTubeBank_GUI(fluidStr, Ti, To, Ts, D, V, ST, SL, NL, NT, tubeType, assumedVar, tol);

        Tf = mean([Ti_out, To_out]); k = getFluidProp(fluidStr,Tf,'k'); h = Nu_D*k/D;

        tbOut.Value = sprintf(['Tube Bank Results\n' ...
            '-----------------\n' ...
            'Arrangement: %s   ST=%.4g m   SL=%.4g m   D=%.4g m\n' ...
            'N_L=%d   N_T=%d   V_{free}=%.4g m/s   V_{max}=%.4g m/s\n' ...
            'Re_{D,max}=%.4g   Pr_s=%.4g   C1=%.4g   C2=%.4g   m=%.4g\n' ...
            'Nu_D=%.4g   h=%.4g W/m^2·K\n' ...
            'T_i=%.2f K   T_o=%.2f K   T_s=%.2f K   T_f=%.2f K\n' ...
            '\x0394T_{lm}=%.4g K   q''=%.4g W/m\n' ...
            'Converged: %s (iter=%d, tol=%s%%)\n'], ...
            ddArrange.Value, ST, SL, D, NL, NT, V, V_max, Re_Dmax, Pr_s, C1, C2, m, ...
            Nu_D, h, Ti_out, To_out, Ts, Tf, DT_lm, q_p, string(converged), iter, num2str(tolUsed));

        % Compact summary into Results page:
        ta.Value = sprintf('Forced External → Tube Bank\nNu_D=%.4g, h=%.4g W/m^2·K, Re_{D,max}=%.4g, V_{max}=%.4g m/s, q''=%.4g W/m, \x0394T_{lm}=%.4g K', ...
                           Nu_D, h, Re_Dmax, V_max, q_p, DT_lm);
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
            [~, Pr, k, ~] = getProps(fluidStr, TfK); L = max(efL_free.Value,eps);
        else
            TsC = efTs_free.Value; TinfC = efTinf_free.Value; L = efL_free.Value;
            [Ra, Pr, k, TfK] = freeRa(fluidStr, TsC, TinfC, L);
        end
        switch ddOrient.Value
            case 'Vertical', Nu = freeConvExtFlatPlateVert(Ra,Pr); caseTitle = 'Free Convection → Flat Plate (Vertical)';
            case 'Horizontal: hot surface up',   Nu = freeConvExtPlateHorizHotUpper(Ra,Pr); caseTitle = 'Free Convection → Flat Plate (Horizontal, hot up)';
            case 'Horizontal: hot surface down', Nu = freeConvExtPlateHorizHotLower(Ra,Pr); caseTitle = 'Free Convection → Flat Plate (Horizontal, hot down)';
            case 'Inclined (cold up/hot down)'
                theta = efTheta.Value; Nu = freeConvExtPlateInc(Ra,Pr,theta); caseTitle = sprintf('Free Convection → Flat Plate (Inclined, \x03B8=%.1f°)',theta);
        end
        if isempty(Nu), uialert(fig,'Inputs outside correlation range for flat plate.','Warning'); return; end
        h = Nu * k / L;
        ta.Value = composeResultsFree(fluidStr, caseTitle, TfK, 'L', L, Ra, Pr, Nu, k, h);
        showPage('results');
    end

    function onFreeCyl()
        [fluidStr, ~] = getFluidAndTempK();
        if cbRaDirect.Value
            Ra = efRaDirect.Value; TfK = (efTs_free.Value + efTinf_free.Value)/2 + 273.15; [~, Pr, k, ~] = getProps(fluidStr,TfK); D = max(efD_free_cyl.Value,eps);
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
            Ra = efRaDirect.Value; TfK = (efTs_free.Value + efTinf_free.Value)/2 + 273.15; [~, Pr, k, ~] = getProps(fluidStr,TfK); D = max(efD_free_sph.Value,eps);
        else
            TsC = efTs_free.Value; TinfC = efTinf_free.Value; D = efD_free_sph.Value;
            [Ra, Pr, k, TfK] = freeRa(fluidStr, TsC, TinfC, D);
        end
        Nu = freeConvExtSphere(Ra,Pr); if isempty(Nu), uialert(fig,'Inputs outside correlation range (sphere).','Warning'); return; end
        h = Nu * k / D;
        ta.Value = composeResultsFree(fluidStr, 'Free Convection → Sphere', TfK, 'D', D, Ra, Pr, Nu, k, h);
        showPage('results');
    end

end % ===== end main GUI function =====

%% ===== Correlation functions (unchanged names/logic) =====
function [Nu_L, convType] = ExtConvFlatPlate(Re_L, Pr)
if Re_L <= 5e5 && Pr >= 0.6
    convType = 'laminar'; Nu_L = 0.664*Re_L.^0.5 .* Pr.^(1/3);
elseif Re_L >= 5e5 && Re_L <= 1e8 && Pr >= 0.6 && Pr <= 60
    convType = 'mixed';   Nu_L = (0.037*Re_L.^0.8 - 871) .* Pr.^(1/3);
else, warning('Re_L or Pr outside acceptable range.'); Nu_L = []; convType = 'RE_OUTSIDE';
end
end
function Nu_D = ExtConvCyl(Re_D, Pr)
if Re_D*Pr >= 0.2
    Nu_D = 0.3 + (0.62*Re_D.^0.5 .* Pr.^(1/3) .* (1+(0.4./Pr).^(2/3)).^(-1/4)) .* (1 + (Re_D/282000).^(5/8)).^(4/5);
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
else, Nu_L = []; disp('freeConvExtFlatPlateVert: Ra_L outside range');
end
end
function Nu_L = freeConvExtPlateInc(Ra_L,Pr,theta)
if theta>=0 && theta<=60, Nu_L = freeConvExtFlatPlateVert(Ra_L, Pr);
else, disp('theta outside range for inclined plate.'); Nu_L=[];
end
end
function Nu_L = freeConvExtPlateHorizHotUpper(Ra_L,Pr)
if Ra_L >= 1e4 && Ra_L <= 1e7 && Pr >= 0.7, Nu_L = 0.54*Ra_L.^(1/4);
elseif Ra_L >= 1e7 && Ra_L <= 1e11, Nu_L = 0.15*Ra_L.^(1/3);
else, disp('Invalid Ra_L/Pr (hot up)'); Nu_L=[];
end
end
function Nu_L = freeConvExtPlateHorizHotLower(Ra_L,Pr)
if Ra_L >= 1e4 && Ra_L <= 1e9 && Pr >= 0.7, Nu_L = 0.52*Ra_L.^(1/5);
else, disp('Invalid Ra_L/Pr (hot down)'); Nu_L=[];
end
end
function Nu_D = freeConvExtHorizCyl(Ra_D,Pr)
if Ra_D <= 1e12, Nu_D = (0.60 + 0.387*Ra_D.^(1/6) ./(1+(0.559/Pr).^(9/16)).^(8/27)).^2;
else, disp('Invalid Ra_D'); Nu_D=[];
end
end
function Nu_D = freeConvExtSphere(Ra_D,Pr)
if Ra_D <= 1e12 && Pr >= 0.7, Nu_D = 2 + 0.589*Ra_D.^(1/4) ./(1+(0.469/Pr).^(9/16)).^(4/9);
else, disp('Invalid Ra_D/Pr'); Nu_D=[];
end
end

%% ===== Tube bank logic (GUI-friendly wrapper around your original) =====
function [Nu_D, Re_Dmax, V_max, T_i, T_o, T_s, D, V, S_T, S_L, ...
          N_L, N_T, C1, C2, m, Pr_s, tubeType, q_p, DT_lm, iter, tol, converged] = ...
          handleTubeBank_GUI(fluid, T_i, T_o, T_s, D, V, S_T, S_L, N_L, N_T, tubeType, assumedVar, tol)

% Defaults for iteration like your CLI
if isnan(tol), tol = 1; end
iter = 0; maxIter = 50; converged = false;

while ~converged && iter < maxIter
    iter = iter + 1;

    % Re and constants
    [Re_Dmax, V_max] = getReTubeBank(fluid, D, V, T_i, T_o, S_T, S_L, tubeType);
    [Nu_D, C1, C2, m, Pr_s] = ExtConvTubeBank(fluid, Re_Dmax, T_i, T_o, T_s, S_T, S_L, N_L, tubeType);

    % h
    T_f = mean([T_i T_o]); k = getFluidProp(fluid, T_f, 'k'); h = Nu_D * k / D;
    N = N_T*N_L;

    % Properties
    rho = getFluidProp(fluid, T_f, 'rho'); cp  = getFluidProp(fluid, T_f, 'cp');

    % Exponential factor
    expo = exp(-pi*D*N*h / (rho*V*N_T*S_T*cp));

    switch assumedVar
        case 2  % assumed T_o
            T_o_calc = T_s - (T_s - T_i) * expo;
            err = abs((T_o_calc - T_o) / max(abs(T_o),1)) * 100;
            if err <= tol, converged = true; else, T_o = T_o + 0.5*(T_o_calc - T_o); end

        case 3  % assumed T_i
            T_i_calc = T_s - (T_s - T_o) / max(expo,1e-12);
            err = abs((T_i_calc - T_i) / max(abs(T_i),1)) * 100;
            if err <= tol, converged = true; else, T_i = T_i + 0.5*(T_i_calc - T_i); end

        otherwise
            converged = true; % no iteration needed
    end
end

% After convergence
DT_lm = ((T_s - T_i) - (T_s - T_o)) / log(max((T_s - T_i)/(T_s - T_o),1e-12));
q_p = N_T * (h * pi * D * DT_lm);
end

function [Re_Dmax, V_max] = getReTubeBank(fluid, D, V, T_i, T_o, S_T, S_L, tubeType)
fluid = lower(string(fluid));
T_f = 0.5*(T_i + T_o);
nu = getFluidProp(fluid, T_f, 'nu');
if tubeType == 1 % aligned
    V_max = V*S_T/(S_T-D);
else             % staggered
    V_max = V*S_T*S_L/((S_T-D)*(0.5*S_L));
end
Re_Dmax = V_max * D / nu;
end

function [Nu_D, C1, C2, m, Pr_s] = ExtConvTubeBank(fluid, Re_Dmax, T_i, T_o, T_s, S_T, S_L, N_L, tubeType)
T_f = mean([T_i, T_o]); Pr = getFluidProp(fluid, T_f, 'pr'); Pr_s = getFluidProp(fluid, T_s, 'pr');
ST_SL = S_T / S_L;

if tubeType == 1  % Aligned
    if Re_Dmax < 1e3, C1 = 0.8; m = 0.4;
    elseif Re_Dmax < 2e5, C1 = 0.27; m = 0.63;
    else, C1 = 0.021; m = 0.84; end
else % Staggered
    if Re_Dmax < 1e3, C1 = 0.9; m = 0.4;
    elseif Re_Dmax < 2e5
        if ST_SL < 2, C1 = 0.35*(ST_SL)^(1/5); m = 0.6; else, C1 = 0.4; m = 0.6; end
    else, C1 = 0.022; m = 0.84; end
end

NL_table = [1 2 3 4 5 7 10 13 16];
if tubeType == 1, C2_table = [0.70 0.80 0.86 0.90 0.92 0.95 0.97 0.98 0.99];
else,              C2_table = [0.64 0.76 0.84 0.89 0.92 0.95 0.97 0.98 0.99];
end
if N_L < 1, C2 = NaN; warning('N_L < 1');
elseif N_L < 20, C2 = interp1(NL_table, C2_table, N_L, 'linear');
else, C2 = 1; end

if Re_Dmax >= 10 && Re_Dmax <= 2e6 && Pr >= 0.68 && Pr <= 500
    Nu_D = C1*C2*Re_Dmax^m*Pr^0.36*(Pr/Pr_s)^(1/4);
else
    warning('Re_Dmax/Pr outside table range; computing anyway.');
    Nu_D = C1*C2*Re_Dmax^m*Pr^0.36*(Pr/Pr_s)^(1/4);
end
end

%% ===== Property tables (same as your originals) =====
function val = getFluidProp(fluid, T_K, prop)
switch lower(fluid)
    case 'air',   val = getAirProp(T_K, prop);
    case 'water', val = getWaterProp(T_K, prop);
    otherwise, error('Fluid not recognized.');
end
end

function val = getAirProp(T_K, prop)
if T_K < 100 || T_K > 3000, warning('Air T outside table.'); val = []; return; end
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
      1.417 1.478 1.558 1.665 2.726]*1000;
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
    otherwise, error('Property not recognized.');
end
end

function val = getWaterProp(T_K, prop)
if T_K < 273.15 || T_K > 647.3, warning('Water T outside table.'); val = []; return; end
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
 2.007 0.941 0.563 0.000]*1000;
mu_water_1e6 = [1750 1654 1422 1225 1080 959 855 765 695 639 577 532 489 453 420 ...
 389 364 343 326 306 289 279 264 237 215 162 138 123 112 104 ...
 95 91 90 81 72 64 45];
k_water_1e3 = [561 574 582 590 598 606 613 619 626 632 640 645 650 656 660 ...
 664 671 679 686 677 688 682 688 688 688 682 673 667 651 628 ...
 594 548 515 497 444 330 238]*1e-3;
Pr_water = [12.8 12.3 11.0 10.2 9.41 8.73 7.57 6.83 6.24 5.73 5.10 4.65 4.26 ...
 3.96 3.71 3.41 3.20 3.02 2.81 2.61 2.44 1.76 1.61 1.35 1.16 ...
 1.04 0.99 0.92 0.87 0.84 0.86 0.94 0.99 1.14 1.52 2.70 Inf];
vf_water   = vf_water_1e3 * 1e-3; rho_water  = 1 ./ vf_water;
mu_water   = mu_water_1e6 * 1e-6; nu_water   = mu_water ./ rho_water;
alpha_water= k_water ./ (rho_water .* (cp_water));
switch lower(prop)
    case 'rho',   val = interp1(T_water, rho_water,   T_K, 'pchip');
    case 'mu',    val = interp1(T_water, mu_water,    T_K, 'pchip');
    case 'cp',    val = interp1(T_water, cp_water,    T_K, 'pchip');
    case 'k',     val = interp1(T_water, k_water,     T_K, 'pchip');
    case 'nu',    val = interp1(T_water, nu_water,    T_K, 'pchip');
    case 'alpha', val = interp1(T_water, alpha_water, T_K, 'pchip');
    case 'pr',    val = interp1(T_water, Pr_water,    T_K, 'pchip');
    otherwise,    error('Property not recognized.');
end
end