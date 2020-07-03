
begin

lonstr = 75
lonend = 105
latstr = 25
latend = 40
yrStrt = 198102
yrLast = 200601

casename = (/"Ctrl_FAMIP","Hist_FAMIP"/)
ncase = dimsizes(casename)

season = "DJF"
OMEGA_X = new((/ncase,17,144/), float)
OMEGA_Y = new((/ncase,17,96/), float)
U_X = new((/ncase,17,144/), float)
U_Y = new((/ncase,17,96/), float)
V_X = new((/ncase,17,144/), float)
V_Y = new((/ncase,17,96/), float)
Z_X = new((/ncase,17,144/), float)
Z_Y = new((/ncase,17,96/), float)

do i = 0, ncase-1
    diri = "/home/yangsong3/data-model/wuyt/TPbias/TPbias_"+casename(i)+"/a20191206/"
    filiw = casename(i)+".cam.h0.OMEGA.1979-2005.nc"
    filiu = casename(i)+".cam.h0.U.1979-2005.nc"
    filiv = casename(i)+".cam.h0.V.1979-2005.nc"
    filiz = casename(i)+".cam.h0.Z3.1979-2005.nc"
    
    fu = addfile(diri+filiu, "r")
    fv = addfile(diri+filiv, "r")
    fw = addfile(diri+filiw, "r")
    fz = addfile(diri+filiz, "r")

    TIME = fu->time
    YYYY = cd_calendar(TIME, -1)
    iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

    OMEGAX = fw->OMEGA(iYYYY,:,:,:)
    OMEGA0 = month_to_season(OMEGAX, season)
    OMEGA = dim_avg_n_Wrap(OMEGA0, 0)
    OMEGA = OMEGA*-200.
    OMEGA@units = "10^-2 Pa/s"

    UX = fu->U(iYYYY,:,:,:)
    U0 = month_to_season(UX, season)
    U = dim_avg_n_Wrap(U0, 0)

    VX = fv->V(iYYYY,:,:,:)
    V0 = month_to_season(VX, season)
    V = dim_avg_n_Wrap(V0, 0)

    ZX = fz->Z3(iYYYY,:,:,:)
    Z0 = month_to_season(ZX, season)
    Z = dim_avg_n_Wrap(Z0, 0)
    Z@units = "gpm"

    OMEGA_X(i,:,:) = dim_avg_n_Wrap(OMEGA(:,{latstr:latend},:), 1)
    OMEGA_Y(i,:,:) = dim_avg_n_Wrap(OMEGA(:,:,{lonstr:lonend}), 2)

    U_X(i,:,:) = dim_avg_n_Wrap(U(:,{latstr:latend},:), 1)
    U_Y(i,:,:) = dim_avg_n_Wrap(U(:,:,{lonstr:lonend}), 2)

    V_X(i,:,:) = dim_avg_n_Wrap(V(:,{latstr:latend},:), 1)
    V_Y(i,:,:) = dim_avg_n_Wrap(V(:,:,{lonstr:lonend}), 2)

    Z_X(i,:,:) = dim_avg_n_Wrap(Z(:,{latstr:latend},:), 1)
    Z_Y(i,:,:) = dim_avg_n_Wrap(Z(:,:,{lonstr:lonend}), 2)

    ; OMEGA_X(i,:,:) = OMEGA(:,{34},:)
    ; OMEGA_Y(i,:,:) = OMEGA(:,:,{90})

    ; U_X(i,:,:) = U(:,{34},:)
    ; U_Y(i,:,:) = U(:,:,{90})

    ; V_X(i,:,:) = V(:,{34},:)
    ; V_Y(i,:,:) = V(:,:,{90})

    ; Z_X(i,:,:) = Z(:,{34},:)
    ; Z_Y(i,:,:) = Z(:,:,{90})

end do

printVarSummary(OMEGA_X)

diff_x = OMEGA_X(0,:,:)
diff_x = OMEGA_X(1,:,:) - OMEGA_X(0,:,:)
diff_y = OMEGA_Y(0,:,:)
diff_y = OMEGA_Y(1,:,:) - OMEGA_Y(0,:,:)

diff_ux = U_X(0,:,:)
diff_ux = U_X(1,:,:) - U_X(0,:,:)
diff_uy = U_Y(0,:,:)
diff_uy = U_Y(1,:,:) - U_Y(0,:,:)

diff_vx = V_X(0,:,:)
diff_vx = V_X(1,:,:) - V_X(0,:,:)
diff_vy = V_Y(0,:,:)
diff_vy = V_Y(1,:,:) - V_Y(0,:,:)

diff_zy = Z_Y(0,:,:)
diff_zy = Z_Y(1,:,:) - Z_Y(0,:,:)
diff_zx = Z_X(0,:,:)
diff_zx = Z_X(1,:,:) - Z_X(0,:,:)

;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)

;=========================================================
; plot 
;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/a20200207/"
figname = season+"_zonal_vertical_cell_U_hgt_profile"

; lonstr = 60
; lonend = 120
; latstr = 10
; latend = 50

wks = gsn_open_wks("pdf", figpath+figname)
plot = new(3, graphic)
plotc = plot
plot_topo = plot

;=====================================================
    res_topo                   = True

    res_topo@gsnDraw  = False           ;先别直接画出来
    res_topo@gsnFrame = False           ;也不要翻页

    res_topo@gsnLeftString  = " "              ;不要标示出变量名和单位等信息
    res_topo@gsnRightString = " "

    ;res_topo@gsnMaximize        = True
    ;res_topo@trYMaxF            = 1000.0
    res_topo@gsnYRefLine        = 1000 ;(/100,1000/)          ; create a reference line 
    res_topo@gsnBelowYRefLineColor = "black" ;(/"black","black"/) ; above ref line fill grey
    ;res_topo@gsnXYBelowFillColors = "goldenrod2"
    res_topo@gsnAboveYRefLineColor = "blue"
    res_topo@xyLineColors = "black"
    res_topo@xyLineThicknesses = 4.
;===================================================================

    res = True
    res@gsnDraw  = False
    res@gsnFrame = False
    res@vpHeightF = 0.4
    res@vpWidthF  = 0.7

    res@cnFillOn = True
    res@cnLinesOn = False
    res@cnLineLabelsOn = False
    res@cnInfoLabelOn = False

    res@gsnStringFontHeightF = 0.03
    res@tmXBLabelFontHeightF = 0.02
    res@tmYLLabelFontHeightF = 0.02

    res@pmLabelBarOrthogonalPosF = -0.06
    ;res@lbLabelBarOn = False
    ;res@lbLabelFontHeightF = 0.02
    res@tiYAxisString = ""

    res@tmYLMode      = "Explicit"
    res@tmYLValues    = (/"1000","850","750","500","300","200","100"/)
    res@tmYLLabels    = (/"1000","850","750","500","300","200","100"/)

 ;------------------------------------------
     res@trYMinF   = 100
     res@trYMaxF   = 1000
         
     res@vcRefAnnoOrthogonalPosF =  -0.35 ;-1.0             ; move ref vector up
     ;res@vcGlyphStyle            = "CurlyVector"  ;"LineArrow"    ; turn on curly vectors
     ;res@vcMinDistanceF          = 0.017
     res@vcGlyphStyle            = "LineArrow"
     res@vcMinDistanceF          = 0.03
     res@vcRefLengthF            = 0.04            ; define length of vec ref
     
     res@vcLineArrowThicknessF   = 3.0
     ;res@vcLineArrowHeadMinSizeF = 0.008
     ;res@vcLineArrowHeadMaxSizeF = 0.3
     res@vcVectorDrawOrder = "PostDraw" 
     
     res@vcRefAnnoString2On = False
     res@vcRefAnnoFontHeightF = 0.02

;=====================================================

    cnres = True
    cnres@gsnDraw = False
    cnres@gsnFrame = False

    cnres@cnFillOn = False
    cnres@cnLinesOn = True
    cnres@cnLineLabelsOn = True
    cnres@cnLineLabelFontColor = "red"
    cnres@cnLineColor = "green4"
    cnres@cnLineThicknessF = 5
    cnres@cnInfoLabelOn = False
    ;cnres@cnLineDrawOrder = "PostDraw"
    cnres@cnLabelDrawOrder = "PostDraw"
    ;cnres@cnLevelFlag = "LineAndLabel"
    cnres@cnLineLabelInterval = 1
    cnres@gsnPresHgtHeightLabelOn = False
    cnres@gsnLeftString = ""
    cnres@gsnRightString = ""
;-----------------------------------------------------
cnres@gsnContourPosLineDashPattern = 0
cnres@gsnContourNegLineDashPattern = 2

res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 1000
res@cnMaxLevelValF = 15000
res@cnMinLevelValF = 1000
res@cnFillPalette = "BlueWhiteOrangeRed" ; "MPL_RdBu"
res@vcRefMagnitudeF         = 20.             ; define vector ref mag

; cnres@cnLevelSelectionMode = "ExplicitLevels"
; cnres@cnLevels = (/-10,-7,-4,-2,-0.5,0.5,2,4,7,10/)
cnres@cnLevelSelectionMode = "ManualLevels"
cnres@cnLevelSpacingF = 5
cnres@cnMaxLevelValF = 50
cnres@cnMinLevelValF = -50

res@gsnLeftString = "Ctrl"; (90~S~o~N~E)
plot(0) = gsn_csm_pres_hgt_vector(wks, Z_Y(0,:,{-20:60}), U_Y(0,:,{-20:60}), OMEGA_Y(0,:,{-20:60}), res)
;plot(0) = gsn_csm_pres_hgt(wks, Z_Y(0,:,{-20:60}), res)
plotc(0) = gsn_csm_contour(wks, U_Y(0,:,{-20:60}), cnres)
overlay(plot(0), plotc(0))

res@gsnLeftString = "Sen"; (90~S~o~N~E)
plot(1) = gsn_csm_pres_hgt_vector(wks, Z_Y(1,:,{-20:60}), U_Y(1,:,{-20:60}), OMEGA_Y(1,:,{-20:60}), res)
;plot(0) = gsn_csm_pres_hgt(wks, Z_Y(0,:,{-20:60}), res)
plotc(1) = gsn_csm_contour(wks, U_Y(1,:,{-20:60}), cnres)
overlay(plot(1), plotc(1))

;==========================================================================
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 5
res@cnMaxLevelValF = 50
res@cnMinLevelValF = -50
res@cnFillPalette = "BlueWhiteOrangeRed" ; "MPL_RdBu"
res@vcRefMagnitudeF         = 2.             ; define vector ref mag

; cnres@cnLevelSelectionMode = "ExplicitLevels"
; cnres@cnLevels = (/-10,-7,-4,-2,-0.5,0.5,2,4,7,10/)
cnres@cnLevelSelectionMode = "ManualLevels"
cnres@cnLevelSpacingF = 1
cnres@cnMaxLevelValF = 10
cnres@cnMinLevelValF = -10

res@gsnLeftString = "Diff"; (90~S~o~N~E)
plot(2) = gsn_csm_pres_hgt_vector(wks, diff_zy(:,{-20:60}), diff_uy(:,{-20:60}), diff_y(:,{-20:60}), res)
plotc(2) = gsn_csm_contour(wks, diff_uy(:,{-20:60}), cnres)
overlay(plot(2), plotc(2))


pres = True
pres@gsnMaximize = True
;pres@gsnPanelLabelBar  = True
pres@gsnPanelMainString = "Z (shading) & U (contour)"
gsn_panel(wks, plot, (/3,1/), pres)

end