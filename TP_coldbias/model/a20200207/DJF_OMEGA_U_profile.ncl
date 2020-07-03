
begin

lonstr = 80
lonend = 105
latstr = 30
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


do i = 0, ncase-1
    diri = "/home/yangsong3/data-model/wuyt/TPbias/TPbias_"+casename(i)+"/a20191206/"
    fili = casename(i)+".cam.h0.OMEGA.1979-2005.nc"
    filiu = casename(i)+".cam.h0.U.1979-2005.nc"
    
    f1 = addfile(diri+fili, "r")
    fu = addfile(diri+filiu, "r")

    TIME = f1->time
    YYYY = cd_calendar(TIME, -1)
    iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

    OMEGAX = f1->OMEGA(iYYYY,:,:,:)
    OMEGA0 = month_to_season(OMEGAX, season)
    OMEGA = dim_avg_n_Wrap(OMEGA0, 0)
    OMEGA = OMEGA*100.
    OMEGA@units = "10^-2 Pa/s"

    UX = fu->U(iYYYY,:,:,:)
    U0 = month_to_season(UX, season)
    U = dim_avg_n_Wrap(U0, 0)

    OMEGA_X(i,:,:) = OMEGA(:,{34},:)
    OMEGA_Y(i,:,:) = OMEGA(:,:,{90})

    U_X(i,:,:) = U(:,{34},:)
    U_Y(i,:,:) = U(:,:,{90})
    
    ; OMEGA_X(i,:,:) = dim_avg_n_Wrap(OMEGA(:,{latstr:latend},:), 1)
    ; OMEGA_Y(i,:,:) = dim_avg_n_Wrap(OMEGA(:,:,{lonstr:lonend}), 2)

    ; U_X(i,:,:) = dim_avg_n_Wrap(U(:,{latstr:latend},:), 1)
    ; U_Y(i,:,:) = dim_avg_n_Wrap(U(:,:,{lonstr:lonend}), 2)

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
;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)

;=========================================================
; plot 
;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/a20200207/"
figname = season+"_OMEGA_U_profile"

; lonstr = 60
; lonend = 120
; latstr = 10
; latend = 50

wks = gsn_open_wks("pdf", figpath+figname)
plot = new(4, graphic)
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

    ;res@pmLabelBarOrthogonalPosF = 0.12
    ;res@lbLabelBarOn = False
    ;res@lbLabelFontHeightF = 0.02
    res@cnFillPalette = "BlueWhiteOrangeRed"
    res@tiYAxisString = ""


 ;------------------------------------------
     res@trYMinF   = 100
     res@trYMaxF   = 1000
      
     res@tmYLMode      = "Explicit"
     res@tmYLValues    = (/"1000","850","750","500","300","200","100"/)
     res@tmYLLabels    = (/"1000","850","750","500","300","200","100"/)

     ; res@vcRefAnnoOrthogonalPosF =  -0.35 ;-1.0             ; move ref vector up
     ; ;res@vcGlyphStyle            = "CurlyVector"  ;"LineArrow"    ; turn on curly vectors
     ; ;res@vcMinDistanceF          = 0.017
     ; res@vcGlyphStyle            = "LineArrow"
     ; res@vcMinDistanceF          = 0.03
     ; res@vcRefLengthF            = 0.04            ; define length of vec ref
     
     ; res@vcLineArrowThicknessF   = 3.0
     ; ;res@vcLineArrowHeadMinSizeF = 0.008
     ; ;res@vcLineArrowHeadMaxSizeF = 0.3
     ; res@vcVectorDrawOrder = "PostDraw" 
     
     ; res@vcRefAnnoString2On = False
     ; res@vcRefAnnoFontHeightF = 0.02

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

res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 1
res@cnMaxLevelValF = 10
res@cnMinLevelValF = -10
res@cnFillPalette = "MPL_RdBu"
;res@vcRefMagnitudeF         = 10.             ; define vector ref mag

; cnres@cnLevelSelectionMode = "ExplicitLevels"
; cnres@cnLevels = (/-10,-7,-4,-2,-0.5,0.5,2,4,7,10/)
cnres@cnLevelSelectionMode = "ManualLevels"
cnres@cnLevelSpacingF = 10
cnres@cnMaxLevelValF = 50
cnres@cnMinLevelValF = -50

cnres@gsnContourPosLineDashPattern = 0
cnres@gsnContourNegLineDashPattern = 2

res@gsnLeftString = "Ctrl (34~S~o~N~N)";
plot(0) = gsn_csm_pres_hgt(wks, OMEGA_X(0,:,{60:120}), res)
plotc(0) = gsn_csm_contour(wks, U_X(0,:,{60:120}), cnres)
overlay(plot(0), plotc(0))

res@gsnLeftString = "Ctrl (90~S~o~N~E)";
plot(1) = gsn_csm_pres_hgt(wks, OMEGA_Y(0,:,{-20:60}), res)
plotc(1) = gsn_csm_contour(wks, U_Y(0,:,{-20:60}), cnres)
overlay(plot(1), plotc(1))

;-----------------------------------------------------
delete(res@cnLevelSpacingF)
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 0.2
res@cnMaxLevelValF = 1
res@cnMinLevelValF = -1
res@cnFillPalette = "MPL_RdBu"
;res@vcRefMagnitudeF         = 10.             ; define vector ref mag

; cnres@cnLevelSelectionMode = "ExplicitLevels"
; cnres@cnLevels = (/-10,-7,-4,-2,-0.5,0.5,2,4,7,10/)
cnres@cnLevelSelectionMode = "ManualLevels"
cnres@cnLevelSpacingF = 3
cnres@cnMaxLevelValF = 15
cnres@cnMinLevelValF = -15

cnres@gsnContourPosLineDashPattern = 0
cnres@gsnContourNegLineDashPattern = 2

res@gsnLeftString = "Diff (34~S~o~N~N)";
plot(2) = gsn_csm_pres_hgt(wks, diff_x(:,{60:120}), res)
plotc(2) = gsn_csm_contour(wks, diff_ux(:,{60:120}), cnres)
overlay(plot(2), plotc(2))

res@gsnLeftString = "Diff (90~S~o~N~E)";
plot(3) = gsn_csm_pres_hgt(wks, diff_y(:,{-20:60}), res)
plotc(3) = gsn_csm_contour(wks, diff_uy(:,{-20:60}), cnres)
overlay(plot(3), plotc(3))


pres = True
pres@gsnMaximize = True
;pres@gsnPanelLabelBar  = True
pres@gsnPanelMainString = "Omega (shading) & U (contour)"
gsn_panel(wks, plot, (/2,2/), pres)

end