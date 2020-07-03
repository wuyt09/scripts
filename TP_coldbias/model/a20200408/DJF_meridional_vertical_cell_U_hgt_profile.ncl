
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
nyear = 25
U_X = new((/ncase,nyear,17,144/), float)
U_Y = new((/ncase,nyear,17,96/), float)
V_X = U_X
V_Y = U_Y
W_X = U_X
W_Y = U_Y
Z_X = U_X
Z_Y = U_Y

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
    OMEGA = month_to_season(OMEGAX, season)
    OMEGA = OMEGA*100.
    OMEGA@units = "10^-2 Pa/s"    
    W_X(i,:,:,:) = dim_avg_n_Wrap(OMEGA(:,:,{latstr:latend},:), 2)
    W_Y(i,:,:,:) = dim_avg_n_Wrap(OMEGA(:,:,:,{lonstr:lonend}), 3)

    UX = fu->U(iYYYY,:,:,:)
    U = month_to_season(UX, season)
    U_X(i,:,:,:) = dim_avg_n_Wrap(U(:,:,{latstr:latend},:), 2)
    U_Y(i,:,:,:) = dim_avg_n_Wrap(U(:,:,:,{lonstr:lonend}), 3)

    VX = fv->V(iYYYY,:,:,:)
    V = month_to_season(VX, season)
    V_X(i,:,:,:) = dim_avg_n_Wrap(V(:,:,{latstr:latend},:), 2)
    V_Y(i,:,:,:) = dim_avg_n_Wrap(V(:,:,:,{lonstr:lonend}), 3)

    ZX = fz->Z3(iYYYY,:,:,:)
    Z = month_to_season(ZX, season)
    Z@units = "gpm"    
    Z_X(i,:,:,:) = dim_avg_n_Wrap(Z(:,:,{latstr:latend},:), 2)
    Z_Y(i,:,:,:) = dim_avg_n_Wrap(Z(:,:,:,{lonstr:lonend}), 3)
end do
printVarSummary(U_X)

;=======================================================================
; ttest
;=======================================================================

ux_ave = dim_avg_n_Wrap(U_X, 1)
ux_var = dim_variance_n_Wrap(U_X, 1)
probux = ttest(ux_ave(0,:,:), ux_var(0,:,:), nyear, ux_ave(1,:,:), ux_var(1,:,:), nyear, False, False)
copy_VarCoords(ux_ave(0,:,:), probux)
printVarSummary(probux)
uy_ave = dim_avg_n_Wrap(U_Y, 1)
uy_var = dim_variance_n_Wrap(U_Y, 1)
probuy = ttest(uy_ave(0,:,:), uy_var(0,:,:), nyear, uy_ave(1,:,:), uy_var(1,:,:), nyear, False, False)
copy_VarCoords(uy_ave(0,:,:), probuy)
printVarSummary(probuy)


vx_ave = dim_avg_n_Wrap(V_X, 1)
vx_var = dim_variance_n_Wrap(V_X, 1)
probvx = ttest(vx_ave(0,:,:), vx_var(0,:,:), nyear, vx_ave(1,:,:), vx_var(1,:,:), nyear, False, False)
copy_VarCoords(vx_ave(0,:,:), probvx)
printVarSummary(probvx)
vy_ave = dim_avg_n_Wrap(V_Y, 1)
vy_var = dim_variance_n_Wrap(V_Y, 1)
probvy = ttest(vy_ave(0,:,:), vy_var(0,:,:), nyear, vy_ave(1,:,:), vy_var(1,:,:), nyear, False, False)
copy_VarCoords(vy_ave(0,:,:), probvy)
printVarSummary(probvy)

wx_ave = dim_avg_n_Wrap(W_X, 1)
wx_var = dim_variance_n_Wrap(W_X, 1)
probwx = ttest(wx_ave(0,:,:), wx_var(0,:,:), nyear, wx_ave(1,:,:), wx_var(1,:,:), nyear, False, False)
copy_VarCoords(wx_ave(0,:,:), probwx)
printVarSummary(probwx)
wy_ave = dim_avg_n_Wrap(W_Y, 1)
wy_var = dim_variance_n_Wrap(W_Y, 1)
probwy = ttest(wy_ave(0,:,:), wy_var(0,:,:), nyear, wy_ave(1,:,:), wy_var(1,:,:), nyear, False, False)
copy_VarCoords(wy_ave(0,:,:), probwy)
printVarSummary(probwy)

zx_ave = dim_avg_n_Wrap(Z_X, 1)
zx_var = dim_variance_n_Wrap(Z_X, 1)
probzx = ttest(zx_ave(0,:,:), zx_var(0,:,:), nyear, zx_ave(1,:,:), zx_var(1,:,:), nyear, False, False)
copy_VarCoords(zx_ave(0,:,:), probzx)
printVarSummary(probzx)
zy_ave = dim_avg_n_Wrap(Z_Y, 1)
zy_var = dim_variance_n_Wrap(Z_Y, 1)
probzy = ttest(zy_ave(0,:,:), zy_var(0,:,:), nyear, zy_ave(1,:,:), zy_var(1,:,:), nyear, False, False)
copy_VarCoords(zy_ave(0,:,:), probzy)
printVarSummary(probzy)


diff_wx = wx_ave(0,:,:)
diff_wx = wx_ave(1,:,:) - wx_ave(0,:,:)
diff_wy = wy_ave(0,:,:)
diff_wy = wy_ave(1,:,:) - wy_ave(0,:,:)

diff_ux = ux_ave(0,:,:)
diff_ux = ux_ave(1,:,:) - ux_ave(0,:,:)
diff_uy = uy_ave(0,:,:)
diff_uy = uy_ave(1,:,:) - uy_ave(0,:,:)

diff_vx = vx_ave(0,:,:)
diff_vx = vx_ave(1,:,:) - vx_ave(0,:,:)
diff_vy = vy_ave(0,:,:)
diff_vy = vy_ave(1,:,:) - vy_ave(0,:,:)

diff_zy = zy_ave(0,:,:)
diff_zy = zy_ave(1,:,:) - zy_ave(0,:,:)
diff_zx = zx_ave(0,:,:)
diff_zx = zx_ave(1,:,:) - zx_ave(0,:,:)

;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)

;=========================================================
; plot 
;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/a20200408/"
figname = season+"_meridional_vertical_cell_U_hgt_profile"

; lonstr = 60
; lonend = 120
; latstr = 10
; latend = 50

wks = gsn_open_wks("pdf", figpath+figname)
plot = new(3, graphic)
plotc = plot
plot_topo = plot
plot_prob = plot

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
     res@gsnPresHgtHeightLabelOn = False
     res@tmYROn = False

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
res@vcRefMagnitudeF         = 10.             ; define vector ref mag

; cnres@cnLevelSelectionMode = "ExplicitLevels"
; cnres@cnLevels = (/-10,-7,-4,-2,-0.5,0.5,2,4,7,10/)
cnres@cnLevelSelectionMode = "ManualLevels"
cnres@cnLevelSpacingF = 10
cnres@cnMaxLevelValF = 50
cnres@cnMinLevelValF = -50

res@gsnLeftString = "Ctrl"; (90~S~o~N~E)
plot(0) = gsn_csm_pres_hgt_vector(wks, zx_ave(0,:,{60:120}), ux_ave(0,:,{60:120}), -wx_ave(0,:,{60:120}), res)
plotc(0) = gsn_csm_contour(wks, ux_ave(0,:,{60:120}), cnres)
overlay(plot(0), plotc(0))

res@gsnLeftString = "Sen"; (90~S~o~N~E)
plot(1) = gsn_csm_pres_hgt_vector(wks, zx_ave(1,:,{60:120}), ux_ave(1,:,{60:120}), -wx_ave(1,:,{60:120}), res)
plotc(1) = gsn_csm_contour(wks, ux_ave(1,:,{60:120}), cnres)
overlay(plot(1), plotc(1))

;==========================================================================
;----------------------------------------------
; draw prob
;----------------------------------------------
 res3  = True
 res3@gsnDraw   = False
 res3@gsnFrame  = False

 ;res3@cnFillOn  = False 
 res3@cnLinesOn = False
 res3@cnLineLabelsOn = False
 res3@cnInfoLabelOn  = False
 ;res3@lbLabelBarOn   = False

 res3@gsnRightString = ""

 res3@cnLevelSelectionMode = "ManualLevels"
 res3@cnMinLevelValF = 0.00
 res3@cnMaxLevelValF = 0.10
 res3@cnLevelSpacingF = 0.01

 ; res3@cnFillPattern = 17 ;(/17/)
 res3@cnFillColor  =  "purple"
 res3@gsnLeftString = ""
 res3@cnFillDotSizeF = 0.003
 res3@cnFillScaleF  = 1.2
 ;----------------------------------------------
 ;----------------------------------------------

res@cnFillPalette = "BlueWhiteOrangeRed" ; "MPL_RdBu"
res@vcRefMagnitudeF         = 2.             ; define vector ref mag

cnres@cnLevelSelectionMode = "ManualLevels"
cnres@cnLevelSpacingF = 2
cnres@cnMaxLevelValF = 10
cnres@cnMinLevelValF = -10

res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 10
res@cnMaxLevelValF = 50
res@cnMinLevelValF = -50

res@gsnLeftString = "Diff"; (90~S~o~N~E)
diff_zx = mask(diff_zx, probzx.le.0.055, True)
; diff_wy = mask(diff_wy, probwy.le.0.055, True)
; diff_uy = mask(diff_uy, probuy.le.0.055, True)
plot(2) = gsn_csm_pres_hgt_vector(wks, diff_zx(:,{60:120}), diff_ux(:,{60:120}), -diff_wx(:,{60:120}), res)
plotc(2) = gsn_csm_contour(wks, diff_ux(:,{60:120}), cnres)
overlay(plot(2), plotc(2))

; res@vcLineArrowColor = "red"
; diff_zy = 0.0001
; diff_wy = mask(diff_wy, probwy.le.0.055, True)
; diff_uy = mask(diff_uy, probuy.le.0.055, True)
; plot_prob(0) = gsn_csm_pres_hgt_vector(wks, diff_zy(:,{-20:60}), diff_uy(:,{-20:60}), -diff_wy(:,{-20:60}), res)
; overlay(plot(2), plot_prob(0))

plot_prob(1) = gsn_csm_contour(wks,probzx,res3)
plot_prob(1) = ShadeLtContour(plot_prob(1),0.055,17)
overlay(plot(2), plot_prob(1))

; diff_uy = mask(diff_uy, probuy.le.0.055, True)
; cnres@cnLineColor = "black"
; plot_prob(2) = gsn_csm_contour(wks,diff_uy(:,{-20:60}),cnres)
; plot_prob(2) = ShadeLtContour(plot_prob(2),0.055,17)
; overlay(plot(2), plot_prob(2))

pres = True
pres@gsnMaximize = True
;pres@gsnPanelLabelBar  = True
pres@gsnPanelMainString = season+" mean: Z (shading) & U (contour), 25~S~o~N~N - 40~S~o~N~N"
gsn_panel(wks, plot, (/1,3/), pres)

end