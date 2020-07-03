
begin

lonstr = 75
lonend = 105
latstr = 25
latend = 40
yrStrt = 198102
yrLast = 200601

casename = (/"Ctrl_FAMIP","Hist_FAMIP"/)
ncase = dimsizes(casename)
nyear = 25

season = "DJF"
Q_X = new((/ncase,nyear,17,144/), float)
Q_Y = new((/ncase,nyear,17,96/), float)
T_X = Q_X
T_Y = Q_Y

do i = 0, ncase-1
    diri = "/home/yangsong3/data-model/wuyt/TPbias/TPbias_"+casename(i)+"/a20191206/"
    filiq = casename(i)+".cam.h0.Q.1979-2005.nc"
    filit = casename(i)+".cam.h0.T.1979-2005.nc"
    
    fq = addfile(diri+filiq, "r")
    ft = addfile(diri+filit, "r")

    TIME = fq->time
    YYYY = cd_calendar(TIME, -1)
    iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

    QX = fq->Q(iYYYY,:,:,:)
    Q = month_to_season(QX, season)
    TX = ft->T(iYYYY,:,:,:)
    T = month_to_season(TX, season)

    Q_X(i,:,:,:) = dim_avg_n_Wrap(Q(:,:,{latstr:latend},:), 2)
    Q_Y(i,:,:,:) = dim_avg_n_Wrap(Q(:,:,:,{lonstr:lonend}), 3)

    T_X(i,:,:,:) = dim_avg_n_Wrap(T(:,:,{latstr:latend},:), 2)
    T_Y(i,:,:,:) = dim_avg_n_Wrap(T(:,:,:,{lonstr:lonend}), 3)

end do
printVarSummary(Q_X)

;=======================================================================
; ttest
;=======================================================================

qx_ave = dim_avg_n_Wrap(Q_X, 1)
qx_var = dim_variance_n_Wrap(Q_X, 1)
probqx = ttest(qx_ave(0,:,:), qx_var(0,:,:), nyear, qx_ave(1,:,:), qx_var(1,:,:), nyear, False, False)
copy_VarCoords(qx_ave(0,:,:), probqx)
printVarSummary(probqx)

qy_ave = dim_avg_n_Wrap(Q_Y, 1)
qy_var = dim_variance_n_Wrap(Q_Y, 1)
probqy = ttest(qy_ave(0,:,:), qy_var(0,:,:), nyear, qy_ave(1,:,:), qy_var(1,:,:), nyear, False, False)
copy_VarCoords(qy_ave(0,:,:), probqy)
printVarSummary(probqy)



tx_ave = dim_avg_n_Wrap(T_X, 1)
tx_var = dim_variance_n_Wrap(T_X, 1)
probtx = ttest(tx_ave(0,:,:), tx_var(0,:,:), nyear, tx_ave(1,:,:), tx_var(1,:,:), nyear, False, False)
copy_VarCoords(tx_ave(0,:,:), probtx)
printVarSummary(probtx)

ty_ave = dim_avg_n_Wrap(T_Y, 1)
ty_var = dim_variance_n_Wrap(T_Y, 1)
probty = ttest(ty_ave(0,:,:), ty_var(0,:,:), nyear, ty_ave(1,:,:), ty_var(1,:,:), nyear, False, False)
copy_VarCoords(ty_ave(0,:,:), probty)
printVarSummary(probty)


diffqx = qx_ave(0,:,:)
diffqx = qx_ave(1,:,:) - qx_ave(0,:,:)
diffqy = qy_ave(0,:,:)
diffqy = qy_ave(1,:,:) - qy_ave(0,:,:)


difftx = tx_ave(0,:,:)
difftx = tx_ave(1,:,:) - tx_ave(0,:,:)
diffty = ty_ave(0,:,:)
diffty = ty_ave(1,:,:) - ty_ave(0,:,:)
printVarSummary(diffqx)

;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)

;=========================================================
; plot 
;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/a20200408/"
figname = season+"_wv_airT_profile"

; lonstr = 60
; lonend = 120
; latstr = 10
; latend = 50

wks = gsn_open_wks("pdf", figpath+figname)
plot = new(4, graphic)
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
    res@tmXBLabelFontHeightF = 0.025
    res@tmYLLabelFontHeightF = 0.025

    ;res@pmLabelBarOrthogonalPosF = 0.12
    ;res@lbLabelBarOn = False
    ;res@lbLabelFontHeightF = 0.02
    res@cnFillPalette = "BlueWhiteOrangeRed"
    res@tiYAxisString = ""
    res@gsnPresHgtHeightLabelOn = False

 ;------------------------------------------
     res@trYMinF   = 100
     res@trYMaxF   = 1000
      
     res@tmYLMode      = "Explicit"
     res@tmYLValues    = (/"1000","850","750","500","300","200","100"/)
     res@tmYLLabels    = (/"1000","850","750","500","300","200","100"/)

;=====================================================

    cnres = True
    cnres@gsnDraw = False
    cnres@gsnFrame = False

    cnres@cnFillOn = False
    cnres@cnLinesOn = True
    cnres@cnLineLabelsOn = True
    cnres@cnLineLabelFontColor = "red"
    cnres@cnLineColor = "grey50"
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
;-----------------------------------------------------

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
 res3@cnFillColor  =  "black"
 res3@gsnLeftString = ""
 res3@cnFillDotSizeF = 0.003
 res3@cnFillScaleF  = 1.2
 ;----------------------------------------------
 ;----------------------------------------------
res@gsnLeftString = "Zonal mean ("+latstr+"-"+latend+"~S~o~N~N)"; (34~S~o~N~N)
res@gsnRightString = "Ctrl"
res@cnFillPalette = "CBR_wet"
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 0.001
res@cnMaxLevelValF = 0.01
res@cnMinLevelValF = 0.

; cnres@cnLevelSelectionMode = "ExplicitLevels"
; cnres@cnLevels = (/-10,-7,-4,-2,-0.5,0.5,2,4,7,10/)
cnres@cnLevelSelectionMode = "ManualLevels"
cnres@cnLevelSpacingF = 10.
cnres@cnMaxLevelValF = 300.
cnres@cnMinLevelValF = 200.
plot(0) = gsn_csm_pres_hgt(wks, qx_ave(0,:,{60:120}), res)
plotc(0) = gsn_csm_contour(wks, tx_ave(0,:,{60:120}), cnres)
overlay(plot(0), plotc(0))


res@gsnRightString = "diff"
res@cnFillPalette = "CBR_drywet"
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 0.0001
res@cnMaxLevelValF = 0.001
res@cnMinLevelValF = -0.001

; cnres@cnLevelSelectionMode = "ExplicitLevels"
; cnres@cnLevels = (/-10,-7,-4,-2,-0.5,0.5,2,4,7,10/)
cnres@cnLevelSelectionMode = "ManualLevels"
cnres@cnLevelSpacingF = 0.5
cnres@cnMaxLevelValF = 5.
cnres@cnMinLevelValF = -5.
plot(1) = gsn_csm_pres_hgt(wks, diffqx(:,{60:120}), res)
plotc(1) = gsn_csm_contour(wks, difftx(:,{60:120}), cnres)
overlay(plot(1), plotc(1))



res@gsnLeftString = "Meridional mean ("+lonstr+"-"+lonend+"~S~o~N~E)";
res@gsnRightString = "Ctrl"
res@cnFillPalette = "CBR_wet"
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 0.001
res@cnMaxLevelValF = 0.01
res@cnMinLevelValF = 0.

; cnres@cnLevelSelectionMode = "ExplicitLevels"
; cnres@cnLevels = (/-10,-7,-4,-2,-0.5,0.5,2,4,7,10/)
cnres@cnLevelSelectionMode = "ManualLevels"
cnres@cnLevelSpacingF = 10.
cnres@cnMaxLevelValF = 300.
cnres@cnMinLevelValF = 200.
plot(2) = gsn_csm_pres_hgt(wks, qy_ave(0,:,{-20:60}), res)
plotc(2) = gsn_csm_contour(wks, ty_ave(0,:,{-20:60}), cnres)
overlay(plot(2), plotc(2))

res@gsnRightString = "diff"
res@cnFillPalette = "CBR_drywet"
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 0.0001
res@cnMaxLevelValF = 0.001
res@cnMinLevelValF = -0.001

cnres@cnLevelSelectionMode = "ManualLevels"
cnres@cnLevelSpacingF = 0.5
cnres@cnMaxLevelValF = 5.
cnres@cnMinLevelValF = -5.
plot(3) = gsn_csm_pres_hgt(wks, diffqy(:,{-20:60}), res)
plotc(3) = gsn_csm_contour(wks, diffty(:,{-20:60}), cnres)
overlay(plot(3), plotc(3))

 


plot_prob(0) = gsn_csm_contour(wks,probqx,res3)
plot_prob(0) = ShadeLtContour(plot_prob(0),0.055,17)
overlay(plot(1), plot_prob(0))
plot_prob(1) = gsn_csm_contour(wks,probqy,res3)
plot_prob(1) = ShadeLtContour(plot_prob(1),0.055,17)
overlay(plot(3), plot_prob(1))

difftx = mask(difftx, probtx.le.0.055, True)
diffty = mask(diffty, probty.le.0.055, True)
cnres@cnLineColor = "red"
plot_prob(2) = gsn_csm_contour(wks, difftx(:,{60:120}), cnres)
overlay(plot(1), plot_prob(2))
plot_prob(3) = gsn_csm_contour(wks, diffty(:,{-20:60}), cnres)
overlay(plot(3), plot_prob(3))

pres = True
pres@gsnMaximize = True
pres@gsnPanelMainString = "Q (Kg/Kg, shading) & T (K, contour)"
gsn_panel(wks, plot, (/2,2/), pres)

end