
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
qu = new((/ncase,nyear,96,144/), float)
qv = qu
qq = qu
dv = qu

LEV = 500

do i = 0, ncase-1
    diri = "/home/yangsong3/data-model/wuyt/TPbias/TPbias_"+casename(i)+"/a20191206/"
    fili1 = casename(i)+".cam.h0.U.1979-2005.nc"
    fili2 = casename(i)+".cam.h0.V.1979-2005.nc"
    fili3 = casename(i)+".cam.h0.Q.1979-2005.nc"
    fili4 = casename(i)+".cam.h0.PS.1979-2005.nc"

    f1 = addfile(diri+fili1, "r")
    f2 = addfile(diri+fili2, "r")
    f3 = addfile(diri+fili3, "r")
    f4 = addfile(diri+fili4, "r")

    TIME = f1->time
    YYYY = cd_calendar(TIME, -1)
    iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

    U = f1->U(iYYYY,:,:,:)
    V = f2->V(iYYYY,:,:,:)
    Q = f3->Q(iYYYY,:,:,:)
    PS = f4->PS(iYYYY,:,:)

    lev = Q&level
    plev = lev
    plev = plev *100
    plev@units = "Pa"

    ptop = min(plev)
    psfc = PS
    dP = dpres_plevel_Wrap(plev, psfc, ptop, 0)

    QUdP = Q
    QUdP = Q*U*dP
    QVdP = Q
    QVdP = Q*V*dP
    QdP = Q
    QdP = Q*dP

    QUX = dim_sum_n_Wrap(QUdP(:,{1000:200},:,:), 1)
    QVX = dim_sum_n_Wrap(QVdP(:,{1000:200},:,:), 1)
    QX = dim_sum_n_Wrap(QdP(:,{1000:200},:,:), 1)
    ; QUX = dim_sum_n_Wrap(QUdP, 1)
    ; QVX = dim_sum_n_Wrap(QVdP, 1)
    QUX = QUX/9.8
    QVX = QVX/9.8
    QX = QX/9.8
    lat = Q&lat
    lon = Q&lon
    DVX = uv2dv_cfd(QUX, QVX, lat, lon, 1)
    copy_VarCoords(QUX, DVX)

    qu(i,:,:,:) = month_to_season(QUX, season)
    qv(i,:,:,:) = month_to_season(QVX, season)
    dv(i,:,:,:) = month_to_season(DVX, season)    
    qq(i,:,:,:) = month_to_season(QX, season)

end do

dv = dv*1000.
dv@units = "g/s/m"
copy_VarCoords(qu, dv)
printVarSummary(dv)

;=======================================================================
; ttest
;=======================================================================

qu_ave = dim_avg_n_Wrap(qu, 1)
qu_var = dim_variance_n_Wrap(qu, 1)
prob_qu = ttest(qu_ave(0,:,:), qu_var(0,:,:), nyear, qu_ave(1,:,:), qu_var(1,:,:), nyear, False, False)
copy_VarCoords(qu_ave(0,:,:), prob_qu)
diff_qu = qu_ave(0,:,:)
diff_qu = qu_ave(1,:,:) - qu_ave(0,:,:)

qv_ave = dim_avg_n_Wrap(qv, 1)
qv_var = dim_variance_n_Wrap(qv, 1)
prob_qv = ttest(qv_ave(0,:,:), qv_var(0,:,:), nyear, qv_ave(1,:,:), qv_var(1,:,:), nyear, False, False)
copy_VarCoords(qv_ave(0,:,:), prob_qv)
diff_qv = qv_ave(0,:,:)
diff_qv = qv_ave(1,:,:) - qv_ave(0,:,:)

qq_ave = dim_avg_n_Wrap(qq, 1)
qq_var = dim_variance_n_Wrap(qq, 1)
prob_qq = ttest(qq_ave(0,:,:), qq_var(0,:,:), nyear, qq_ave(1,:,:), qq_var(1,:,:), nyear, False, False)
copy_VarCoords(qq_ave(0,:,:), prob_qq)
diff_qq = qq_ave(0,:,:)
diff_qq = qq_ave(1,:,:) - qq_ave(0,:,:)

dv_ave = dim_avg_n_Wrap(dv, 1)
dv_var = dim_variance_n_Wrap(dv, 1)
prob_dv = ttest(dv_ave(0,:,:), dv_var(0,:,:), nyear, dv_ave(1,:,:), dv_var(1,:,:), nyear, False, False)
copy_VarCoords(dv_ave(0,:,:), prob_dv)
diff_dv = dv_ave(0,:,:)
diff_dv = dv_ave(1,:,:) - dv_ave(0,:,:)


;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)

;=========================================================
; plot 
;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/a20200408/"
figname = season+"_wv_flux_whole_column"

wks = gsn_open_wks("pdf", figpath+figname)
plot = new(2, graphic)
plotv = plot
plot_topo = plot
plot_prob = plot

 res_topo            = True
 res_topo@cnFillOn   = False
 res_topo@cnLinesOn  = True
 res_topo@gsnFrame   = False
 res_topo@gsnDraw    = False
 res_topo@cnLineLabelsOn         = False
 res_topo@cnLineThicknessF       = 3
 res_topo@cnLineColor            = "grey"
 res_topo@cnInfoLabelOn          = False
 res_topo@cnLevelSelectionMode   ="ExplicitLevels"
 res_topo@cnLevels               = (/1500/)
 res_topo@gsnLeftString = ""
 res_topo@gsnRightString = ""

 ;----------------------------------------------------------

  res   = True
  res@gsnDraw  = False
  res@gsnFrame = False
  res@gsnAddCyclic = False

  res@cnFillOn      = True
  res@cnLinesOn     = False
  res@cnInfoLabelOn = False
  res@cnLineLabelsOn = False

  res@mpOutlineOn   = True
  res@mpFillOn      = False
  res@mpMinLatF     = latstr
  res@mpMaxLatF     = latend
  res@mpMinLonF     = lonstr
  res@mpMaxLonF     = lonend
  res@mpGeophysicalLineThicknessF = 4
  res@gsnRightStringFontHeightF = 0.035

  res@gsnStringFontHeightF = 0.04

  res@tmXBLabelFontHeightF = 0.025
  res@tmYLLabelFontHeightF = 0.025

  res@lbLabelFontHeightF   = 0.015
  res@pmLabelBarHeightF    = 0.1
  res@pmLabelBarWidthF     = 0.75
  res@pmLabelBarOrthogonalPosF = 0.12
  ;res@o=pmLabelBarParallelPosF = 

  res@gsnSpreadColors = True
  res@gsnSpreadColorEnd   = 0
  res@gsnSpreadColorStart = 127
  res@cnFillPalette = "MPL_BrBG" ;CBR_drywet
  res@cnLevelSelectionMode = "ManualLevels"
  res@gsnRightString = "g/m~S~2~N~/s"
  res@cnMaxLevelValF  = 0.05
  res@cnMinLevelValF  = -0.05
  res@cnLevelSpacingF = 0.005
  ; res@cnLevelSelectionMode = "ExplicitLevels"
  ; res@cnLevels = (/-0.01,-0.007,-0.005,-0.003,-0.001,0.001,0.003,0.005,0.007,0.01/)
  ; res@cnFillColors = (/10,18,26,34,42,65,82,90,98,106,114/)
  
  ;==========================================================
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
 res3@cnFillColor  =  "red"
 res3@gsnLeftString = ""
 res3@cnFillDotSizeF = 0.003
 res3@cnFillScaleF  = 1.2
 ;----------------------------------------------
 ;----------------------------------------------

  res@gsnLeftString = "Ctrl"
  plot(0) = gsn_csm_contour_map(wks, dv_ave(0,:,:), res)
  plot_topo(0) = gsn_csm_contour(wks, topo, res_topo)
  overlay(plot(0), plot_topo(0))
  
  res@gsnLeftString = "diff"
  res@cnFillPalette = "MPL_BrBG"
  res@cnLevelSelectionMode = "ExplicitLevels"
  res@cnLevels = (/-0.01,-0.007,-0.005,-0.003,-0.001,0.001,0.003,0.005,0.007,0.01/)
  ;res@cnFillColors = (/10,18,26,34,42,65,82,90,98,106,114/)
  res@cnFillColors = (/114,106,98,90,82,65,42,34,26,18,10/)
  plot(1) = gsn_csm_contour_map(wks, diff_dv, res)
  plot_topo(1) = gsn_csm_contour(wks, topo, res_topo)
  overlay(plot(1), plot_topo(1))
  plot_prob(0) = gsn_csm_contour(wks,prob_dv,res3)
  plot_prob(0) = ShadeLtContour(plot_prob(0),0.055,17)
  overlay(plot(1), plot_prob(0))

;-----------------------------------------------------------
;===========================================================
 vecres  = True
 vecres@gsnDraw  = False
 vecres@gsnFrame = False

 vecres@vcGlyphStyle     = "LineArrow"   ; curly vectors"Fillarrow" 
 vecres@vcRefLengthF     = 0.07            ; define length of vec ref
 vecres@vcMinDistanceF   = 0.02
 vecres@vcLineArrowThicknessF = 3.0
 vecres@tiXAxisString    = " "             ; turn off axis label
 vecres@vcRefAnnoOn      = True            ; draw the reference vector annotation.
 vecres@vcRefAnnoString2On =  False          ;display a string below or to the right
 vecres@vcRefAnnoOrientation = "horizontal"
 vecres@vcRefAnnoString2  = " "
 vecres@vcRefAnnoOrthogonalPosF = -1.0  ; move ref vector into plot
 vecres@vcFillArrowsOn       = True
 vecres@vcLineArrowHeadMinSizeF = 0.008
 vecres@vcLineArrowHeadMaxSizeF = 0.3
 vecres@gsnLeftString  = ""
 ;vecres@gsnRightString  = "g*s/Kg"
 vecres@gsnRightString  = ""
 vecres@gsnRightStringFontHeightF = 0.04

 ;vecres@vcLineArrowColor = "gray52"
 vecres@vcRefAnnoFontHeightF = 0.01
 vecres@vcRefMagnitudeF  = 50             ; define vector ref mag
 vecres@vcRefAnnoString1  = "50 Kg/m/s"
 plotv(0) = gsn_csm_vector(wks, qu_ave(0,:,:), qv_ave(0,:,:), vecres)
 overlay(plot(0), plotv(0))

 vecres@vcRefLengthF     = 0.07            ; define length of vec ref
 vecres@vcRefMagnitudeF  = 10             ; define vector ref mag
 vecres@vcRefAnnoString1  = "10 Kg/m/s"
 plotv(1) = gsn_csm_vector(wks, diff_qu, diff_qv, vecres)
 overlay(plot(1), plotv(1))
 
 diff_qu = mask(diff_qu, prob_qu.le.0.055 .and. prob_qv.le.0.055, True)
 diff_qv = mask(diff_qv, prob_qu.le.0.055 .and. prob_qv.le.0.055, True)
 vecres@vcLineArrowColor = "red"
 plot_prob(1) = gsn_csm_vector(wks, diff_qu, diff_qv, vecres)
 overlay(plot(1), plot_prob(1))


 pres = True
 pres@gsnMaximize = True
 pres@gsnPanelMainString = season+" mean: Wv flux and its divergence (whole column)"

 gsn_panel(wks, plot, (/2,1/), pres)

end
