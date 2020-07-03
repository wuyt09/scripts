
begin

lonstr = 75
lonend = 105
latstr = 25
latend = 40
yrStrt = 198102
yrLast = 200601

casename = (/"Ctrl_FAMIP","Hist_FAMIP"/)
ncase = dimsizes(casename)
LEV = (/600,500,400,300/)
nlev = dimsizes(LEV)

season = "DJF"
qu = new((/ncase,nlev,96,144/), float)
qv = qu
dv = qu


do i = 0, ncase-1
    diri = "/home/yangsong3/data-model/wuyt/TPbias/TPbias_"+casename(i)+"/a20191206/"
    fili1 = casename(i)+".cam.h0.U.1979-2005.nc"
    fili2 = casename(i)+".cam.h0.V.1979-2005.nc"
    fili3 = casename(i)+".cam.h0.Q.1979-2005.nc"

    f1 = addfile(diri+fili1, "r")
    f2 = addfile(diri+fili2, "r")
    f3 = addfile(diri+fili3, "r")

    TIME = f1->time
    YYYY = cd_calendar(TIME, -1)
    iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

    U = f1->U(iYYYY,{LEV},:,:)
    V = f2->V(iYYYY,{LEV},:,:)
    Q = f3->Q(iYYYY,{LEV},:,:)

    QU = Q
    QU = Q * U
    QV = Q
    QV = Q * V

    lat = Q&lat
    lon = Q&lon
    DV = uv2dv_cfd(QU, QV, lat, lon, 1)
    copy_VarCoords(QU, DV)

    QUX = month_to_season(QU, season)
    QVX = month_to_season(QV, season)
    DVX = month_to_season(DV, season)

    qu(i,:,:,:) = dim_avg_n_Wrap(QUX, 0)
    qv(i,:,:,:) = dim_avg_n_Wrap(QVX, 0)
    dv(i,:,:,:) = dim_avg_n_Wrap(DVX, 0)

end do
qu = qu*1e+3
qv = qv*1e+3
dv = dv*1e+6
dv@units = "10^-6 s^-1"
copy_VarCoords(qu, dv)
printVarSummary(dv)

; printMinMax(qu, 0)
; printMinMax(qv, 0)
; printMinMax(qq, 0)
; exit

qu_diff = qu(0,:,:,:)
qu_diff = qu(1,:,:,:) - qu(0,:,:,:)
qv_diff = qv(0,:,:,:)
qv_diff = qv(1,:,:,:) - qv(0,:,:,:)
dv_diff = dv(0,:,:,:)
dv_diff = dv(1,:,:,:) - dv(0,:,:,:)

; printMinMax(qu_diff, 0)
; printMinMax(qv_diff, 0)
; printMinMax(qq_diff, 0)
; exit
;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)

;=========================================================
; plot 
;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/a20200207/"
figname = season+"_wv_flux_multilevel"

wks = gsn_open_wks("pdf", figpath+figname)
plot = new(8, graphic)
plotv = plot
plot_topo = plot

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


  dv_diff = dv_diff*1e+3

  do i = 0,nlev-1
  res@cnFillPalette = "MPL_BrBG"
  res@cnLevelSelectionMode = "ExplicitLevels"
  res@gsnRightString = "Ctrl"
  res@cnMaxLevelValF  = 0.05
  res@cnMinLevelValF  = -0.05
  res@cnLevelSpacingF = 0.005
  
  res@gsnLeftString = "WV flux "+LEV(i)+"hPa"
  plot(i) = gsn_csm_contour_map(wks, dv(0,i,:,:), res)
  plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
  overlay(plot(i), plot_topo(i))
  
  res@gsnRightString = "diff (10~S~-3~N~)"
  ; res@cnLevelSelectionMode = "ExplicitLevels"
  ; res@cnLevels = (/-0.01,-0.007,-0.005,-0.003,-0.001,0.001,0.003,0.005,0.007,0.01/)
  ; ;res@cnFillColors = (/10,18,26,34,42,65,82,90,98,106,114/)
  ; res@cnFillColors = (/114,106,98,90,82,65,42,34,26,18,10/)
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnMaxLevelValF  = 2
  res@cnMinLevelValF  = -2
  res@cnLevelSpacingF = 0.4 
  plot(i+4) = gsn_csm_contour_map(wks, dv_diff(i,:,:), res)
  plot_topo(i+4) = gsn_csm_contour(wks, topo, res_topo)
  overlay(plot(i+4), plot_topo(i+4))
 end do
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
 vecres@vcRefAnnoFontHeightF = 0.02

 do i = 0,nlev-1
 vecres@vcRefMagnitudeF  = 20             ; define vector ref mag
 vecres@vcRefAnnoString1  = "20 g/m/s"
 plotv(i) = gsn_csm_vector(wks, qu(0,i,:,:), qv(0,i,:,:), vecres)
 overlay(plot(i), plotv(i))

 vecres@vcRefMagnitudeF  = 3             ; define vector ref mag
 vecres@vcRefAnnoString1  = "3 g/m/s"
 plotv(i+4) = gsn_csm_vector(wks, qu_diff(i,:,:), qv_diff(i,:,:), vecres)
 overlay(plot(i+4), plotv(i+4))
 end do

 pres = True
 pres@gsnMaximize = True

 gsn_panel(wks, plot, (/2,4/), pres)

end
