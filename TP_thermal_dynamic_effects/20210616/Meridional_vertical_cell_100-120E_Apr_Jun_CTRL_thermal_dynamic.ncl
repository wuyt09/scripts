
begin

nlat = 96
nlon = 144
nlev = 37
season = "JJA"

lonstr = 0
lonend = 180
latstr = 0 ;-10
latend = 60

lonW = 100
lonE = 120

casename = (/"CTRL","NOQ1", "NOTP", "NOTP_NOQ1"/) ;
ncase = dimsizes(casename)
year = ispan(6, 30, 1)
nyear = dimsizes(year)

u = new((/ncase,3,nyear,nlev,nlat/), float)
v = u

do i = 0, ncase-1
  diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(i)+"_F2000_F19G16_CAM4/pro/"
  uname = "F2000_"+casename(i)+".plevel.OMEGA.mon.0601-3012.nc"
  vname = "F2000_"+casename(i)+".plevel.V.mon.0601-3012.nc"
  fu = addfile(diri+uname, "r")
  fv = addfile(diri+vname, "r")
  U_med = fu->OMEGA
  V_med = fv->V

  U = dim_avg_n_Wrap(U_med(:,:,:,{lonW:lonE}), 3)
  V = dim_avg_n_Wrap(V_med(:,:,:,{lonW:lonE}), 3)

  u(i,0,:,:,:) = U(3::12,:,:)
  u(i,1,:,:,:) = U(4::12,:,:)
  u(i,2,:,:,:) = U(5::12,:,:) ;month_to_season(U, season)
  copy_VarCoords(U(0,:,:), u(i,0,0,:,:))

  v(i,0,:,:,:) = V(3::12,:,:)
  v(i,1,:,:,:) = V(4::12,:,:)
  v(i,2,:,:,:) = V(5::12,:,:) ;month_to_season(V, season)
  ; copy_VarCoords(V(0,:,:), v(i,0,0,:,:))

end do
u!0 = "case"
u&case = casename
u!1 = "time"
title = (/"April","May","June"/)
u&time = title
u!2 = "year"
u&year = year
u = u*(-50.)
printVarSummary(u)
copy_VarCoords(u, v)

;=========================ttest================================================
u_var = dim_variance_n_Wrap(u, 2)
v_var = dim_variance_n_Wrap(v, 2)

u_ave = dim_avg_n_Wrap(u, 2)
v_ave = dim_avg_n_Wrap(v, 2)

diff_u = u_ave(0:ncase-2,:,:,:)
diff_v = v_ave(0:ncase-2,:,:,:)
prob_u = u_ave(0:ncase-2,:,:,:)
prob_v = v_ave(0:ncase-2,:,:,:)

diff_u(0,:,:,:) = u_ave(0,:,:,:) - u_ave(1,:,:,:)
diff_u(1,:,:,:) = u_ave(0,:,:,:) - u_ave(2,:,:,:)
diff_v(0,:,:,:) = v_ave(0,:,:,:) - v_ave(1,:,:,:)
diff_v(1,:,:,:) = v_ave(0,:,:,:) - v_ave(2,:,:,:)

prob_u(0,:,:,:) = ttest(u_ave(0,:,:,:), u_var(0,:,:,:), nyear, u_ave(1,:,:,:), u_var(1,:,:,:), nyear, False, False)
prob_u(1,:,:,:) = ttest(u_ave(0,:,:,:), u_var(0,:,:,:), nyear, u_ave(2,:,:,:), u_var(2,:,:,:), nyear, False, False)
prob_v(0,:,:,:) = ttest(v_ave(0,:,:,:), v_var(0,:,:,:), nyear, v_ave(1,:,:,:), v_var(1,:,:,:), nyear, False, False)
prob_v(1,:,:,:) = ttest(v_ave(0,:,:,:), v_var(0,:,:,:), nyear, v_ave(2,:,:,:), v_var(2,:,:,:), nyear, False, False)

sig = 0.05
; diff_u = mask(diff_u, (prob_u.le.sig) .and. (prob_v.le.sig), True)
; diff_v = mask(diff_v, (prob_u.le.sig) .and. (prob_v.le.sig), True)
; diff_z = mask(diff_z, prob_z.le.sig, True)
printVarSummary(diff_u)

;----------------------------------------------------------------
tpmaskpath = "/home/yangsong3/wuyt/data/"
tpmaskname = "topo.sfc.1x1.nc"
ftp = addfile(tpmaskpath+tpmaskname, "r")
topo = ftp->topo(::-1,:)
xlat = topo&lat
; tpmask = smth9_Wrap(tpmask, 0.5, 0.5, True)

 ; ff = "~/wuyt/data/TOPO_1500.nc"
 ; topo_in = addfile(ff,"r")
 ; topo = topo_in->TP
 ; xlat = topo_in->LATIXY(:,0)
 ; xlon = topo_in->LONGXY(0,:)
 ; topo!0 = "lsmlat"
 ; topo!1 = "lsmlon"
 ; topo&lsmlat = xlat
 ; topo&lsmlon = xlon
 ; printVarSummary(topo)
 
topo_m = dim_avg_n_Wrap(topo(:,{lonW:lonE}), 1)
topo_merid = 1013.25*(1-topo_m*0.0065/288.15)^5.25145
copy_VarCoords(topo(:,:), topo_merid)
;-------------------------------------------------------------------
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_thermal_dynamic_effects/20210616/"
figname = "Meridional_vertical_cell_100-120E_Apr_Jun_CTRL_thermal_dynamic"

wks = gsn_open_wks("pdf", figpath+figname)

nplot = 9
plot = new(nplot, graphic)
plotv = plot
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
;=====================================================

    res = True
    res@gsnDraw  = False
    res@gsnFrame = False
    res@vpHeightF = 0.4
    res@vpWidthF  = 0.7

    res@cnFillOn = True
    res@cnLinesOn = False
    res@cnLineLabelsOn = False
    res@cnInfoLabelOn = False

    res@gsnStringFontHeightF = 0.035
    res@tmXBLabelFontHeightF = 0.025
    res@tmYLLabelFontHeightF = 0.025

    ;res@pmLabelBarOrthogonalPosF = 0.12
    res@lbLabelBarOn = False
    ;res@lbLabelFontHeightF = 0.02
    res@cnFillPalette = "BlueWhiteOrangeRed"
    res@tiYAxisString = ""
    ; res@gsnpre
    res@cnFillDrawOrder = "PreDraw"
    res@gsnPresHgtHeightLabelOn = False


 ;------------------------------------------
     res@trYMinF   = 100
     res@trYMaxF   = 1000
      
     res@tmYLMode      = "Explicit"
     res@tmYLValues    = (/"1000","850","750","500","300","200","100"/)
     res@tmYLLabels    = (/"1000","850","750","500","300","200","100"/)

     res@vcRefAnnoOrthogonalPosF =  -0.5 ;-1.0             ; move ref vector up
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
     res@vcRefAnnoFontHeightF = 0.025

;------------------------------------------

res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 1
res@cnMaxLevelValF = 5
res@cnMinLevelValF = -5
res@cnFillPalette = "BlueWhiteOrangeRed"
res@vcRefMagnitudeF         = 5.  

  res@gsnRightString = "CTRL"
  res@gsnLeftString = title(0)
  plot(0) = gsn_csm_pres_hgt_vector(wks, u_ave(0,0,{1000:100},{0:60}),\
 v_ave(0,0,{1000:100},{0:60}), u_ave(0,0,{1000:100},{0:60}), res)

  plot_topo(0) = gsn_csm_xy(wks, xlat, topo_merid, res_topo)
  overlay(plot(0), plot_topo(0))

  res@gsnLeftString = title(1)
  plot(3) = gsn_csm_pres_hgt_vector(wks, u_ave(0,1,{1000:100},{0:60}),\
 v_ave(0,1,{1000:100},{0:60}), u_ave(0,1,{1000:100},{0:60}), res)

  plot_topo(3) = gsn_csm_xy(wks, xlat, topo_merid, res_topo)
  overlay(plot(3), plot_topo(3))

  res@gsnLeftString = title(2)
  plot(6) = gsn_csm_pres_hgt_vector(wks, u_ave(0,2,{1000:100},{0:60}),\
 v_ave(0,2,{1000:100},{0:60}), u_ave(0,2,{1000:100},{0:60}), res)

  plot_topo(6) = gsn_csm_xy(wks, xlat, topo_merid, res_topo)
  overlay(plot(6), plot_topo(6))

  res@gsnRightString = "Thermal effects"
  ; res@cnFillPalette = "CBR_coldhot"
  ; res@cnLevelSpacingF   = 1.
  ; res@cnMaxLevelValF    = 10.
  ; res@cnMinLevelValF    = -10.
  res@vcRefMagnitudeF         = 2.  
  res@gsnLeftString = title(0)
  plot(1) = gsn_csm_pres_hgt_vector(wks, diff_u(0,0,{1000:100},{0:60}),\
 diff_v(0,0,{1000:100},{0:60}), diff_u(0,0,{1000:100},{0:60}), res)
  plot_topo(1) = gsn_csm_xy(wks, xlat, topo_merid, res_topo)
  overlay(plot(1), plot_topo(1))

  res@gsnLeftString = title(1)
  plot(4) = gsn_csm_pres_hgt_vector(wks, diff_u(0,1,{1000:100},{0:60}),\
 diff_v(0,1,{1000:100},{0:60}), diff_u(0,1,{1000:100},{0:60}), res)
  plot_topo(4) = gsn_csm_xy(wks, xlat, topo_merid, res_topo)
  overlay(plot(4), plot_topo(4))

  res@gsnLeftString = title(2)
  plot(7) = gsn_csm_pres_hgt_vector(wks, diff_u(0,2,{1000:100},{0:60}),\
 diff_v(0,2,{1000:100},{0:60}), diff_u(0,2,{1000:100},{0:60}), res)
  plot_topo(7) = gsn_csm_xy(wks, xlat, topo_merid, res_topo)
  overlay(plot(7), plot_topo(7))


  res@gsnRightString = "Dynamical effects"
  ; vcres@vcRefMagnitudeF = 10.0
  ; res@cnLevelSpacingF   = 1.
  ; res@cnMaxLevelValF    = 10.
  ; res@cnMinLevelValF    = -10.
  res@gsnLeftString = title(0)
  plot(2) = gsn_csm_pres_hgt_vector(wks, diff_u(1,0,{1000:100},{0:60}),\
 diff_v(1,0,{1000:100},{0:60}), diff_u(1,0,{1000:100},{0:60}), res)
  plot_topo(2) = gsn_csm_xy(wks, xlat, topo_merid, res_topo)
  overlay(plot(2), plot_topo(2))

  res@gsnLeftString = title(1)
  plot(5) = gsn_csm_pres_hgt_vector(wks, diff_u(1,1,{1000:100},{0:60}),\
 diff_v(1,1,{1000:100},{0:60}), diff_u(1,1,{1000:100},{0:60}), res)
  plot_topo(5) = gsn_csm_xy(wks, xlat, topo_merid, res_topo)
  overlay(plot(5), plot_topo(5))

  res@gsnLeftString = title(2)
  plot(8) = gsn_csm_pres_hgt_vector(wks, diff_u(1,2,{1000:100},{0:60}),\
 diff_v(1,2,{1000:100},{0:60}), diff_u(1,2,{1000:100},{0:60}), res)
  plot_topo(8) = gsn_csm_xy(wks, xlat, topo_merid, res_topo)
  overlay(plot(8), plot_topo(8))


  pres = True
  pres@gsnMaximize = True
  pres@gsnPanelLabelBar = True
  gsn_panel(wks, plot, (/3,3/), pres)
end