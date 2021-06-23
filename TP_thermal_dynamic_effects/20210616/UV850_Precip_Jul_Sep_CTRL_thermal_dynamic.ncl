
begin

nlat = 96
nlon = 144
LEV = 850
season = "JJA"

lonstr = 30
lonend = 150
latstr = 0
latend = 60

casename = (/"CTRL","NOQ1", "NOTP", "NOTP_NOQ1"/) ;
ncase = dimsizes(casename)
year = ispan(6, 30, 1)
nyear = dimsizes(year)

z = new((/ncase,3,nyear,nlat,nlon/), float)
u = z
v = z

do i = 0, ncase-1
  diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(i)+"_F2000_F19G16_CAM4/pro/"
  zname = "F2000_"+casename(i)+".PREC.mon.0601-3012.nc"
  uname = "F2000_"+casename(i)+".plevel.U.mon.0601-3012.nc"
  vname = "F2000_"+casename(i)+".plevel.V.mon.0601-3012.nc"
  fz = addfile(diri+zname, "r")
  fu = addfile(diri+uname, "r")
  fv = addfile(diri+vname, "r")
  Z = fz->PREC
  Z = Z*1000*86400
  Z@units = "mm/day"
  U = fu->U(:,{LEV},:,:)
  V = fv->V(:,{LEV},:,:)

  z(i,0,:,:,:) = Z(6::12,:,:)
  z(i,1,:,:,:) = Z(7::12,:,:)
  z(i,2,:,:,:) = Z(8::12,:,:) ;month_to_season(Z, season)
  copy_VarCoords(Z(0,:,:), z(i,0,0,:,:))

  u(i,0,:,:,:) = U(6::12,:,:)
  u(i,1,:,:,:) = U(7::12,:,:)
  u(i,2,:,:,:) = U(8::12,:,:) ;month_to_season(U, season)
  ; copy_VarCoords(U(0,:,:), u(i,0,0,:,:))

  v(i,0,:,:,:) = V(6::12,:,:)
  v(i,1,:,:,:) = V(7::12,:,:)
  v(i,2,:,:,:) = V(8::12,:,:) ;month_to_season(V, season)
  ; copy_VarCoords(V(0,:,:), v(i,0,0,:,:))

end do
z!0 = "case"
z&case = casename
z!1 = "time"
title = (/"July","August","September"/)
z&time = title
z!2 = "year"
z&year = year
printVarSummary(z)
copy_VarCoords(z, u)
copy_VarCoords(z, v)
printVarSummary(u)

;=========================ttest================================================
u_var = dim_variance_n_Wrap(u, 2)
v_var = dim_variance_n_Wrap(v, 2)
z_var = dim_variance_n_Wrap(z, 2)

u_ave = dim_avg_n_Wrap(u, 2)
v_ave = dim_avg_n_Wrap(v, 2)
z_ave = dim_avg_n_Wrap(z, 2)

diff_u = u_ave(0:ncase-2,:,:,:)
diff_v = v_ave(0:ncase-2,:,:,:)
diff_z = z_ave(0:ncase-2,:,:,:)
prob_u = u_ave(0:ncase-2,:,:,:)
prob_v = v_ave(0:ncase-2,:,:,:)
prob_z = z_ave(0:ncase-2,:,:,:)

; do i = 0, ncase-2
;     diff_u(i,:,:,:) = u_ave(i+1,:,:,:) - u_ave(0,:,:,:)
;     diff_v(i,:,:,:) = v_ave(i+1,:,:,:) - v_ave(0,:,:,:)
;     diff_z(i,:,:,:) = z_ave(i+1,:,:,:) - z_ave(0,:,:,:)

;     prob_u(i,:,:,:) = ttest(u_ave(0,:,:,:), u_var(0,:,:,:), nyear, u_ave(i+1,:,:,:), u_var(i+1,:,:,:), nyear, False, False)
;     prob_v(i,:,:,:) = ttest(v_ave(0,:,:,:), v_var(0,:,:,:), nyear, v_ave(i+1,:,:,:), v_var(i+1,:,:,:), nyear, False, False)
;     prob_z(i,:,:,:) = ttest(z_ave(0,:,:,:), z_var(0,:,:,:), nyear, z_ave(i+1,:,:,:), z_var(i+1,:,:,:), nyear, False, False)
; end do

diff_u(0,:,:,:) = u_ave(0,:,:,:) - u_ave(1,:,:,:)
diff_u(1,:,:,:) = u_ave(0,:,:,:) - u_ave(2,:,:,:)
diff_v(0,:,:,:) = v_ave(0,:,:,:) - v_ave(1,:,:,:)
diff_v(1,:,:,:) = v_ave(0,:,:,:) - v_ave(2,:,:,:)
diff_z(0,:,:,:) = z_ave(0,:,:,:) - z_ave(1,:,:,:)
diff_z(1,:,:,:) = z_ave(0,:,:,:) - z_ave(2,:,:,:)

prob_u(0,:,:,:) = ttest(u_ave(0,:,:,:), u_var(0,:,:,:), nyear, u_ave(1,:,:,:), u_var(1,:,:,:), nyear, False, False)
prob_u(1,:,:,:) = ttest(u_ave(0,:,:,:), u_var(0,:,:,:), nyear, u_ave(2,:,:,:), u_var(2,:,:,:), nyear, False, False)
prob_v(0,:,:,:) = ttest(v_ave(0,:,:,:), v_var(0,:,:,:), nyear, v_ave(1,:,:,:), v_var(1,:,:,:), nyear, False, False)
prob_v(1,:,:,:) = ttest(v_ave(0,:,:,:), v_var(0,:,:,:), nyear, v_ave(2,:,:,:), v_var(2,:,:,:), nyear, False, False)
prob_z(0,:,:,:) = ttest(z_ave(0,:,:,:), z_var(0,:,:,:), nyear, z_ave(1,:,:,:), z_var(1,:,:,:), nyear, False, False)
prob_z(1,:,:,:) = ttest(z_ave(0,:,:,:), z_var(0,:,:,:), nyear, z_ave(2,:,:,:), z_var(2,:,:,:), nyear, False, False)

sig = 0.05
; diff_u = mask(diff_u, (prob_u.le.sig) .and. (prob_v.le.sig), True)
; diff_v = mask(diff_v, (prob_u.le.sig) .and. (prob_v.le.sig), True)
; diff_z = mask(diff_z, prob_z.le.sig, True)
printVarSummary(diff_z)

;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo ;(::-1,:)
printVarSummary(topo)
tpmask = topo

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
 ; tpmask = topo
;-------------------------------------------------------------------
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_thermal_dynamic_effects/20210616/"
figname = "UV850_Precip_Jul_Sep_CTRL_thermal_dynamic"

wks = gsn_open_wks("pdf", figpath+figname)

nplot = 9
plot = new(nplot, graphic)
plotv = plot
plot_topo = plot

;----------------------------------------------------
 res_topo            = True
 res_topo@cnFillOn   = False
 res_topo@cnLinesOn  = True
 res_topo@gsnFrame   = False
 res_topo@gsnDraw    = False
 res_topo@cnLineLabelsOn         = False
 res_topo@cnLineThicknessF       = 4
 res_topo@cnLineDashPattern      = 0
 res_topo@cnLineColor            = "purple"
 res_topo@cnInfoLabelOn          = False
 res_topo@cnLevelSelectionMode   ="ExplicitLevels"
 res_topo@cnLevels               = 2000 ;(/0.9/)
 res_topo@gsnLeftString = ""
 res_topo@gsnRightString = ""
;----------------------------------------------------

res = True
res@gsnDraw  = False
res@gsnFrame = False

res@mpOutlineOn = True
res@mpGeophysicalLineThicknessF = 2
res@mpGeophysicalLineColor = "grey25" ;"red" 
res@mpFillDrawOrder = "PreDraw"
; res@mpFillOn = True

res@mpMinLatF = latstr
res@mpMaxLatF = latend
res@mpMinLonF = lonstr
res@mpMaxLonF = lonend

res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
; res@lbLabelBarOn = False
res@cnInfoLabelOn = False
;res@cnLineDrawOrder = "PostDraw"
; res@cnLabelDrawOrder = "PostDraw"
;res@cnLevelFlag = "LineAndLabel"
;res@cnLineLabelInterval = 1

res@gsnStringFontHeightF = 0.03
res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025

;------------------------------------------

  vcres                         = True             ; plot mods desired
  vcres@gsnFrame                = False            ; so we can draw time stamp
  vcres@gsnDraw                 = False
  vcres@vcRefAnnoOrthogonalPosF = -0.8          ; move ref vector up

  vcres@vcGlyphStyle            = "LineArrow" ; "CurlyVector"    ; turn on curly vectors
  ;vcres@vcMinDistanceF          = 0.017
  ;vcres@vcGlyphStyle            = "LineArrow"
  vcres@vcMinDistanceF          = 0.025
  vcres@vcRefLengthF            = 0.05            ; define length of vec ref
  
  vcres@vcLineArrowThicknessF   = 3.0
  ;vecres@vcLineArrowHeadMinSizeF = 0.008
  ;vecres@vcLineArrowHeadMaxSizeF = 0.3

  vcres@vcRefAnnoString2On = False
  vcres@vcRefAnnoFontHeightF = 0.03
;------------------------------------------
;res@lbLabelBarOn = False

;res@gsnRightString = "mm/day"

  res@cnLevelSelectionMode = "ManualLevels"
  res@cnFillPalette = "CBR_wet"
  res@pmLabelBarOrthogonalPosF = 0.12

  
  res@gsnRightString = "CTRL"
  vcres@vcRefMagnitudeF = 8.0
  res@cnLevelSpacingF   = 2.
  res@cnMaxLevelValF    = 20.
  res@cnMinLevelValF    = 2.

  res@gsnLeftString = title(0);+" UV850 & Precip"
  plot(0) = gsn_csm_contour_map(wks, z_ave(0,0,:,:), res)
  plotv(0) = gsn_csm_vector(wks, u_ave(0,0,:,:), v_ave(0,0,:,:), vcres)
  plot_topo(0) = gsn_csm_contour(wks, tpmask, res_topo)
  overlay(plot(0), plotv(0))
  overlay(plot(0), plot_topo(0))
  res@gsnLeftString = title(1);+" UV850 & Precip"
  plot(3) = gsn_csm_contour_map(wks, z_ave(0,1,:,:), res)
  plotv(3) = gsn_csm_vector(wks, u_ave(0,1,:,:), v_ave(0,1,:,:), vcres)
  plot_topo(3) = gsn_csm_contour(wks, tpmask, res_topo)
  overlay(plot(3), plotv(3))
  overlay(plot(3), plot_topo(3))
  res@gsnLeftString = title(2);+" UV850 & Precip"
  plot(6) = gsn_csm_contour_map(wks, z_ave(0,2,:,:), res)
  plotv(6) = gsn_csm_vector(wks, u_ave(0,2,:,:), v_ave(0,2,:,:), vcres)
  plot_topo(6) = gsn_csm_contour(wks, tpmask, res_topo)
  overlay(plot(6), plotv(6))
  overlay(plot(6), plot_topo(6))

  res@gsnRightString = "Thermal effects"
  vcres@vcRefMagnitudeF = 5.0
  res@cnFillPalette = "CBR_drywet"
  res@cnLevelSpacingF   = 1.
  res@cnMaxLevelValF    = 10.
  res@cnMinLevelValF    = -10.
  res@gsnLeftString = title(0);+" UV850 & Precip"
  plot(1) = gsn_csm_contour_map(wks, diff_z(0,0,:,:), res)
  plotv(1) = gsn_csm_vector(wks, diff_u(0,0,:,:), diff_v(0,0,:,:), vcres)
  plot_topo(1) = gsn_csm_contour(wks, tpmask, res_topo)
  overlay(plot(1), plotv(1))
  overlay(plot(1), plot_topo(1))
  res@gsnLeftString = title(1);+" UV850 & Precip"
  plot(4) = gsn_csm_contour_map(wks, diff_z(0,1,:,:), res)
  plotv(4) = gsn_csm_vector(wks, diff_u(0,1,:,:), diff_v(0,1,:,:), vcres)
  plot_topo(4) = gsn_csm_contour(wks, tpmask, res_topo)
  overlay(plot(4), plotv(4))
  overlay(plot(4), plot_topo(4))
  res@gsnLeftString = title(2);+" UV850 & Precip"
  plot(7) = gsn_csm_contour_map(wks, diff_z(0,2,:,:), res)
  plotv(7) = gsn_csm_vector(wks, diff_u(0,2,:,:), diff_v(0,2,:,:), vcres)
  plot_topo(7) = gsn_csm_contour(wks, tpmask, res_topo)
  overlay(plot(7), plotv(7))
  overlay(plot(7), plot_topo(7))


  res@gsnRightString = "Dynamical effects"
  vcres@vcRefMagnitudeF = 5.0
  res@cnLevelSpacingF   = 1.
  res@cnMaxLevelValF    = 10.
  res@cnMinLevelValF    = -10.
  res@gsnLeftString = title(0);+" UV850 & Precip"
  plot(2) = gsn_csm_contour_map(wks, diff_z(1,0,:,:), res)
  plotv(2) = gsn_csm_vector(wks, diff_u(1,0,:,:), diff_v(1,0,:,:), vcres)
  plot_topo(2) = gsn_csm_contour(wks, tpmask, res_topo)
  overlay(plot(2), plotv(2))
  overlay(plot(2), plot_topo(2))
  res@gsnLeftString = title(1);+" UV850 & Precip"
  plot(5) = gsn_csm_contour_map(wks, diff_z(1,1,:,:), res)
  plotv(5) = gsn_csm_vector(wks, diff_u(1,1,:,:), diff_v(1,1,:,:), vcres)
  plot_topo(5) = gsn_csm_contour(wks, tpmask, res_topo)
  overlay(plot(5), plotv(5))
  overlay(plot(5), plot_topo(5))
  res@gsnLeftString = title(2);+" UV850 & Precip"
  plot(8) = gsn_csm_contour_map(wks, diff_z(1,2,:,:), res)
  plotv(8) = gsn_csm_vector(wks, diff_u(1,2,:,:), diff_v(1,2,:,:), vcres)
  plot_topo(8) = gsn_csm_contour(wks, tpmask, res_topo)
  overlay(plot(8), plotv(8))
  overlay(plot(8), plot_topo(8))


  pres = True
  pres@gsnMaximize = True
  pres@gsnPanelMainString = "UV850 & Precip"
  pres@gsnPanelMainFontHeightF = 0.02
  gsn_panel(wks, plot, (/3,3/), pres)
end