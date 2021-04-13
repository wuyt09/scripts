
begin

nlat = 96
nlon = 144
LEV = 200
season = "JJA"

lonstr = 0
lonend = 180
latstr = 0 ;-10
latend = 60

casename = (/"CTRL","TPQ1_double", "NOTP"/) ;
ncase = dimsizes(casename)
year = ispan(6, 30, 1)
nyear = dimsizes(year)

z = new((/ncase,2,nyear,nlat,nlon/), float)
u = z
v = z

do i = 0, ncase-1
  diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(i)+"_F2000_F19G16_CAM4/pro/"
  zname = "F2000_"+casename(i)+".plevel.T.mon.0601-3012.nc"
  uname = "F2000_"+casename(i)+".plevel.U.mon.0601-3012.nc"
  vname = "F2000_"+casename(i)+".plevel.V.mon.0601-3012.nc"
  fz = addfile(diri+zname, "r")
  fu = addfile(diri+uname, "r")
  fv = addfile(diri+vname, "r")
  Z = fz->T(:,{500},:,:)
  U = fu->U(:,{LEV},:,:)
  V = fv->V(:,{LEV},:,:)

  z(i,0,:,:,:) = Z(4::12,:,:)
  z(i,1,:,:,:) = month_to_season(Z, season)
  copy_VarCoords(Z(0,:,:), z(i,0,0,:,:))

  u(i,0,:,:,:) = U(4::12,:,:)
  u(i,1,:,:,:) = month_to_season(U, season)
  ; copy_VarCoords(U(0,:,:), u(i,0,0,:,:))

  v(i,0,:,:,:) = V(4::12,:,:)
  v(i,1,:,:,:) = month_to_season(V, season)
  ; copy_VarCoords(V(0,:,:), v(i,0,0,:,:))

end do
z!0 = "case"
z&case = casename
z!1 = "time"
z&time = (/"May",season/)
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

diff_u(0,:,:,:) = u_ave(1,:,:,:) - u_ave(0,:,:,:)
diff_u(1,:,:,:) = u_ave(0,:,:,:) - u_ave(2,:,:,:)
diff_v(0,:,:,:) = v_ave(1,:,:,:) - v_ave(0,:,:,:)
diff_v(1,:,:,:) = v_ave(0,:,:,:) - v_ave(2,:,:,:)
diff_z(0,:,:,:) = z_ave(1,:,:,:) - z_ave(0,:,:,:)
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

;----------------------------------------------------------------
; tpmaskpath = "/home/yangsong3/wuyt/sysu/scripts/TP_CESM_F/model_need_input_data/"
; tpmaskname = "TPmask_1000m_C20190410.nc"
; ftp = addfile(tpmaskpath+tpmaskname, "r")

; tpmask = ftp->tp
; printVarSummary(tpmask)

; tpmaskpath = "/home/yangsong3/wuyt/data/"
; tpmaskname = "topo.sfc.1x1.nc"
; ftp = addfile(tpmaskpath+tpmaskname, "r")
; tpmask = ftp->topo(::-1,:)
; tpmask = smth9_Wrap(tpmask, 0.5, 0.5, True)

 ff = "~/wuyt/data/TOPO_1500.nc"
 topo_in = addfile(ff,"r")
 topo = topo_in->TP
 xlat = topo_in->LATIXY(:,0)
 xlon = topo_in->LONGXY(0,:)
 topo!0 = "lsmlat"
 topo!1 = "lsmlon"
 topo&lsmlat = xlat
 topo&lsmlon = xlon
 printVarSummary(topo)
 tpmask = topo
;-------------------------------------------------------------------
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_thermal_dynamic_effects/"
figname = "UV200_T500_CTRL_thermal_dynamic"

wks = gsn_open_wks("pdf", figpath+figname)

nplot = 6
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
 res_topo@cnLevels               = (/1500/) ;(/0.9/)
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

res@gsnStringFontHeightF = 0.025
res@tmXBLabelFontHeightF = 0.02
res@tmYLLabelFontHeightF = 0.02

;------------------------------------------

  vcres                         = True             ; plot mods desired
  vcres@gsnFrame                = False            ; so we can draw time stamp
  vcres@gsnDraw                 = False
  vcres@vcRefAnnoOrthogonalPosF = -0.3          ; move ref vector up

  vcres@vcGlyphStyle            = "LineArrow" ; "CurlyVector"    ; turn on curly vectors
  ;vcres@vcMinDistanceF          = 0.017
  ;vcres@vcGlyphStyle            = "LineArrow"
  vcres@vcMinDistanceF          = 0.025
  vcres@vcRefLengthF            = 0.05            ; define length of vec ref
  
  vcres@vcLineArrowThicknessF   = 3.0
  ;vecres@vcLineArrowHeadMinSizeF = 0.008
  ;vecres@vcLineArrowHeadMaxSizeF = 0.3

  vcres@vcRefAnnoString2On = False
  vcres@vcRefAnnoFontHeightF = 0.02
;------------------------------------------
;res@lbLabelBarOn = False

;res@gsnRightString = "mm/day"

  res@cnLevelSelectionMode = "ManualLevels"
  res@cnFillPalette = "CBR_coldhot"
  res@pmLabelBarOrthogonalPosF = 0.2

  
  res@gsnRightString = "CTRL"
  vcres@vcRefMagnitudeF = 30.0
  res@cnLevelSpacingF   = 1.
  res@cnMaxLevelValF    = 270.
  res@cnMinLevelValF    = 260.

  res@gsnLeftString = "May UV200 & T500"
  plot(0) = gsn_csm_contour_map(wks, z_ave(0,0,:,:), res)
  plotv(0) = gsn_csm_vector(wks, u_ave(0,0,:,:), v_ave(0,0,:,:), vcres)
  plot_topo(0) = gsn_csm_contour(wks, tpmask, res_topo)
  overlay(plot(0), plotv(0))
  overlay(plot(0), plot_topo(0))
  res@gsnLeftString = "JJA UV200 & T500"
  plot(3) = gsn_csm_contour_map(wks, z_ave(0,1,:,:), res)
  plotv(3) = gsn_csm_vector(wks, u_ave(0,1,:,:), v_ave(0,1,:,:), vcres)
  plot_topo(3) = gsn_csm_contour(wks, tpmask, res_topo)
  overlay(plot(3), plotv(3))
  overlay(plot(3), plot_topo(3))

  res@gsnRightString = "Thermal effects"
  vcres@vcRefMagnitudeF = 20.0
  res@cnFillPalette = "CBR_coldhot"
  res@cnLevelSpacingF   = 1.
  res@cnMaxLevelValF    = 10.
  res@cnMinLevelValF    = -10.
  res@gsnLeftString = "May UV200 & T500"
  plot(1) = gsn_csm_contour_map(wks, diff_z(0,0,:,:), res)
  plotv(1) = gsn_csm_vector(wks, diff_u(0,0,:,:), diff_v(0,0,:,:), vcres)
  plot_topo(1) = gsn_csm_contour(wks, tpmask, res_topo)
  overlay(plot(1), plotv(1))
  overlay(plot(1), plot_topo(1))
  res@gsnLeftString = "JJA UV200 & T500"
  plot(4) = gsn_csm_contour_map(wks, diff_z(0,1,:,:), res)
  plotv(4) = gsn_csm_vector(wks, diff_u(0,1,:,:), diff_v(0,1,:,:), vcres)
  plot_topo(4) = gsn_csm_contour(wks, tpmask, res_topo)
  overlay(plot(4), plotv(4))
  overlay(plot(4), plot_topo(4))


  res@gsnRightString = "Dynamical effects"
  vcres@vcRefMagnitudeF = 10.0
  res@cnLevelSpacingF   = 1.
  res@cnMaxLevelValF    = 10.
  res@cnMinLevelValF    = -10.
  res@gsnLeftString = "May UV200 & T500"
  plot(2) = gsn_csm_contour_map(wks, diff_z(1,0,:,:), res)
  plotv(2) = gsn_csm_vector(wks, diff_u(1,0,:,:), diff_v(1,0,:,:), vcres)
  plot_topo(2) = gsn_csm_contour(wks, tpmask, res_topo)
  overlay(plot(2), plotv(2))
  overlay(plot(2), plot_topo(2))

  res@gsnLeftString = "JJA UV200 & T500"
  plot(5) = gsn_csm_contour_map(wks, diff_z(1,1,:,:), res)
  plotv(5) = gsn_csm_vector(wks, diff_u(1,1,:,:), diff_v(1,1,:,:), vcres)
  plot_topo(5) = gsn_csm_contour(wks, tpmask, res_topo)
  overlay(plot(5), plotv(5))
  overlay(plot(5), plot_topo(5))


  pres = True
  pres@gsnMaximize = True
  gsn_panel(wks, plot, (/2,3/), pres)
end