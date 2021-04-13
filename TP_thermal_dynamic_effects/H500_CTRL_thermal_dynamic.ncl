
begin

nlat = 96
nlon = 144
LEV = 500
season = "JJA"

lonstr = 0
lonend = 360
latstr = 0 ;-10
latend = 80

casename = (/"CTRL","TPQ1_double", "NOTP"/) ;
ncase = dimsizes(casename)
year = ispan(6, 30, 1)
nyear = dimsizes(year)

z = new((/ncase,2,nyear,nlat,nlon/), float)
do i = 0, ncase-1
  diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(i)+"_F2000_F19G16_CAM4/pro/"
  zname = "F2000_"+casename(i)+".plevel.Z3.mon.0601-3012.nc"
  fz = addfile(diri+zname, "r")
  Z = fz->Z3(:,{LEV},:,:)

  z(i,0,:,:,:) = Z(4::12,:,:)
  z(i,1,:,:,:) = month_to_season(Z, season)
  copy_VarCoords(Z(0,:,:), z(i,0,0,:,:))
end do
z!0 = "case"
z&case = casename
z!1 = "time"
z&time = (/"May",season/)
z!2 = "year"
z&year = year
printVarSummary(z)

;=========================ttest================================================
z_var = dim_variance_n_Wrap(z, 2)
z_ave = dim_avg_n_Wrap(z, 2)

diff_z = z_ave(0:ncase-2,:,:,:)
prob_z = z_ave(0:ncase-2,:,:,:)

diff_z(0,:,:,:) = z_ave(1,:,:,:) - z_ave(0,:,:,:)
diff_z(1,:,:,:) = z_ave(0,:,:,:) - z_ave(2,:,:,:)

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
figname = "H500_CTRL_thermal_dynamic"

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
res@mpCenterLonF = (lonend-lonstr)/2

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

  res@cnLevelSelectionMode = "ManualLevels"
  res@cnFillPalette = "BlueWhiteOrangeRed"
  res@pmLabelBarOrthogonalPosF = 0.3

  
  res@gsnRightString = "CTRL"
  res@cnLevelSpacingF   = 50.
  res@cnMaxLevelValF    = 5940.
  res@cnMinLevelValF    = 5500.

  res@gsnLeftString = "May H500"
  plot(0) = gsn_csm_contour_map(wks, z_ave(0,0,:,:), res)
  plot_topo(0) = gsn_csm_contour(wks, tpmask, res_topo)

  overlay(plot(0), plot_topo(0))
  res@gsnLeftString = "JJA H500"
  plot(3) = gsn_csm_contour_map(wks, z_ave(0,1,:,:), res)
  plot_topo(3) = gsn_csm_contour(wks, tpmask, res_topo)
  overlay(plot(3), plot_topo(3))

  res@gsnRightString = "Thermal effects"
  res@cnFillPalette = "BlueWhiteOrangeRed"
  res@cnLevelSpacingF   = 20.
  res@cnMaxLevelValF    = 100.
  res@cnMinLevelValF    = -100.
  res@gsnLeftString = "May H500"
  plot(1) = gsn_csm_contour_map(wks, diff_z(0,0,:,:), res)
  plot_topo(1) = gsn_csm_contour(wks, tpmask, res_topo)
  overlay(plot(1), plot_topo(1))
  res@gsnLeftString = "JJA H500"
  plot(4) = gsn_csm_contour_map(wks, diff_z(0,1,:,:), res)
  plot_topo(4) = gsn_csm_contour(wks, tpmask, res_topo)
  overlay(plot(4), plot_topo(4))


  res@gsnRightString = "Dynamical effects"
  res@cnLevelSpacingF   = 20.
  res@cnMaxLevelValF    = 100.
  res@cnMinLevelValF    = -100.
  res@gsnLeftString = "May H500"
  plot(2) = gsn_csm_contour_map(wks, diff_z(1,0,:,:), res)
  plot_topo(2) = gsn_csm_contour(wks, tpmask, res_topo)
  overlay(plot(2), plot_topo(2))

  res@gsnLeftString = "JJA H500"
  plot(5) = gsn_csm_contour_map(wks, diff_z(1,1,:,:), res)
  plot_topo(5) = gsn_csm_contour(wks, tpmask, res_topo)
  overlay(plot(5), plot_topo(5))


  pres = True
  pres@gsnMaximize = True
  gsn_panel(wks, plot, (/2,3/), pres)
end