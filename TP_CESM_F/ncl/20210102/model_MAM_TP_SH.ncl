
begin

LEV = 200
nlat = 96
nlon = 144

; lonstr = 108
; lonend = 120
; latstr = 20
; latend = 31
lonstr = 60
lonend = 120
latstr = 10
latend = 50

casename = (/"CTRL","TPQ1_double"/) ;"NOTP"
ncase = dimsizes(casename)
year = ispan(6, 30, 1)
nyear = dimsizes(year)
month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
nmonth = dimsizes(month)
season = "MAM"

LEV = 850
z = new((/ncase,nyear,96,144/), float)
u = z
v = z

do i = 0, ncase-1
  diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(i)+"_F2000_F19G16_CAM4/pro/"
  zname = "F2000_"+casename(i)+".SHFLX.mon.0601-3012.nc"
  fz = addfile(diri+zname, "r")
  Z = fz->SHFLX
  z(i,:,:,:) = month_to_season(Z, season)
end do
z!0 = "case"
z&case = casename
; z!1 = "year"
; z&year = year
; z!2 = "month"
; z&month = month
printVarSummary(z)
copy_VarCoords(z, u)
copy_VarCoords(z, v)

;------------ttest--------------------------------------
z_var = dim_variance_n_Wrap(z, 1)

z_ave = dim_avg_n_Wrap(z, 1)
probz = ttest(z_ave(0,:,:), z_var(0,:,:), nyear, z_ave(1,:,:), z_var(1,:,:), nyear, False, False)

diff_z = z_ave(0,:,:)
diff_z = z_ave(1,:,:) - z_ave(0,:,:)
sig = 0.05
; diff_z = mask(diff_z, probz.lt.sig, True)
printVarSummary(diff_z)
;----------------------------------------------------------------
tpmaskpath = "/home/yangsong3/wuyt/sysu/scripts/TP_CESM_F/model_need_input_data/"
tpmaskname = "TPmask_1500m_C20200713.nc"
ftp = addfile(tpmaskpath+tpmaskname, "r")

tpmask = ftp->tp
printVarSummary(tpmask)
;-------------------------------------------------------------------
;-------------------------------------------------------------------


figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20210102/"
figname = "model_"+season+"_TP_SH"

wks = gsn_open_wks("pdf", figpath+figname)

nplot = 3
plot = new(nplot, graphic)
plotv = plot
plotc = plot
plot_topo = plot
;----------------------------------------------------
 res_topo            = True
 res_topo@cnFillOn   = False
 res_topo@cnLinesOn  = True
 res_topo@gsnFrame   = False
 res_topo@gsnDraw    = False
 res_topo@cnLineLabelsOn         = False
 res_topo@cnLineThicknessF       = 4
 res_topo@cnLineDashPattern      = 16
 res_topo@cnLineColor            = "gray"
 res_topo@cnInfoLabelOn          = False
 res_topo@cnLevelSelectionMode   ="ExplicitLevels"
 res_topo@cnLevels               = (/0.9/)
 res_topo@gsnLeftString = ""
 res_topo@gsnRightString = ""
;----------------------------------------------------

res = True
res@gsnDraw  = False
res@gsnFrame = False
res@gsnMaximize = True

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

res@gsnStringFontHeightF = 0.035
res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025

;res@lbLabelBarOn = False

  res@cnLevelSelectionMode = "ManualLevels"
  res@cnLevelSpacingF = 20.
  res@cnMaxLevelValF = 200.
  res@cnMinLevelValF = -200.
  res@cnFillPalette = "cmp_flux"
  res@gsnLeftString = "(a) CTRL"
plot(0) = gsn_csm_contour_map(wks, z_ave(0,:,:), res)
plot_topo(0) = gsn_csm_contour(wks, tpmask, res_topo)
overlay(plot(0), plot_topo(0))

  res@gsnLeftString = "(b) TPQ1_double"
plot(1) = gsn_csm_contour_map(wks, z_ave(1,:,:), res)
plot_topo(1) = gsn_csm_contour(wks, tpmask, res_topo)
overlay(plot(1), plot_topo(1))





  res@cnLevelSelectionMode = "ManualLevels"
  res@cnLevelSpacingF = 5.
  res@cnMaxLevelValF = 30.
  res@cnMinLevelValF = -30.
  res@cnFillPalette = "cmp_flux"
  res@gsnLeftString = "(c) b - a"
plot(2) = gsn_csm_contour_map(wks, diff_z, res)
plot_topo(2) = gsn_csm_contour(wks, tpmask, res_topo)
overlay(plot(2), plot_topo(2))

pres = True
pres@gsnMaximize = True
pres@gsnPanelMainString = season+"_sensible_flux"
; pres@gsnPanelLabelBar = True
;pres@pmLabelBarOrthogonalPosF = -0.01

gsn_panel(wks, plot, (/1,3/), pres)

end