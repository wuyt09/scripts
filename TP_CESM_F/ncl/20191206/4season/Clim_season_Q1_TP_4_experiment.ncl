
begin

casename = (/"CTRL","NOQ1"/) ;"NOTP",,"NOTP_NOQ1"
ncase = dimsizes(casename)
q1 = new((/ncase,4,96,144/), float)

season = (/"MAM","JJA","SON","DJF"/)
Cpair = 1004.07
G = 9.8 

do i = 0, ncase-1
  diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(i)+"_F2000_F19G16_CAM4/pro/"
  filiu = "F2000_"+casename(i)+".plevel.DTCOND.mon.0601-3012.nc"
  filiv = "F2000_"+casename(i)+".plevel.DTV.mon.0601-3012.nc"
  filis = "F2000_"+casename(i)+".plevel.QRS.mon.0601-3012.nc"
  filil = "F2000_"+casename(i)+".plevel.QRL.mon.0601-3012.nc"
  filips = "F2000_"+casename(i)+".PS.mon.0601-3012.nc"

  fu = addfile(diri+filiu, "r")
  fv = addfile(diri+filiv, "r")
  fs = addfile(diri+filis, "r")
  fl = addfile(diri+filil, "r")
  fps = addfile(diri+filips, "r")

  U = fu->DTCOND
  V = fv->DTV
  S = fs->QRS
  L = fl->QRL
  PS = fps->PS

  Q1 = U
  Q1 = (U+V+S+L)*Cpair/G

  plev = U&lev*100
  plev!0 = "lev"
  plev&lev = U&lev
  plev@units = "Pa"
  psfc = PS
  ptop = min(plev)
  dp = dpres_plevel(plev, psfc, ptop, 0)

  Q1dp = Q1
  Q1dp = Q1*dp
  q1dp = dim_sum_n_Wrap(Q1dp, 1)
  copy_VarCoords(Q1(:,0,:,:), q1dp)
  q1x = month_to_seasonN(q1dp, season)
  q1(i,:,:,:) = dim_avg_n_Wrap(q1x, 1)
  print("end of "+casename(i))
end do
printVarSummary(q1)

; u_h = uu1 - uu3
; u_m = uu1 - uu2
; u_h = u_t - u_m
; copy_VarCoords(uu1, u_m)
; copy_VarCoords(uu1, u_t)
; copy_VarCoords(uu1, u_h)

;==================================================================

tpmaskpath = "/home/yangsong3/wuyt/sysu/scripts/TP_CESM_F/model_need_input_data/"
tpmaskname = "TPmask_1000m_C20190410.nc"
ftp = addfile(tpmaskpath+tpmaskname, "r")

tpmask = ftp->tp
printVarSummary(tpmask)
;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20191206/4season/"
figname = "Clim_season_Q1_TP_4_experiment"

; lonstr = 0
; lonend = 360
; latstr = -30
; latend = 80

lonstr = 50
lonend = 120
latstr = 20
latend = 50


wks = gsn_open_wks("pdf", figpath+figname)
nplot = 4*ncase
plot = new(nplot, graphic)
plot_topo = plot

 res_topo            = True
 res_topo@cnFillOn   = False
 res_topo@cnLinesOn  = True
 res_topo@gsnFrame   = False
 res_topo@gsnDraw    = False
 res_topo@cnLineLabelsOn         = False
 res_topo@cnLineThicknessF       = 2
 res_topo@cnLineColor            = "grey80"
 res_topo@cnLineDashPattern      = 12
 res_topo@cnInfoLabelOn          = False
 res_topo@cnLevelSelectionMode   ="ExplicitLevels"
 res_topo@cnLevels               = (/0.9/)
 res_topo@gsnLeftString = ""
 res_topo@gsnRightString = ""

res = True
res@gsnDraw  = False
res@gsnFrame = False

res@mpOutlineOn = True
res@mpFillOn = True
res@mpGeophysicalLineThicknessF = 2
res@mpGeophysicalLineColor = "grey45" ;"grey25"
res@mpFillDrawOrder = "PreDraw"

res@mpMinLatF = latstr
res@mpMaxLatF = latend
res@mpMinLonF = lonstr
res@mpMaxLonF = lonend
res@mpCenterLonF = (lonstr+lonend)/2

res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False

res@gsnStringFontHeightF = 0.035
res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025

res@lbLabelBarOn = False
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 30.
res@cnMaxLevelValF = 300.
res@cnMinLevelValF = -300.
res@cnFillPalette = "BlueWhiteOrangeRed"

leftsti = (/"Ctrl run","NH run","NTNH run"/)
do i = 0, ncase-1
  do j = 0, dimsizes(season)-1
    res@gsnLeftString = leftsti(i)
    res@gsnRightString = season(j)
    plot(j+i*4) = gsn_csm_contour_map(wks, q1(i,j,:,:), res)
    plot_topo(j+i*4) = gsn_csm_contour(wks, tpmask, res_topo)
    overlay(plot(j+i*4), plot_topo(j+i*4))
  end do
end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelMainString = "Q1_climatology"
pres@gsnPanelLabelBar = True
;pres@pmLabelBarOrthogonalPosF = -0.01
pres@lbLabelFontHeightF = 0.01

gsn_panel(wks, plot, (/ncase,4/), pres)
exit
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 5.
res@cnMaxLevelValF = 50.
res@cnMinLevelValF = -50.
res@cnFillPalette = "BlueWhiteOrangeRed"

thermal = q1(0,:,:,:)
dynamic = q1(0,:,:,:)
thermal = q1(0,:,:,:) - q1(1,:,:,:)
dynamic = q1(1,:,:,:) - q1(2,:,:,:)

do j = 0, dimsizes(season)-1
  res@gsnRightString = season(j)
  res@gsnLeftString = "Thermal"
  plot(2*j) = gsn_csm_contour_map(wks, thermal(j,:,:), res)
  plot_topo(2*j) = gsn_csm_contour(wks, tpmask, res_topo)
  overlay(plot(2*j), plot_topo(2*j))

  res@gsnLeftString = "Dynamical"
  plot(2*j+1) = gsn_csm_contour_map(wks, dynamic(j,:,:), res)
  plot_topo(2*j+1) = gsn_csm_contour(wks, tpmask, res_topo)
  overlay(plot(2*j+1), plot_topo(2*j+1))
end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelMainString = "Q1_difference"
pres@gsnPanelLabelBar = True
;pres@pmLabelBarOrthogonalPosF = -0.01


gsn_panel(wks, plot, (/4,2/), pres)

end