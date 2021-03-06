
begin

lonstr = 0
lonend = 360
latstr = -30
latend = 60

casename = (/"CTRL","NOQ1","NOTP_NOQ1"/) ;,"NOTP"
ncase = dimsizes(casename)
month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
nm = dimsizes(month)
season = (/"MAM","JJA","SON","DJF"/)

LEV = (/200/)
dv = new((/ncase,12,96,144/), float)
dv_season = new((/ncase,4,96,144/), float)


do i = 0, dimsizes(casename)-1
    diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(i)+"_F2000_F19G16_CAM4/pro/"
    filiu = "F2000_"+casename(i)+".plevel.U.mon.0601-3012.nc"
    filiv = "F2000_"+casename(i)+".plevel.V.mon.0601-3012.nc"

    fu = addfile(diri+filiu, "r")
    fv = addfile(diri+filiv, "r")
    u1 = fu->U(:,{LEV},:,:)
    v1 = fv->V(:,{LEV},:,:)
    dv1 = uv2sfvpF(u1, v1)
    dv1 = dv1 * 1.0e-6
    copy_VarCoords(u1, dv1(0,:,:,:))
    printVarSummary(dv1)

    ; dv_Season = month_to_seasonN(dv1(1,:,:,:), season)
    ; dv_season(i,:,:,:) = dim_avg_n_Wrap(dv_Season, 1)
    
    dv(i,:,:,:) = clmMonTLL(dv1(1,:,:,:))

end do
printVarSummary(dv_season)

; dvx = dv_season(0:1,:,:,:)
; dvx(0,:,:,:) = dv_season(0,:,:,:) - dv_season(1,:,:,:)
; dvx(1,:,:,:) = dv_season(1,:,:,:) - dv_season(2,:,:,:)
; printVarSummary(dvx)

dv_warm = dv(:,0,:,:)
dv_warm = dim_avg_n_Wrap(dv(:,3:8,:,:), 1)
dv_cold = dv(:,0,:,:)
dv_cold = dim_avg_n_Wrap(dv(:,9:11,:,:), 1) + dim_avg_n_Wrap(dv(:,0:2,:,:), 1)

;----------------------------------------------------------------
tpmaskpath = "/home/yangsong3/wuyt/sysu/scripts/TP_CESM_F/model_need_input_data/"
tpmaskname = "TPmask_1000m_C20190410.nc"
ftp = addfile(tpmaskpath+tpmaskname, "r")

tpmask = ftp->tp
printVarSummary(tpmask)
;==================================================================
;==================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20191231/"

figname = "clim_season_velocity_potential"

wks = gsn_open_wks("pdf", figpath+figname)
nplot = 4*2
plot = new(nplot, graphic)
plot_topo = plot
plotv = plot
plot_prob = plot


   res_topo            = True
   res_topo@cnFillOn   = False
   res_topo@cnLinesOn  = True
   res_topo@gsnFrame   = False
   res_topo@gsnDraw    = False
   res_topo@cnLineLabelsOn         = False
   res_topo@cnLineThicknessF       = 3
   res_topo@cnLineColor            = "blue"
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
  res@pmLabelBarOrthogonalPosF = 0.15
  

res@lbLabelBarOn = False
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 0.3
res@cnMaxLevelValF = 3.
res@cnMinLevelValF = -3.
res@cnFillPalette = "BlueWhiteOrangeRed"

res@gsnLeftString = "Thermal"
res@gsnRightString = "Warm season"
plot(0) = gsn_csm_contour_map(wks, dv_warm(0,:,:), res)
plot_topo(0) = gsn_csm_contour(wks, tpmask, res_topo)
overlay(plot(0), plot_topo(0))

res@gsnLeftString = "Dynamical"
plot(1) = gsn_csm_contour_map(wks, dv_warm(1,:,:), res)
plot_topo(1) = gsn_csm_contour(wks, tpmask, res_topo)
overlay(plot(1), plot_topo(1))

res@gsnLeftString = "Thermal"
res@gsnRightString = "Cold season"
plot(2) = gsn_csm_contour_map(wks, dv_cold(0,:,:), res)
plot_topo(2) = gsn_csm_contour(wks, tpmask, res_topo)
overlay(plot(2), plot_topo(2))

res@gsnLeftString = "Dynamical"
plot(3) = gsn_csm_contour_map(wks, dv_cold(1,:,:), res)
plot_topo(3) = gsn_csm_contour(wks, tpmask, res_topo)
overlay(plot(3), plot_topo(3))

; do i = 0, dimsizes(season)-1
;   res@gsnRightString = season(i)
;   res@gsnLeftString = "Thermal"
;   plot(i) = gsn_csm_contour_map(wks, dvx(0,i,:,:), res)
;   plot_topo(i) = gsn_csm_contour(wks, tpmask, res_topo)
;   overlay(plot(i), plot_topo(i))
;   res@gsnLeftString = "Dynamical"
;   plot(i+4) = gsn_csm_contour_map(wks, dvx(1,i,:,:), res)
;   plot_topo(i+4) = gsn_csm_contour(wks, tpmask, res_topo)
;   overlay(plot(i+4), plot_topo(i+4))
; end do


pres = True
pres@gsnMaximize = True
pres@gsnPanelMainString = LEV+"hPa_Velocity_Potential"
pres@gsnPanelLabelBar = True
;pres@pmLabelBarOrthogonalPosF = -0.01
pres@lbLabelFontHeightF = 0.01
;gsn_panel(wks, plot, (/2,2/), pres)
gsn_panel(wks, plot, (/2,4/), pres)
end
