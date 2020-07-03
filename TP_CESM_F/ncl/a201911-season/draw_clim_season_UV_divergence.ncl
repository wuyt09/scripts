
begin

LEV = 200

diri1 = "/home/yangsong3/data-model/wuyt/TP_CTRL_F2000_F19G16_CAM4/pro/"
diri2 = "/home/yangsong3/data-model/wuyt/TP_NOQ1_F2000_F19G16_CAM4/pro/"
diri3 = "/home/yangsong3/data-model/wuyt/TP_NOTP_NOQ1_F2000_F19G16_CAM4/pro/"

filiu1 = "F2000_CTRL.plevel.U.mon.0601-3012.nc"
filiu2 = "F2000_NOQ1.plevel.U.mon.0601-3012.nc"
filiu3 = "F2000_NOTP_NOQ1.plevel.U.mon.0601-3012.nc"

filiv1 = "F2000_CTRL.plevel.V.mon.0601-3012.nc"
filiv2 = "F2000_NOQ1.plevel.V.mon.0601-3012.nc"
filiv3 = "F2000_NOTP_NOQ1.plevel.V.mon.0601-3012.nc"

fu1 = addfile(diri1+filiu1, "r")
fu2 = addfile(diri2+filiu2, "r")
fu3 = addfile(diri3+filiu3, "r")
fv1 = addfile(diri1+filiv1, "r")
fv2 = addfile(diri2+filiv2, "r")
fv3 = addfile(diri3+filiv3, "r")
u1 = fu1->U(:,{LEV},:,:)
u2 = fu2->U(:,{LEV},:,:)
u3 = fu3->U(:,{LEV},:,:)
v1 = fv1->V(:,{LEV},:,:)
v2 = fv2->V(:,{LEV},:,:)
v3 = fv3->V(:,{LEV},:,:)
printVarSummary(u1)

lat = fu1->lat
lon = fu1->lon
;=========================================
div1 = uv2dv_cfd(u1, v1, lat, lon, 1)
div2 = uv2dv_cfd(u2, v2, lat, lon, 1)
div3 = uv2dv_cfd(u3, v3, lat, lon, 1)
div1 = div1 * 1e+6 
div2 = div2 * 1e+6 
div3 = div3 * 1e+6 
copy_VarCoords(u1, div1)
copy_VarCoords(u2, div2)
copy_VarCoords(u3, div3)
printVarSummary(div1)

;=========================================

Season = (/"MAM","JJA","SON","DJF"/)
x1 = month_to_seasonN(div1, Season)
x11 = dim_avg_n_Wrap(x1, 1)
x2 = month_to_seasonN(div2, Season)
x22 = dim_avg_n_Wrap(x2, 1)
x3 = month_to_seasonN(div3, Season)
x33 = dim_avg_n_Wrap(x3, 1)

x_h = x11 - x22
x_m = x22 - x33
copy_VarCoords(x11, x_h)
copy_VarCoords(x22, x_m)
printVarSummary(x_h)

;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/"
figname = "draw_clim_season_UV_divergence"

lonstr = 40 ;20
lonend = 180 ;220
latstr = -10
latend = 80

wks = gsn_open_wks("pdf", figpath+figname)
nplot = 12
plot = new(nplot, graphic)


res = True
res@gsnDraw  = False
res@gsnFrame = False

res@mpOutlineOn = True
res@mpFillOn = True
res@mpGeophysicalLineThicknessF = 2
res@mpGeophysicalLineColor = "grey45" ;"grey25"
;res@mpFillDrawOrder = "PreDraw"

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
res@cnFillPalette = "BlueWhiteOrangeRed"
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 1.
res@cnMaxLevelValF = 10.
res@cnMinLevelValF = -10.

do i = 0, dimsizes(Season)-1
res@gsnRightString = Season(i)
res@gsnLeftString = "CTRL"
plot(i) = gsn_csm_contour_map(wks, x11(i,:,:), res)
;plot_topo(i) = gsn_csm_contour(wks, tpmask, res_topo)
;plotv(i) = gsn_csm_vector(wks, uu1(i,:,:), vv1(i,:,:), vcres)
;overlay(plot(i), plot_topo(i))
;overlay(plot(i), plotv(i))

res@gsnLeftString = "NOQ1"
plot(i+4) = gsn_csm_contour_map(wks, x22(i,:,:), res)


res@gsnLeftString = "NOTP_NOQ1"
plot(i+8) = gsn_csm_contour_map(wks, x33(i,:,:), res)

end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelMainString = "UV"+LEV+"_divergence_climatology"
pres@gsnPanelLabelBar = True
;pres@pmLabelBarOrthogonalPosF = -0.01

gsn_panel(wks, plot, (/3,4/), pres)


res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 0.5
res@cnMaxLevelValF = 3.
res@cnMinLevelValF = -3.

do i = 0, dimsizes(Season)-1

res@gsnRightString = Season(i)
res@gsnLeftString = "thermal"
plot(i) = gsn_csm_contour_map_ce(wks, x_h(i,:,:), res)

res@gsnLeftString = "dynamical"
plot(i+4) = gsn_csm_contour_map_ce(wks, x_m(i,:,:), res)

end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelMainString = "UV"+LEV+"_divergence_difference"
pres@gsnPanelLabelBar = True
pres@pmLabelBarOrthogonalPosF = -0.01

gsn_panel(wks, plot, (/2,4/), pres)


end