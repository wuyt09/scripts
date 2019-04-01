
begin

yrStat = 260
yrLast = 304 

diri = "/home/yangsong3/data-model/wuyt/TP_NoTOPO_heat_B2000_F19G16_CAM4/pro/"

fili1 = "TP_HEAT.plevel.DTCOND.mon.260-304.nc"
fili2 = "TP_HEAT.plevel.DTV.mon.260-304.nc"

f1 = addfile(diri+fili1, "r")
f2 = addfile(diri+fili2, "r")
x1 = f1->DTCOND
x2 = f2->DTV
printVarSummary(x1)
x1 = x1*86400
x2 = x2*86400

total = x1
total = x1+x2
printVarSummary(total)

x11 = dim_avg_n_Wrap(x1, 0)
x22 = dim_avg_n_Wrap(x2, 0)
totall = dim_avg_n_Wrap(total, 0)

;=========================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_B/ncl/"
figname = "draw_clim_ann_NOTP_HEAT_DTCOND_DTV_spatial"

wks = gsn_open_wks("pdf", figpath+figname)

nplot = 6
plot = new(nplot, graphic)

res = True
res@gsnDraw  = False
res@gsnFrame = False

res@mpOutlineOn = True
res@mpFillOn = True
res@mpGeophysicalLineThicknessF = 2

res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False

res@gsnStringFontHeightF = 0.025


res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 0.5
res@cnMaxLevelValF = 5.
res@cnMinLevelValF = -5.

res@gsnLeftString = "DTCOND_850hPa"
plot(0) = gsn_csm_contour_map(wks, x11({850},:,:), res)
res@gsnLeftString = "DTV_850hPa"
plot(1) = gsn_csm_contour_map(wks, x22({850},:,:), res)

res@gsnLeftString = "DTCOND_500hPa"
plot(2) = gsn_csm_contour_map(wks, x11({500},:,:), res)
res@gsnLeftString = "DTV_850hPa"
plot(3) = gsn_csm_contour_map(wks, x22({500},:,:), res)

res@gsnLeftString = "DTCOND_200hPa"
plot(4) = gsn_csm_contour_map(wks, x11({200},:,:), res)
res@gsnLeftString = "DTV_200hPa"
plot(5) = gsn_csm_contour_map(wks, x22({200},:,:), res)

;res@cnLevelSpacingF = 0.5
;res@cnMaxLevelValF = 4.
;res@cnMinLevelValF = -4.

pres = True
pres@gsnMaximize = True
pres@gsnPanelMainString = " test "

gsn_panel(wks, plot, (/3,2/), pres)

end