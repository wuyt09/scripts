
begin

season = (/"MAM","JJA","SON","DJF"/)

diri1 = "/home/yangsong3/data-model/wuyt/TP_CTRL_B2000_F19G16_CAM4/pro/"
diri2 = "/home/yangsong3/data-model/wuyt/TP_NOSH_B2000_F19G16_CAM4/pro/"

fili1 = "B2000_f19g16_CP_CTRL.cam.h0.PRECL.260-304.nc"
fili2 = "TP_NOSH.PRECL.mon.26001-30712.nc"

f1 = addfile(diri1+fili1, "r")
f2 = addfile(diri2+fili2, "r")

p1 = f1->PRECL
p2 = f2->PRECL
printVarSummary(p1)

p1 = p1*1000*86400
p2 = p2*1000*86400 ; m/s to mm/day

p1Sea = month_to_seasonN(p1, season)
p2Sea = month_to_seasonN(p2, season)
printVarSummary(p1Sea)

p1sea = dim_avg_n_Wrap(p1Sea, 1)
p2sea = dim_avg_n_Wrap(p2Sea, 1)
printVarSummary(p1sea)

;=====================================================================================

lonstr = 60
lonend = 120
latstr = 20
latend = 45

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_B/ncl/"
figname = "draw_clim_season_CTRL_NOSH_PRECL"

wks = gsn_open_wks("pdf", figpath+figname)
nplot = 4
plot = new(nplot, graphic)
plot_topo = 4

res = True
res@gsnDraw  = False
res@gsnFrame = False

res@mpOutlineOn = True
res@mpFillOn = True
res@mpGeophysicalLineThicknessF = 2

res@mpMinLatF = latstr
res@mpMaxLatF = latend
res@mpMinLonF = lonstr
res@mpMaxLonF = lonend

res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False

res@gsnStringFontHeightF = 0.035
res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025

res@lbLabelBarOn = False
res@cnFillPalette = "CBR_wet"
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 1
res@cnMaxLevelValF = 11.
res@cnMinLevelValF = 1.

do i = 0, dimsizes(season)-1
    res@gsnLeftString = season(i)
    plot(i) = gsn_csm_contour_map(wks, p1sea(i,:,:), res)
end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = "CTRL"
gsn_panel(wks, plot, (/2,2/), pres)

do i = 0, dimsizes(season)-1
    res@gsnLeftString = season(i)
    plot(i) = gsn_csm_contour_map(wks, p2sea(i,:,:), res)
end do
pres@gsnPanelMainString = "NOSH"
gsn_panel(wks, plot, (/2,2/), pres)


end