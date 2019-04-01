
begin

yrStat = 260
yrLast = 304 

diri1 = "/home/yangsong3/data-model/wuyt/TP_CTRL_B2000_F19G16_CAM4/pro/"
diri2 = "/home/yangsong3/data-model/wuyt/TP_NoTOPO_B2000_F19G16_CAM4/pro/"
diri3 = "/home/yangsong3/data-model/wuyt/TP_NoTOPO_heat_B2000_F19G16_CAM4/pro/"

fili1 = "B2000_f19g16_CP_CTRL.plevel.DTCOND.mon.260-304.nc"
fili2 = "TP_NOTOPO.plevel.DTCOND.mon.260-304.nc"
fili3 = "TP_HEAT.plevel.DTCOND.mon.260-304.nc"

f1 = addfile(diri1+fili1, "r")
f2 = addfile(diri2+fili2, "r")
f3 = addfile(diri3+fili3, "r")

x1 = f1->DTCOND
x2 = f2->DTCOND
x3 = f3->DTCOND
printVarSummary(x1)
printVarSummary(x2)
printVarSummary(x3)

x_m = x1 - x3 ; mechanical effect
x_t = x1 - x2 ; total effect
x_h = x_t - x_m ; heating effect
copy_VarCoords(x1, x_m)
copy_VarCoords(x1, x_t)
copy_VarCoords(x1, x_h)

x11 = dim_avg_n_Wrap(x1, 0)
x22 = dim_avg_n_Wrap(x2, 0)
x33 = dim_avg_n_Wrap(x3, 0)
x_mm = dim_avg_n_Wrap(x_m, 0)
x_tt = dim_avg_n_Wrap(x_t, 0)
x_hh = dim_avg_n_Wrap(x_h, 0)


;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_B/ncl/"
figname = "draw_clim_ann_Ts_heating_mechanical_effect"


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
res@cnLevelSpacingF = 10.
res@cnMaxLevelValF = 320.
res@cnMinLevelValF = 240.

res@gsnLeftString = "CTRL"
plot(0) = gsn_csm_contour_map(wks, x11, res)

res@gsnLeftString = "NO_TP"
plot(2) = gsn_csm_contour_map(wks, x22, res)

res@gsnLeftString = "NOTP_HEAT"
plot(4) = gsn_csm_contour_map(wks, x33, res)

res@cnLevelSpacingF = 0.5
res@cnMaxLevelValF = 4.
res@cnMinLevelValF = -4.

res@gsnLeftString = "Total (CTRL - NO_TP)"
plot(1) = gsn_csm_contour_map(wks, x_tt, res)

res@gsnLeftString = "Mecha (CTRL - NOTP_HEAT)"
plot(3) = gsn_csm_contour_map(wks, x_mm, res)

res@gsnLeftString = "Heat (Total - Mecha)"
plot(5) = gsn_csm_contour_map(wks, x_hh, res)


pres = True
pres@gsnMaximize = True
pres@gsnPanelMainString = "Surface Temperature"

gsn_panel(wks, plot, (/3,2/), pres)
end