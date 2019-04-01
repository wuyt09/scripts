
begin

yrStat = 260
yrLast = 304 

diri1 = "/home/yangsong3/data-model/wuyt/TP_CTRL_B2000_F19G16_CAM4/pro/"
diri2 = "/home/yangsong3/data-model/wuyt/TP_NoTOPO_B2000_F19G16_CAM4/pro/"
diri3 = "/home/yangsong3/data-model/wuyt/TP_NoTOPO_heat_B2000_F19G16_CAM4/pro/"

fili1 = "B2000_f19g16_CP_CTRL.cam.h0.TS.200-349.nc"
fili2 = "TP_NOTOPO.TS.mon.26001-30508.nc"
fili3 = "TP_HEAT.TS.mon.26001-34912.nc"

f1 = addfile(diri1+fili1, "r")
f2 = addfile(diri2+fili2, "r")
f3 = addfile(diri3+fili3, "r")

time1 = f1->time
time2 = f2->time
time3 = f3->time
TIME1 = cd_calendar(time1, -2)/100
TIME2 = cd_calendar(time2, -2)/100
TIME3 = cd_calendar(time3, -2)/100
;print(TIME1(0:40))
;print(TIME2)

iYYYY1 = ind(TIME1.ge.26002 .and. TIME1.le.30501)
iYYYY2 = ind(TIME2.ge.26002 .and. TIME2.le.30501)
iYYYY3 = ind(TIME3.ge.26002 .and. TIME3.le.30501)

t1 = f1->TS(iYYYY1,:,:)
t2 = f2->TS(iYYYY2,:,:)
t3 = f3->TS(iYYYY3,:,:)
printVarSummary(t1)
printVarSummary(t2)
printVarSummary(t3)

Season = (/"DJF","MAM","JJA","SON"/)
tt1 = month_to_seasonN(t1, Season)
tt2 = month_to_seasonN(t2, Season)
tt3 = month_to_seasonN(t3, Season)

t11 = dim_avg_n_Wrap(tt1, 1)
t22 = dim_avg_n_Wrap(tt2, 1)
t33 = dim_avg_n_Wrap(tt3, 1)

t_m = t11 - t33 ; mechanical effect
t_t = t11 - t22 ; total effect
t_h = t_t - t_m ; heating effect
copy_VarCoords(t11, t_m)
copy_VarCoords(t11, t_t)
copy_VarCoords(t11, t_h)
printVarSummary(t_m)

;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_B/ncl/"
figname = "draw_clim_season_TS_heating_mechanical_effect"

lonstr = 40
lonend = 180
latstr = 0
latend = 70

wks = gsn_open_wks("pdf", figpath+figname)
nplot = 12
plot = new(nplot, graphic)

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
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 0.5
res@cnMaxLevelValF = 4.
res@cnMinLevelValF = -4.

do i = 0,dimsizes(Season)-1
res@gsnRightString = Season(i)

res@gsnLeftString = "Total (CTRL - NO_TP)"
plot(i) = gsn_csm_contour_map(wks, t_t(i,:,:), res)

res@gsnLeftString = "Mecha (CTRL - NOTP_HEAT)"
plot(i+4) = gsn_csm_contour_map(wks, t_m(i,:,:), res)

res@gsnLeftString = "Heat (Total - Mecha)"
plot(i+8) = gsn_csm_contour_map(wks, t_h(i,:,:), res)
end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelMainString = "Surface Temperature"
pres@gsnPanelLabelBar = True
pres@pmLabelBarOrthogonalPosF = -0.01

gsn_panel(wks, plot, (/3,4/), pres)
end