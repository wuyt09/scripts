
begin

yrStat = 300
yrLast = 304 

diri1 = "/home/yangsong3/data-model/wuyt/TP_CTRL_B2000_F19G16_CAM4/pro/"
diri2 = "/home/yangsong3/data-model/wuyt/TP_NoTOPO_B2000_F19G16_CAM4/pro/"
diri3 = "/home/yangsong3/data-model/wuyt/TP_NOSH_B2000_F19G16_CAM4/pro/"

fili1 = "B2000_f19g16_CP_CTRL.cam.h0.TEMP.200-349.nc"
fili2 = "TP_NOTOPO.TEMP.mon.25001-30508.nc"
fili3 = "TP_NOSH.TEMP.mon.26001-30712.nc"

f1 = addfile(diri1+fili1, "r")
f2 = addfile(diri2+fili2, "r")
f3 = addfile(diri3+fili3, "r")

time1 = f1->time
time2 = f2->time
time3 = f3->time
TIME1 = cd_calendar(time1, -2)/100
TIME2 = cd_calendar(time2, -2)/100
TIME3 = cd_calendar(time3, -2)/100
;iYYYY1 = ind(TIME1.ge.26002 .and. TIME1.le.30501)
;iYYYY2 = ind(TIME2.ge.26002 .and. TIME2.le.30501)
;iYYYY3 = ind(TIME3.ge.26002 .and. TIME3.le.30501)
iYYYY1 = ind(TIME1.ge.30002 .and. TIME1.le.30501)
iYYYY2 = ind(TIME2.ge.30002 .and. TIME2.le.30501)
iYYYY3 = ind(TIME3.ge.30002 .and. TIME3.le.30501)

sst1 = f1->TEMP(iYYYY1,:,:)
sst2 = f2->TEMP(iYYYY2,:,:)
sst3 = f3->TEMP(iYYYY3,:,:)
sst1@lon2d = f1->ULONG
sst1@lat2d = f1->ULAT
sst2@lon2d = f2->ULONG
sst2@lat2d = f2->ULAT
sst3@lon2d = f3->ULONG
sst3@lat2d = f3->ULAT
printVarSummary(sst1)

Season = (/"DJF","MAM","JJA","SON"/)
tt1 = month_to_seasonN(sst1, Season)
tt2 = month_to_seasonN(sst2, Season)
tt3 = month_to_seasonN(sst3, Season)

t11 = dim_avg_n_Wrap(tt1, 1)
t22 = dim_avg_n_Wrap(tt2, 1)
t33 = dim_avg_n_Wrap(tt3, 1)

t_m = t11 - t33 ; mechanical effect
t_t = t11 - t22 ; total effect
t_h = t_t - t_m ; heating effect
copy_VarMeta(t11, t_m)
copy_VarMeta(t11, t_t)
copy_VarMeta(t11, t_h)
printVarSummary(t_m)

;============================================================

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_B/ncl/"
figname = "draw_clim_season_SST"

lonstr = 40
lonend = 120
latstr = -20
latend = 30

wks = gsn_open_wks("pdf", figpath+figname)
nplot = 12
plot = new(nplot, graphic)


res = True
res@gsnDraw  = False
res@gsnFrame = False
res@gsnAddCyclic = True

res@mpOutlineOn = True
res@mpFillOn = True
res@mpGeophysicalLineThicknessF = 2
;res@mpGeophysicalLineColor = "red" ;"grey25"
;res@mpFillDrawOrder = "PreDraw"

res@mpCenterLonF = 180
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
res@cnLevelSpacingF = 0.1
res@cnMaxLevelValF = 0.6
res@cnMinLevelValF = -0.6
res@cnFillPalette = "temp_19lev"

do i = 0,dimsizes(Season)-1
res@gsnRightString = Season(i)

res@gsnLeftString = "CTRL - NOTP"
plot(i) = gsn_csm_contour_map(wks, t_t(i,:,:), res)

res@gsnLeftString = "CTRL - NOSH"
plot(i+4) = gsn_csm_contour_map(wks, t_m(i,:,:), res)

res@gsnLeftString = "NOSH - NOTP"
plot(i+8) = gsn_csm_contour_map(wks, t_h(i,:,:), res)
end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelMainString = "Surface Sea Potential Temperature"
pres@gsnPanelLabelBar = True
pres@pmLabelBarOrthogonalPosF = -0.01

gsn_panel(wks, plot, (/3,4/), pres)

;----------------------------------------------------------------

res@cnLevelSelectionMode = "ExplicitLevels"
;res@cnLevelSpacingF = 0.2
;res@cnMaxLevelValF = 1.
;res@cnMinLevelValF = -1.
;res@cnFillPalette = "temp_19lev"

do i = 0,dimsizes(Season)-1
res@gsnRightString = Season(i)

res@gsnLeftString = "CTRL"
plot(i) = gsn_csm_contour_map(wks, t11(i,:,:), res)

res@gsnLeftString = "NOTP"
plot(i+4) = gsn_csm_contour_map(wks, t22(i,:,:), res)

res@gsnLeftString = "NOSH"
plot(i+8) = gsn_csm_contour_map(wks, t33(i,:,:), res)
end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelMainString = "Surface Sea Potential Temperature"
pres@gsnPanelLabelBar = True
pres@pmLabelBarOrthogonalPosF = -0.01

gsn_panel(wks, plot, (/3,4/), pres)

end
