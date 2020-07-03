
begin

;======================================================================
diriera = "/home/yangsong3/data-observation/wuyt/ERAI/monthly/1x1/"
filiera = "t2m.mon.mean.7901.1804.nc"
fera = addfile(diriera+filiera, "r")
TIME1 = fera->time
YYYY1 = cd_calendar(TIME1, -1)/100
iYYYY1 = ind(YYYY1.ge.1979 .and. YYYY1.le.2017)
TS_era = short2flt(fera->t2m(iYYYY1,::-1,:))
lon = TS_era&longitude
lat = TS_era&latitude
;TS_era = linint2_Wrap(lon, lat, TS_ERA, True, LON, LAT, 0)
;printVarSummary(TS_era)
ts_era = clmMonTLL(TS_era)

;======================================================================
diriobs = "/home/yangsong3/wuyt/sysu/scripts/TP_coldbias/"
filiobs = "station2grid_TP_TAS.nc"
fobs = addfile(diriobs+filiobs, "r")
TS_obs = fobs->TAS
TS_obs = TS_obs+273.15
ts_obs = dim_avg_n_Wrap(TS_obs, 0) ; year,month,lat,lon

;======================================================================
dirictrl =  "/home/yangsong3/data-model/wuyt/TPbias_Ctrl_FAMIP/"
filictrl = "Ctrl.esm.cam.h0.TREFHT.1979-2005.nc"

fctrl = addfile(dirictrl+filictrl, "r")
TS_CRTL = fctrl->TREFHT
LON = TS_CRTL&lon
LAT = TS_CRTL&lat
TS_ctrl = linint2_Wrap(LON, LAT, TS_CRTL, True, lon, lat, 0)
ts_ctrl = clmMonTLL(TS_ctrl)
;======================================================================

differa = ts_ctrl
differa = ts_ctrl - ts_era
printVarSummary(differa)

diffobs = ts_ctrl
diffobs = ts_ctrl - ts_obs
printVarSummary(diffobs)


ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)
;topo = linint2_Wrap(TOPO&lon, TOPO&lat, TOPO, True, LON, LAT, 0)
topo = smth9_Wrap(topo, 0.25, 0.25, True)

differa = mask(differa, topo.ge.1500, True)
diffobs = mask(diffobs, topo.ge.1500, True)
;======================================================================

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/"
figname = "12mn_TP_TS_CTRL_minus_OBS"

lonstr = 60
lonend = 120
latstr = 20
latend = 50

wks = gsn_open_wks("pdf", figpath+figname)
plot = new(12, graphic)

res = True
res@gsnMaximize = True
res@gsnDraw  = False
res@gsnFrame = False

res@mpFillOn = True
res@mpMinLonF = lonstr
res@mpMaxLonF = lonend
res@mpMinLatF = latstr
res@mpMaxLatF = latend

res@cnFillOn = True
res@cnFillPalette = "hotcold_18lev"
res@lbLabelBarOn = False
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False


res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
; res@cnMinLevelValF       = -0.1               ; set min contour level
; res@cnMaxLevelValF       = 0.1                 ; set max contour level
; res@cnLevelSpacingF      = 0.01              ; set contour spacing
res@cnMinLevelValF       = -8               ; set min contour level
res@cnMaxLevelValF       = 8                 ; set max contour level
res@cnLevelSpacingF      = 0.5              ; set contour spacing
res@gsnRightString       = "K"
res@gsnStringFontHeightF = 0.035
res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025

Mth = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
do i = 0, 11
    res@gsnLeftString = Mth(i)
    plot(i) = gsn_csm_contour_map_ce(wks, differa(i,:,:), res)
end do

pres = True
pres@gsnDraw = True
pres@gsnFrame = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = "Tas diff: Ctrl - ERA"
pres@pmLabelBarOrthogonalPosF = -0.02
pres@lbLabelFontHeightF = 0.015

gsn_panel(wks, plot, (/3,4/), pres)

do i = 0, 11
    res@gsnLeftString = Mth(i)
    plot(i) = gsn_csm_contour_map_ce(wks, diffobs(i,:,:), res)
end do

pres = True
pres@gsnDraw = True
pres@gsnFrame = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = "Tas diff: Ctrl - OBS(station)"
pres@pmLabelBarOrthogonalPosF = -0.02
pres@lbLabelFontHeightF = 0.015

gsn_panel(wks, plot, (/3,4/), pres)

end