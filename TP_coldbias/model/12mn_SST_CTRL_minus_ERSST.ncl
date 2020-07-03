begin

yrStrt = 1950
yrLast = 1999
year = ispan(yrStrt, yrLast, 1)

;======================================================================
dirictrl = "/home/yangsong3/wuyt/sysu/scripts/TP_coldbias/data/"
filictrl = "sst_HadOIBl_bc_1.9x2.5_1850_2013_c140701.nc"
fctrl = addfile(dirictrl+filictrl, "r")
TIMEc = fctrl->time
YYYYc = cd_calendar(TIMEc, -1)/100
iYYYYc = ind(YYYYc.ge.yrStrt .and. YYYYc.le.yrLast)

sst_CTRL = fctrl->SST_cpl(iYYYYc,:,:)
LAT = sst_CTRL&lat
LON = sst_CTRL&lon
sst_ctrl = clmMonTLL(sst_CTRL)

;======================================================================
diriobs = "/home/yangsong3/data-observation/SST/NOAA_ERSST/"
filiobs = "sst.mnmean.v5.nc"
fobs = addfile(diriobs+filiobs, "r")
TIMEo = fobs->time
YYYYo = cd_calendar(TIMEo, -1)/100
iYYYYo = ind(YYYYo.ge.yrStrt .and. YYYYo.le.yrLast)

SST_OBS = fobs->sst(iYYYYo,::-1,:)
lat = SST_OBS&lat
lon = SST_OBS&lon

sst_OBS = linint2_Wrap(lon, lat, SST_OBS, True, LON, LAT, 0)
sst_obs = clmMonTLL(sst_OBS)

diff = sst_ctrl
diff = sst_ctrl - sst_obs
printMinMax(diff, 0)

;======================================================================

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/"
figname = "12mn_SST_CTRL_minus_ERSST"

lonstr = 0
lonend = 360
latstr = -30
latend = 30

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
res@cnMinLevelValF       = -1               ; set min contour level
res@cnMaxLevelValF       = 1                 ; set max contour level
res@cnLevelSpacingF      = 0.1              ; set contour spacing
res@gsnRightString       = "K"
res@gsnStringFontHeightF = 0.035
res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025

Mth = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
do i = 0, 11
    res@gsnLeftString = Mth(i)
    plot(i) = gsn_csm_contour_map_ce(wks, diff(i,:,:), res)
end do

pres = True
pres@gsnDraw = True
pres@gsnFrame = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = "SST diff: Ctrl - ERSST"
pres@pmLabelBarOrthogonalPosF = -0.02
pres@lbLabelFontHeightF = 0.015

gsn_panel(wks, plot, (/4,3/), pres)
end