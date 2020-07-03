
begin

yrStrt = 1979
yrLast = 2017

lonstr = 75
lonend = 105
latstr = 25
latend = 40

;===================================================================
; Reanalysis ERAI data
;===================================================================
diriera = "/home/yangsong3/data-observation/wuyt/ERAI/monthly/1x1/"
;filir = "t2m.mon.mean.7901.1804.nc"
filiera = "skt.mon.mean.7917.nc"
fera = addfile(diriera+filiera, "r")
; timer = fera->time
; TIME = cd_calendar(timer, -1)
; print(TIME)
;TS_era = short2flt(fr->ts(0:467,:,:))
TS_era = short2flt(fera->skt(:,::-1,:))
printVarSummary(TS_era)
;print(TS_era&latitude)
; ts_era = clmMonTLL(TS_era)
; printVarSummary(ts_era)
TP_era = wgt_areaave(TS_era(:,{latstr:latend},{lonstr:lonend}), 1, 1, 0) 
printVarSummary(TP_era)


;===================================================================
; ERSST v5 data
;===================================================================
dirisst = "/home/yangsong3/data-observation/SST/NOAA_ERSST/"
filisst = "sst.mnmean.v5.nc"
fsst = addfile(dirisst+filisst, "r")

TIME = fsst->time
YYYY = cd_calendar(TIME, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

Time = fsst->time(iYYYY)
time = ut_calendar(Time, -1)
lat = fsst->lat(::-1)
lon = fsst->lon

sst = fsst->sst(iYYYY,::-1,:)
sst&time = time
printVarSummary(sst)

;===================================================================
; regression data
;===================================================================

rc = regCoef_n(TP_era, sst, 0, 0)
copy_VarCoords(sst(0,:,:), rc)
printVarSummary(rc)

cc = escorc_n(TP_era, sst, 0, 0)
copy_VarCoords(sst(0,:,:), cc)
printVarSummary(cc)


; mxlag = 6
; Sst = sst(lat|:, lon|:, time|:)
; ccx = esccr(Sst, TP_era, mxlag)
; cc = ccx(:,:,0:mxlag:-1)
; copy_VarCoords(sst(0,:,:), cc(:,:,0))
; printVarSummary(cc)

;===================================================================
; Plot
;===================================================================
filipath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/"
filiname = "TP_ts_reg_on_SST_OBS"
wks = gsn_open_wks("pdf", filipath+filiname)
plot = new(7, graphic)

res = True
res@gsnDraw  = False
res@gsnFrame = False
res@gsnMaximize = True

res@mpOutlineOn = True
res@mpFillOn = True
res@mpGeophysicalLineThicknessF = 2
res@mpGeophysicalLineColor = "grey25" ;"red" 
res@mpFillDrawOrder = "PreDraw"
res@mpCenterLonF = 180

res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False

res@gsnStringFontHeightF = 0.025
res@tmXBLabelFontHeightF = 0.02
res@tmYLLabelFontHeightF = 0.02

;res@lbLabelBarOn = False
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 0.1
res@cnMaxLevelValF = 0.5
res@cnMinLevelValF = -0.5
res@cnFillPalette = "CBR_coldhot"
;res@gsnRightString = "K"
;res@cnLevelSelectionMode = "ExplicitLevels"

res@gsnLeftString = "regression coefficient"
plot(0) = gsn_csm_contour_map(wks, rc, res)
; draw(plot)
; frame(wks)
res@gsnLeftString = "correlation coefficient"
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 0.2
res@cnMaxLevelValF = 1.
res@cnMinLevelValF = -1.
plot(1) = gsn_csm_contour_map(wks, cc, res)

; do i = 0, 6
;     res@gsnLeftString = "Lag = "+(-1*i)
;     plot(i) = gsn_csm_contour_map(wks, cc(:,:,i), res)
; end do

pres = True
pres@gsnMaximize = True
;pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = "ERAI Ts & ERSST.v5"
pres@gsnPanelFigureStringsFontHeightF = 0.035

gsn_panel(wks, plot, (/2,1/), pres)


end