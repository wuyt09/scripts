
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
TS_era = dtrend_msg_n(TS_era&time, TS_era, False, False, 0)
printVarSummary(TS_era)
;print(TS_era&latitude)
TS_era_SON = month_to_season(TS_era, "SON")
TS_era_DJF = month_to_season(TS_era, "DJF")
ts_era = (TS_era_SON + TS_era_DJF)/2
copy_VarCoords(TS_era(0,:,:), ts_era(0,:,:))
printVarSummary(ts_era)

TP_era = wgt_areaave(ts_era(:,{latstr:latend},{lonstr:lonend}), 1, 1, 0) 
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
sst = dtrend_msg_n(sst&time, sst, False, False, 0)
printVarSummary(sst)

sst_SON = month_to_season(sst, "SON")
sst_DJF = month_to_season(sst, "DJF")
SST = sst_SON
SST = (sst_SON + sst_DJF)/2
printVarSummary(SST)

;===================================================================
; regression data
;===================================================================

; rc = regCoef_n(TP_era, SST, 0, 0)
; copy_VarCoords(sst(0,:,:), rc)
; printVarSummary(rc)

cc = escorc_n(TP_era, SST, 0, 0)
copy_VarCoords(sst(0,:,:), cc)
printVarSummary(cc)
Nr = dimsizes(SST&time)
prob = cc
prob = rtest(cc, Nr, 0)
printVarSummary(prob)
prob = mask(prob, cc.eq.0, False)

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
filiname = "TP_ts_corr_SST_OBS_SONDJF"
wks = gsn_open_wks("pdf", filipath+filiname)
plot = new(1, graphic)
plot_prob = plot

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

res@cnFillPalette = "CBR_coldhot"
res@gsnLeftString = "SONDJF"
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 0.1
res@cnMaxLevelValF = 0.6
res@cnMinLevelValF = -0.6
res@tiMainString = "Correlation coefficient: ERAI TP Ts & ERSST.v5 "

sig = 0.95
res@gsnRightString = sig*100+"%"

plot(0) = gsn_csm_contour_map(wks, cc, res)

;------------------------
; draw prob
;------------------------
res3  = True
res3@gsnDraw   = False
res3@gsnFrame  = False

;res3@cnFillOn  = False 
res3@cnLinesOn = False
res3@cnLineLabelsOn = False
res3@cnInfoLabelOn  = False
;res3@lbLabelBarOn   = False

res3@gsnRightString = ""

res3@cnLevelSelectionMode = "ManualLevels"
res3@cnMinLevelValF = 0.00
res3@cnMaxLevelValF = 0.10
res3@cnLevelSpacingF = 0.01

;res3@cnFillPattern = 17 ;(/17/)
res3@cnFillColor  =  "black" ;"green4"
res3@gsnLeftString = ""
res3@cnFillDotSizeF = 0.003
res3@cnFillScaleF  = 1 ;1.2

plot_prob(0) = gsn_csm_contour(wks, prob, res3)
plot_prob(0) = ShadeLtContour(plot_prob(0),(1-sig),6)
overlay(plot(0), plot_prob(0))
draw(plot)
frame(wks)

end