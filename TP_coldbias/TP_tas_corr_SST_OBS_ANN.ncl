
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
filiera = "t2m.mon.mean.7901.1804.nc"
;filiera = "skt.mon.mean.7917.nc"
fera = addfile(diriera+filiera, "r")
timer = fera->time
YYYYR = cd_calendar(timer, -1)/100
iYYYYR = ind(YYYYR.ge.1979 .and. YYYYR.le.2017)
TS_era = short2flt(fera->t2m(iYYYYR,:,:))
printVarSummary(TS_era)
TS_era = dtrend_msg_n(TS_era&time, TS_era, False, False, 0)
ts_era = month_to_annual(TS_era, 1)
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
SST = month_to_annual(sst, 1)
SST&year = ispan(yrStrt, yrLast, 1)
printVarSummary(SST)

;===================================================================
; regression data
;===================================================================
Nr = dimsizes(SST&year)
cc = escorc_n(TP_era, SST, 0, 0)
prob = rtest(cc, Nr, 0)
prob = mask(prob, cc.eq.0, False)
copy_VarCoords(sst(0,:,:), cc)
copy_VarCoords(sst(0,:,:), prob)
printVarSummary(cc)

;===================================================================
; Plot
;===================================================================
filipath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/"
filiname = "TP_tas_corr_SST_OBS_ANN"
wks = gsn_open_wks("pdf", filipath+filiname)
plot = new(7, graphic)
plot_prob = plot

res = True
res@gsnDraw  = False
res@gsnFrame = False
res@gsnMaximize = True

res@mpOutlineOn = True
res@mpFillOn = True
res@mpGeophysicalLineThicknessF = 2
res@mpGeophysicalLineColor = "grey25" ;"red" 
res@mpFillDrawOrder = "PostDraw"
res@mpCenterLonF = 180

res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False

res@gsnStringFontHeightF = 0.025
res@tmXBLabelFontHeightF = 0.02
res@tmYLLabelFontHeightF = 0.02

res@cnFillPalette = "CBR_coldhot"
res@gsnLeftString = "Annual mean"
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 0.1
res@cnMaxLevelValF = 0.6
res@cnMinLevelValF = -0.6
res@tiMainString = "Correlation coefficient: ERAI TP Tas & ERSST.v5 "

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
res3@cnFillColor  =  "green4" ;"black"
res3@gsnLeftString = ""
res3@cnFillDotSizeF = 0.003
res3@cnFillScaleF  = 1 ;1.2

sig = 0.95
res@gsnRightString = sig*100+"%"

plot(0) = gsn_csm_contour_map(wks, cc, res)
plot_prob(0) = gsn_csm_contour(wks, prob, res3)
plot_prob(0) = ShadeLtContour(plot_prob(0),(1-sig),6)
overlay(plot(0), plot_prob(0))
draw(plot)
frame(wks)

end