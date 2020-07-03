
begin

yrStrt = 1961
yrLast = 2018
year = ispan(yrStrt, yrLast, 1)
nyr = dimsizes(year)
month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
nm = dimsizes(month)

lonstr = 75
lonend = 105
latstr = 25
latend = 40

;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)
;======================================================================
;===================================================================
ts_obsx = new((/nyr*nm,72,128/), float,-99.)

diriobs = "/home/yangsong3/data-observation/stationChina2000/Grid/t2m/"

do i = 0, nyr-1
    yr = yrStrt+i
    do j = 0, 11
        mn = sprinti("%0.2i", j+1)
        filiobs = "SURF_CLI_CHN_TEM_MON_GRID_0.5-MEAN-"+yr+mn+".txt"
        x = asciiread(diriobs+filiobs, -1, "float")
        ts_obsx(i*12+j,:,:) = reshape(x(6:), (/72,128/))
    end do
end do
printVarSummary(ts_obsx)

ts_obsx = ts_obsx(:,::-1,:)
time = ispan(1, nyr*12, 1)
ts_obs_dtrend = dtrend_msg_n(time, ts_obsx, False, False, 0)
ts_obs = reshape(ts_obs_dtrend, (/nyr,nm,72,128/))

lat = fspan(18.25, 53.75, 72)
lon = fspan(72.25, 135.75, 128)

ts_obs!0 = "year"
ts_obs!1 = "month"
ts_obs!2 = "lat"
ts_obs!3 = "lon"
ts_obs&year = year
ts_obs&month = month
ts_obs&lat = lat
ts_obs&lon = lon
ts_obs&lat@units = "degrees_north"
ts_obs&lon@units = "degrees_east"
printVarSummary(ts_obs)

TS_obs = dim_avg_n_Wrap(ts_obs, 1)

topox1 = ftopo->topo({18:54},{72:136})
lonx = topox1&lon
latx = topox1&lat
printVarSummary(topox1)
topox2 = linint2_Wrap(lonx, latx, topox1, False, lon, lat, 0)
print(topox2({35},{90}))

TS_obs = mask(TS_obs, topox2.ge.1500, True)
TP_obs = wgt_areaave_Wrap(TS_obs(:,{latstr:latend},{lonstr:lonend}), 1, 1, 0)
printVarSummary(TP_obs)
;===================================================================

diriobs = "/home/yangsong3/data-observation/SST/NOAA_ERSST/"
filiobs = "sst.mnmean.v5.nc"
fobs = addfile(diriobs+filiobs, "r")
TIMEo = fobs->time
YYYYo = cd_calendar(TIMEo, -1)/100
iYYYYo = ind(YYYYo.ge.yrStrt .and. YYYYo.le.yrLast)

SST_OBS = fobs->sst(iYYYYo,::-1,:)

sst_OBS = dtrend_msg_n(SST_OBS&time, SST_OBS, False, False, 0)
sst_obs = month_to_annual(SST_OBS, 1)
printVarSummary(sst_obs)

; ;===================================================================
; ; land sea mask
; ;===================================================================
; a      = addfile ("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")
; lsdata = a->LSMASK
; lsm = landsea_mask(lsdata, SST_GCM&lat, SST_GCM&lon)
; ;lmask = where(lsm.eq.1 .or. lsm.eq.3 .or. lsm.eq.4, 0, 1)
; SST_GCM = mask(SST_GCM, lsm.eq.1 .or. lsm.eq.3 .or. lsm.eq.4, False)
;===================================================================
; regression data
;===================================================================
cc = sst_obs(0,:,:)
prob = cc
printVarSummary(TP_obs)
printVarSummary(sst_obs)
cc = escorc_n(TP_obs, sst_obs, 0, 0)
prob = rtest(cc, nyr, 0)
prob = mask(prob, cc.eq.0, False)
printVarSummary(cc)

;===================================================================
; Plot
;===================================================================
filipath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/a20191206/"
filiname = "Ann_TP_ts_corr_SST_OBS"
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
res@mpFillDrawOrder = "PostDraw"
res@mpCenterLonF = 180

res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False

res@gsnStringFontHeightF = 0.035
res@tmXBLabelFontHeightF = 0.02
res@tmYLLabelFontHeightF = 0.02

;res@lbLabelBarOn = False
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 0.1
res@cnMaxLevelValF = 0.5
res@cnMinLevelValF = -0.5
res@cnFillPalette = "CBR_coldhot"
;res@gsnRightString = "K"

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

sig = 0.95
res@gsnRightString = sig*100+"%"

res@gsnLeftString = "OBS"
plot(0) = gsn_csm_contour_map(wks, cc, res)
plot_prob(0) = gsn_csm_contour(wks, prob, res3)
plot_prob(0) = ShadeLtContour(plot_prob(0),(1-sig),6)
overlay(plot(0), plot_prob(0))

draw(plot(0))
frame(wks)
end