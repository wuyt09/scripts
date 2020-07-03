
begin

yrStrt = 1961
yrLast = 2000
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
;===================================================================
; station obs data
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
; GCM 24 data
;===================================================================
dirigcm = "/home/yangsong3/data-model/wuyt/CMIP5/aaa-post/tas/"

models = (/"ACCESS1-0","ACCESS1-3","bcc-csm1-1","CanESM2","CCSM4","CESM1-BGC","CESM1-CAM5",\
        "CMCC-CM","CMCC-CMS","CNRM-CM5","CSIRO-Mk3-6-0","GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M",\
        "GISS-E2-H","GISS-E2-R","HadGEM2-CC","HadGEM2-ES","inmcm4","IPSL-CM5A-LR","IPSL-CM5A-MR",\
        "IPSL-CM5B-LR","MIROC5","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-M","NorESM1-ME"/)

; models = (/"bcc-csm1-1","CanESM2","CCSM4","CESM1-BGC","CESM1-CAM5",\
;         "CMCC-CM","CMCC-CMS","CNRM-CM5","CSIRO-Mk3-6-0","GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M",\
;         "GISS-E2-H","GISS-E2-R","inmcm4","IPSL-CM5A-LR","IPSL-CM5A-MR",\
;         "IPSL-CM5B-LR","MIROC5","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-M","NorESM1-ME"/)

nm = dimsizes(models)
TS_GCM = new((/nm,nyr,181,360/), float)
TS_GCM!0 = "model"
TS_GCM&model = models
do i = 0, nm-1
    filigcm = systemfunc("ls "+dirigcm+"tas_mon_"+models(i)+"_historical_r1i1p1_195101-200012.nc")
    fgcm = addfile(filigcm, "r")
    TIMEm = fgcm->time
    YYYYm = cd_calendar(TIMEm, -1)/100
    iYYYYm = ind(YYYYm.ge.yrStrt .and. YYYYm.le.yrLast)
    TS_GCMs = fgcm->tas(iYYYYm,:,:)
    TS_GCMs = dtrend_msg_n(TS_GCMs&time, TS_GCMs, False, False, 0)
    TS_GCM(i,:,:,:) = month_to_annual(TS_GCMs, 1)
end do
printVarSummary(TS_GCM)
;TP_GCMX = TS_GCM
TS_GCM = mask(TS_GCM, topo.ge.1500, True)
TP_GCM = wgt_areaave_Wrap(TS_GCM(:,:,{latstr:latend},{lonstr:lonend}), 1, 1, 0)
printVarSummary(TP_GCM)

ts_bias = TP_GCM 
ts_bias = TP_GCM - conform_dims(dimsizes(TP_GCM), TP_obs, 1)
;===================================================================
; ersst obs sst
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
;===================================================================
dirigcm = "/home/yangsong3/data-model/wuyt/CMIP5/aaa-post/ts/"
SST_GCM = new((/nm,nyr,181,360/), float)
SST_GCM!0 = "model"
SST_GCM&model = models
do i = 0, nm-1
    filigcm = systemfunc("ls "+dirigcm+"ts_mon_"+models(i)+"_historical_r1i1p1_195101-200012.nc")
    fgcm = addfile(filigcm, "r")
    TIMEm = fgcm->time
    YYYYm = cd_calendar(TIMEm, -1)/100
    iYYYYm = ind(YYYYm.ge.yrStrt .and. YYYYm.le.yrLast)
    SST_GCMs = fgcm->ts(iYYYYm,:,:)
    SST_GCMs = dtrend_msg_n(SST_GCMs&time, SST_GCMs, False, False, 0)
    SST_GCM(i,:,:,:) = month_to_annual(SST_GCMs, 1)
end do
printVarSummary(SST_GCM)

sst_GCM = linint2_Wrap(SST_GCM&lon, SST_GCM&lat, SST_GCM, True, sst_obs&lon, sst_obs&lat, 0)

sst_bias = sst_GCM
sst_bias = sst_GCM - conform_dims(dimsizes(sst_GCM), sst_obs, (/1,2,3/))
;===================================================================
; land sea mask
;===================================================================
a      = addfile ("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")
lsdata = a->LSMASK
lsm = landsea_mask(lsdata, sst_bias&lat, sst_bias&lon)
;lmask = where(lsm.eq.1 .or. lsm.eq.3 .or. lsm.eq.4, 0, 1)
sst_bias = mask(sst_bias, lsm.eq.1 .or. lsm.eq.3 .or. lsm.eq.4, False)

;===================================================================
; regression data
;===================================================================
cc = new((/nm,89,180/), float)
prob = cc
do i = 0, nm-1
    cc(i,:,:) = escorc_n(ts_bias(i,:), sst_bias(i,:,:,:), 0, 0)
    prob(i,:,:) = rtest(cc(i,:,:), nm, 0)
end do
copy_VarCoords(sst_bias(:,0,:,:), cc)
copy_VarCoords(sst_bias(:,0,:,:), prob)
prob = mask(prob, cc.eq.0, False)
printVarSummary(cc)

;===================================================================
; Plot
;===================================================================
filipath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/a20191206/"
filiname = "Ann_TP_ts_bias_corr_SST_bias_GCMs"
wks = gsn_open_wks("pdf", filipath+filiname)
plot = new(nm, graphic)
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

res@lbLabelBarOn = False
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 0.1
res@cnMaxLevelValF = 0.8
res@cnMinLevelValF = -0.8
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

do i = 0, nm-1
    res@gsnLeftString = models(i)
    plot(i) = gsn_csm_contour_map(wks, cc(i,:,:), res)
    plot_prob(i) = gsn_csm_contour(wks, prob(i,:,:), res3)
    plot_prob(i) = ShadeLtContour(plot_prob(i),(1-sig),6)
    overlay(plot(i), plot_prob(i))
end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = "Correlation coefficient of annual mean: model Ts bias & model SST bias "
gsn_panel(wks, plot, (/7, 4/), pres)

cc_mme = dim_avg_n_Wrap(cc, 0)
prob_mme = dim_avg_n_Wrap(prob, 0)
res@gsnLeftString = "MME"
res@lbLabelBarOn = True
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 0.05
res@cnMaxLevelValF = 0.3
res@cnMinLevelValF = -0.3
plot(0) = gsn_csm_contour_map(wks, cc_mme, res)
plot_prob(0) = gsn_csm_contour(wks, prob_mme, res3)
plot_prob(0) = ShadeLtContour(plot_prob(0),(1-sig),6)
overlay(plot(0), plot_prob(0))
draw(plot(0))
frame(wks)
end