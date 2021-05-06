
begin
lonstr = 75 ; 72 ;75
lonend = 105 ;136 ;105
latstr = 25 ;18 ;25
latend = 40 ;54 ; 40

season = "Ann"
;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)
printVarSummary(topo)

;===================================================================
; OBS grid data
;===================================================================
yrStrt = 1961
yrLast = 2018
year = ispan(yrStrt, yrLast, 1)
nyr = dimsizes(year)
month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
nm = dimsizes(month)

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
ts_obs = dtrend_msg_n(time, ts_obsx, False, False, 0)
; ts_obs = reshape(ts_obs_dtrend, (/nyr,nm,72,128/))

lat = fspan(18.25, 53.75, 72)
lon = fspan(72.25, 135.75, 128)
print(lat)
print(lon)

ts_obs!0 = "time"
ts_obs!1 = "lat"
ts_obs!2 = "lon"
ts_obs&time = ispan(1, nyr*nm, 1)
ts_obs&lat = lat
ts_obs&lon = lon
ts_obs&lat@units = "degrees_north"
ts_obs&lon@units = "degrees_east"
printVarSummary(ts_obs)

TS_obs = month_to_annual(ts_obs, 1)
; TS_obs = dim_avg_n_Wrap(TS_obsx, 0)
printVarSummary(TS_obs)

topox1 = ftopo->topo({18:54},{72:136})
lonx = topox1&lon
latx = topox1&lat
printVarSummary(topox1)
topox2 = linint2_Wrap(lonx, latx, topox1, False, lon, lat, 0)
print(topox2({35},{90}))

TS_obs = mask(TS_obs, topox2.ge.1500, True)
TP_obs = wgt_areaave_Wrap(TS_obs(:,{latstr:latend},{lonstr:lonend}), 1, 1, 0)
printVarSummary(TS_obs)

;===================================================================
; obs sst
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
sst_med = dim_avg_n_Wrap(sst_obs, 0)
sst_obs_multiyear = linint2_Wrap(sst_med&lon, sst_med&lat, sst_med, True, fspan(0, 359, 360), fspan(-90, 90, 181), 0)
;===================================================================
; regression data
;===================================================================
; cc_obs = sst_obs(0,:,:)
; prob_obs = cc_obs
; printVarSummary(TP_obs)
; printVarSummary(sst_obs)
; cc_obs = escorc_n(TP_obs, sst_obs, 0, 0)
; prob_obs = rtest(cc_obs, nyr, 0)
; prob_obs = mask(prob_obs, cc_obs.eq.0, False)
; printVarSummary(cc_obs)
;===================================================================
; GCM 28 tas data, calc the cold bias
;===================================================================
dirigcm = "/home/yangsong3/data-model/wuyt/CMIP5/aaa-post/tas/"

models = (/"ACCESS1-0","ACCESS1-3","bcc-csm1-1","CanESM2","CCSM4","CESM1-BGC","CESM1-CAM5",\
        "CMCC-CM","CMCC-CMS","CNRM-CM5","CSIRO-Mk3-6-0","GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M",\
        "GISS-E2-H","GISS-E2-R","HadGEM2-CC","HadGEM2-ES","inmcm4","IPSL-CM5A-LR","IPSL-CM5A-MR",\
        "IPSL-CM5B-LR","MIROC5","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-M","NorESM1-ME"/)

nm = dimsizes(models)
TS_GCM = new((/nm,50,181,360/), float)
TS_GCM!0 = "model"
TS_GCM&model = models
do i = 0, nm-1
    filigcm = systemfunc("ls "+dirigcm+"tas_mon_"+models(i)+"_historical_r1i1p1_195101-200012.nc")
    fgcm = addfile(filigcm, "r")
    TS_GCMs = fgcm->tas
    TS_GCMs = dtrend_msg_n(TS_GCMs&time, TS_GCMs, False, False, 0)
    TS_GCM(i,:,:,:) = month_to_annual(TS_GCMs, 1)
    ; TS_GCMx(i,:,:,:) = dim_avg_n_Wrap(TS_GCMx, 0)
end do

; TS_GCM = dim_avg_n_Wrap(TS_GCMx, 1)
printVarSummary(TS_GCM)
TS_GCM = TS_GCM - 273.15
TS_GCM = mask(TS_GCM, topo.ge.1500, True)
TP_GCM = wgt_areaave_Wrap(TS_GCM(:,:,{latstr:latend},{lonstr:lonend}), 1, 1, 0)
; TP_GCM_mme = dim_avg_n_Wrap(TP_GCM, 0)
; printVarSummary(TP_GCM_mme)

; TP_GCM_tm = dim_avg_n_Wrap(TP_GCM, 1)
; printVarSummary(TP_GCM_tm)
TP_GCM_multiyear = dim_avg_n_Wrap(TP_GCM, 1)
TP_obs_multiyear = dim_avg_n_Wrap(TP_obs, 0)
TP_bias = TP_GCM_multiyear - TP_obs_multiyear
TP_bias!0 = "models"
TP_bias&models = models
printVarSummary(TP_bias)
; ;==================================================================
; ; sort by bias
; ;==================================================================
; t_bias = TP_GCM_tm - TP_obs
; copy_VarCoords(TP_GCM_tm, t_bias)
; printVarSummary(t_bias)
; abs_tbias = abs(t_bias)
; ip_tbias = dim_pqsort(abs_tbias,1)
; print(ip_tbias)
; print(models(ip_tbias)+"  "+sprintf("%4.2f",t_bias(ip_tbias)))
; write_table("./bias_list_"+season+".txt", "w", [/ip_tbias/], "%2i")
; exit
;===================================================================
; model sst
;===================================================================
dirigcm = "/home/yangsong3/data-model/wuyt/CMIP5/aaa-post/ts/"
SST_GCM = new((/nm,181,360/), float)
SST_GCM!0 = "model"
SST_GCM&model = models
do i = 0, nm-1
    filigcm = systemfunc("ls "+dirigcm+"ts_mon_"+models(i)+"_historical_r1i1p1_195101-200012.nc")
    fgcm = addfile(filigcm, "r")
    SST_GCMs = fgcm->ts
    SST_GCMs = dtrend_msg_n(SST_GCMs&time, SST_GCMs, False, False, 0)
    x_med = month_to_annual(SST_GCMs, 1)
    x_med_multiyear = dim_avg_n_Wrap(x_med, 0)
    SST_GCM(i,:,:) = x_med_multiyear - sst_obs_multiyear
end do
copy_VarCoords(x_med_multiyear, SST_GCM(0,:,:))
printVarSummary(SST_GCM)
;===================================================================
; land sea mask
;===================================================================
a      = addfile ("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")
lsdata = a->LSMASK
lsm = landsea_mask(lsdata, SST_GCM&lat, SST_GCM&lon)
;lmask = where(lsm.eq.1 .or. lsm.eq.3 .or. lsm.eq.4, 0, 1)
SST_GCM = mask(SST_GCM, lsm.eq.1 .or. lsm.eq.3 .or. lsm.eq.4, False)
;===================================================================
; regression data
;===================================================================
cc = new((/181,360/), float)
prob = cc
cc = escorc_n(TP_bias, SST_GCM, 0, 0)
prob = rtest(cc, nm, 0)
copy_VarCoords(SST_GCM(0,:,:), cc)
copy_VarCoords(cc, prob)
;===================================================================
; Plot
;===================================================================
filipath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/a20201012/"
filiname = "Corr_TP_bias_SST_bias"
wks = gsn_open_wks("pdf", filipath+filiname)
plot = new(2, graphic)
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

res@gsnStringFontHeightF = 0.03
res@tmXBLabelFontHeightF = 0.018
res@tmYLLabelFontHeightF = 0.018

res@lbLabelBarOn = False
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 0.05
res@cnMaxLevelValF = 0.5
res@cnMinLevelValF = -0.5
res@cnFillPalette = "BlueWhiteOrangeRed"
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

; do i = 0, nm-1
;     res@gsnLeftString = models(ip_tbias(i))
;     res@gsnCenterString = sprintf("%4.2f", t_bias(ip_tbias(i)))
;     plot(i) = gsn_csm_contour_map(wks, cc(ip_tbias(i),:,:), res)
;     plot_prob(i) = gsn_csm_contour(wks, prob(ip_tbias(i),:,:), res3)
;     plot_prob(i) = ShadeLtContour(plot_prob(i),(1-sig),6)
;     overlay(plot(i), plot_prob(i))
; end do

res@gsnLeftString = " "
plot(0) = gsn_csm_contour_map(wks, cc, res)
plot_prob(0) = gsn_csm_contour(wks, prob, res3)
plot_prob(0) = ShadeLtContour(plot_prob(0),(1-sig),6)
overlay(plot(0), plot_prob(0))
draw(plot(0))

; res@gsnLeftString = "MME"
; res@lbLabelBarOn = True
; res@cnLevelSpacingF = 0.05
; res@cnMaxLevelValF = 0.3
; res@cnMinLevelValF = -0.3
; plot(0) = gsn_csm_contour_map(wks, cc_mme, res)
; plot_prob(0) = gsn_csm_contour(wks, prob_mme, res3)
; plot_prob(0) = ShadeLtContour(plot_prob(0),(1-sig),6)
; overlay(plot(0), plot_prob(0))
; draw(plot(0))
; frame(wks)
end
