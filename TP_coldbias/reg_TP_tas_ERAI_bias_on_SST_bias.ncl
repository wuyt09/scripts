begin

lonstr = 75
lonend = 105
latstr = 25
latend = 40

;===================================================================
; GCM 24 data
;===================================================================
dirigcm = "/home/yangsong3/data-model/wuyt/CMIP5/aaa-post/tas/"

models = (/"ACCESS1-0","ACCESS1-3","bcc-csm1-1","CanESM2","CCSM4","CESM1-BGC","CESM1-CAM5",\
        "CMCC-CM","CMCC-CMS","CNRM-CM5","CSIRO-Mk3-6-0","GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M",\
        "GISS-E2-H","GISS-E2-R","HadGEM2-CC","HadGEM2-ES","inmcm4","IPSL-CM5A-LR","IPSL-CM5A-MR",\
        "IPSL-CM5B-LR","MIROC5","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-M"/) ;,"NorESM1-ME"

; models = (/"bcc-csm1-1","CanESM2","CCSM4","CESM1-BGC","CESM1-CAM5",\
;         "CMCC-CM","CMCC-CMS","CNRM-CM5","CSIRO-Mk3-6-0","GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M",\
;         "GISS-E2-H","GISS-E2-R","inmcm4","IPSL-CM5A-LR","IPSL-CM5A-MR",\
;         "IPSL-CM5B-LR","MIROC5","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-M","NorESM1-ME"/)

nm = dimsizes(models)
TS_GCM = new((/nm,12,181,360/), float)
TS_GCM!0 = "model"
TS_GCM&model = models
do i = 0, nm-1
    filigcm = systemfunc("ls "+dirigcm+"tas_mon_"+models(i)+"_historical_r1i1p1_195101-200012.nc")
    fgcm = addfile(filigcm, "r")
    TS_GCMs = fgcm->tas
    TS_gcms = dtrend_msg_n(TS_GCMs&time, TS_GCMs, False, False, 0)
    copy_VarCoords(TS_GCMs, TS_gcms)
    TS_GCM(i,:,:,:) = clmMonTLL(TS_gcms)
end do
printVarSummary(TS_GCM)
TP_GCM = wgt_areaave(TS_GCM(:,:,{latstr:latend},{lonstr:lonend}), 1, 1, 0)
printVarSummary(TP_GCM)
; TP_GCMA = dim_avg_n_Wrap(TP_GCM, 0)
; printVarSummary(TP_GCMA)

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
;TS_era = short2flt(fera->skt(:,::-1,:))
printVarSummary(TS_era)
TS_era_dt = dtrend_msg_n(TS_era&time, TS_era, False, False, 0)
copy_VarCoords(TS_era, TS_era_dt)
;print(TS_era&latitude)
ts_era = clmMonTLL(TS_era_dt)
printVarSummary(ts_era)
TP_era = wgt_areaave(ts_era(:,{latstr:latend},{lonstr:lonend}), 1, 1, 0) 
printVarSummary(TP_era)

ts_bias = TP_GCM
do i = 0, nm-1
    ts_bias(i,:) = TP_GCM(i,:) - TP_era
end do
printVarSummary(ts_bias)

;===================================================================
; SST historical bias data
;===================================================================
yrStrtmodel1 = 1950
yrLastmodel1 = 1999

dirh="ls /home/yangsong3/jyt/ensocro/sst/ssth/"
;read the name of all of models
dirhname = asciiread("/home/yangsong3/jyt/name/dirssth28.txt",(/28/),"string")

ssth = new((/600,91,180/), float)
ssthclm = new((/nm,12,91,180/), float)
time = ispan(1, 600, 1)

do i = 0, nm-1
    ;==============read SST==============================
    filh = systemfunc(dirh+dirhname(i))
    ssth = fbindirread(filh, 0, (/600,91,180/), "float")
    ssth!0 = "time"
    ssth!1 = "lat"
    ssth!2 = "lon"
    ssth&time = time
    ssth&lat  = fspan(-90, 90, 91)
    ssth&lon  = fspan(0, 358, 180)
    ssth&lat@units = "degrees_north"
    ssth&lon@units = "degrees_east"

    if (any(isnan_ieee(ssth))) then
        value = 1.e+20
        replace_ieeenan(ssth, value, 0)
        ssth@_FillValue = value
    end if

    ;------detrend---------
    ssthdt  =  dtrend_msg_n(time,ssth,False,False,0);去趋势
    copy_VarAtts(ssth,ssthdt)
    copy_VarCoords(ssth,ssthdt)

    ;------anomaly---------
    ssthclm(i,:,:,:) = clmMonTLL(ssthdt)
end do
printVarSummary(ssthclm)
;===============================================
;bias of models
;===============================================
;----------------read reanalysis----------------
  dir1 = "/home/yangsong3/data-observation/SST/NOAA_ERSST/"
  fName = addfile(dir1+"sst.mnmean.v5.nc","r")

  TIME      = fName->time
  YYYY      = cd_calendar(TIME,-1)/100
  iYYYY     = ind(YYYY.ge.yrStrtmodel1 .and. YYYY.le.yrLastmodel1)
  Time      = fName->time(iYYYY)
  time      = ut_calendar(Time,-1)
  lon       = fName->lon
  lat       = fName->lat(::-1)

  sstobs    = fName->sst(iYYYY,::-1,:)
  sstobs!0  = "time"
  sstobs!1  = "lat"
  sstobs!2  = "lon"
  sstobs&time = time

 ;------interpolate sstobs to same resolution to SST from CMIP5---------
  LON = fspan(0,358,180)
  LAT = fspan(-90,90,91)
  copy_VarAtts(lon, LON)
  copy_VarAtts(lat, LAT)
  SSTOBS = linint2_Wrap(lon,lat,sstobs,True,LON,LAT,0)

  ;------detrend---------
  sstodt  =  dtrend_msg_n(time,SSTOBS,False,False,0);去趋势
  copy_VarAtts(SSTOBS,sstodt)
  copy_VarCoords(SSTOBS,sstodt)

  ;------anomaly---------
  sstoclm = clmMonTLL(sstodt)

  ;-------model bias------
  sstOBSclm = conform(ssthclm, sstoclm, (/1,2,3/))
  biasclm = ssthclm-sstOBSclm-273.15
  copy_VarMeta(ssthclm, biasclm)
  biasclm!0 = "time"
  printVarSummary(biasclm)

;===================================================================
; reg TP ts bias on historical SST bias
;===================================================================
rc = new((/12, 91, 180/), float)
rc = rc@_FillValue
cc = rc
do i = 0, 11
    rc(i,:,:) = regCoef_n(ts_bias(:,i), biasclm(:,i,:,:), 0, 0)
    cc(i,:,:) = escorc_n(ts_bias(:,i), biasclm(:,i,:,:), 0, 0)
end do
copy_VarCoords(biasclm(0,:,:,:), cc)
copy_VarCoords(biasclm(0,:,:,:), rc)
printVarSummary(rc)

ts_bias_ann = dim_avg_n_Wrap(ts_bias, 1)
biasclm_ann = dim_avg_n_Wrap(biasclm, 1)
rc_ann = regCoef_n(ts_bias_ann, biasclm_ann, 0, 0)
cc_ann = escorc_n(ts_bias_ann, biasclm_ann, 0, 0)
copy_VarCoords(biasclm_ann(0,:,:), cc_ann)
copy_VarCoords(biasclm_ann(0,:,:), rc_ann)
printVarSummary(rc_ann)

;===================================================================
; Plot
;===================================================================
filipath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/"
filiname = "reg_TP_tas_ERAI_bias_on_SST_bias"
wks = gsn_open_wks("pdf", filipath+filiname)
plot = new(12, graphic)

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
; res@mpMinLatF = -30
; res@mpMaxLatF = 30

res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False

res@gsnStringFontHeightF = 0.025
res@tmXBLabelFontHeightF = 0.02
res@tmYLLabelFontHeightF = 0.02

res@lbLabelBarOn = False
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 0.05
res@cnMaxLevelValF = 0.5
res@cnMinLevelValF = -0.5
res@cnFillPalette = "CBR_coldhot"

do i = 0, 11
    res@gsnLeftString = "Month = "+(i+1)
    plot(i) = gsn_csm_contour_map(wks, rc(i,:,:), res)
end do

pres = True
pres@gsnPanelMainString = "regression coefficient of TP tas bias on historical SST"
pres@pmLabelBarOrthogonalPosF = -0.02
pres@gsnPanelLabelBar = True

gsn_panel(wks, plot(0::2), (/3,2/), pres)
gsn_panel(wks, plot(1::2), (/3,2/), pres)

res@lbLabelBarOn = True
res@gsnLeftString = "Annual mean"
res@pmLabelBarOrthogonalPosF = 0.4
plot(0) = gsn_csm_contour_map(wks, rc_ann, res)
draw(plot(0))
frame(wks)

end