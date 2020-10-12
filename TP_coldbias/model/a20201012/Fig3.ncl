;===========================================
; draw ensmble sst bias in CMIP5 models
; 12 month maps
;===========================================

begin

yrStrt = 1950
yrLast = 1999

;===============================================
; OBS 
;===============================================
;----------------read reanalysis----------------
  diri = "/home/yangsong3/data-observation/SST/NOAA_ERSST/"
  fName = addfile(diri+"sst.mnmean.v5.nc","r")

  TIME      = fName->time
  YYYY      = cd_calendar(TIME,-1)/100
  iYYYY     = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
  Time      = fName->time(iYYYY)
  time      = ut_calendar(Time,-1)
  lon       = fName->lon
  lat       = fName->lat(::-1)
  nlat      = dimsizes(lat)
  nlon      = dimsizes(lon)

  sstobs    = fName->sst(iYYYY,::-1,:)
  sstobs    = sstobs+273.15
  sstobs!0  = "time"
  sstobs!1  = "lat"
  sstobs!2  = "lon"
  sstobs&time = time
  sstobs@units = "K"
  printVarSummary(sstobs)


;===================================================================
; CMIP5 28 GCMs sst data
;===================================================================
dirigcm = "/home/yangsong3/data-model/wuyt/CMIP5/"

models = (/"ACCESS1-0","ACCESS1-3","bcc-csm1-1","CanESM2","CCSM4","CESM1-BGC","CESM1-CAM5",\
        "CMCC-CM","CMCC-CMS","CNRM-CM5","CSIRO-Mk3-6-0","GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M",\
        "GISS-E2-H","GISS-E2-R","HadGEM2-CC","HadGEM2-ES","inmcm4","IPSL-CM5A-LR","IPSL-CM5A-MR",\
        "IPSL-CM5B-LR","MIROC5","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-M","NorESM1-ME"/)

nm = dimsizes(models)
bias_sst5 = new((/nm,600,nlat,nlon/), float)
bias_sst5!0 = "model"
bias_sst5&model = models
do i = 0, nm-1
    filigcm = systemfunc("ls "+dirigcm+models(i)+"/ts_Amon_"+models(i)+"_historical_r1i1p1_18*-20*.nc")
    print(filigcm)
    fgcm = addfile(filigcm, "r")
    TIME := fgcm->time
    YYYY := cd_calendar(TIME, -1)/100
    iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
    TS_GCMs := fgcm->ts(iYYYY,:,:)
    TS_GCMs!1 = "lat"
    TS_GCMs!2 = "lon"

    TS_GCM := linint2_Wrap(TS_GCMs&lon, TS_GCMs&lat, TS_GCMs, True, lon, lat, 0)
    TS_GCM = mask(TS_GCM, ismissing(sstobs), False)
    bias_sst5(i,:,:,:) = TS_GCM - sstobs
    print("end of "+models(i))
end do
copy_VarCoords(sstobs, bias_sst5(0,:,:,:))
bias_sst5A = dim_avg_n_Wrap(bias_sst5, 0)
printVarSummary(bias_sst5A)

; ;===================================================================
; ; CMIP6 20 GCMs sst data
; ;===================================================================
; dirigcm = "/home/yangsong3/data-model/wuyt/CMIP6/"

; models := (/"BCC-CSM2-MR","BCC-ESM1","CanESM5","CAS-ESM2-0","CESM2","CESM2-FV2","CESM2-WACCM",\
;           "CESM2-WACCM-FV2","CIESM","E3SM-1-0","E3SM-1-1","E3SM-1-1-ECA","FGOALS-f3-L",\
;           "FGOALS-g3","KIOST-ESM","MCM-UA-1-0","MIROC6","MRI-ESM2-0","NESM3","SAM0-UNICON"/)
; ; models = (/"B" /)

; nm := dimsizes(models)
; bias_sst6 = new((/nm,600,nlat,nlon/), float)
; bias_sst6!0 = "model"
; bias_sst6&model = models
; do i = 0, nm-1
;     filigcm = systemfunc("ls "+dirigcm+"ts_Amon_"+models(i)+"_historical_r1i1p1f1_*_18*-20*.nc")
;     print(filigcm)
;     fgcm = addfile(filigcm, "r")
;     TIME := fgcm->time
;     YYYY := cd_calendar(TIME, -1)/100
;     iYYYY := ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
;     TS_GCMs := fgcm->ts(iYYYY,:,:)
;     TS_GCMs!1 = "lat"
;     TS_GCMs!2 = "lon"

;     TS_GCM := linint2_Wrap(TS_GCMs&lon, TS_GCMs&lat, TS_GCMs, True, lon, lat, 0)
;     printVarSummary(TS_GCM)
;     TS_GCM = mask(TS_GCM, ismissing(sstobs), False)
;     bias_sst6(i,:,:,:) = TS_GCM - sstobs
;     print("end of "+models(i))
; end do
; copy_VarCoords(sstobs, bias_sst6(0,:,:,:))
; bias_sst6A = dim_avg_n_Wrap(bias_sst6, 0)
; printVarSummary(bias_sst6A)


bias_clm = clmMonTLL(bias_sst5A)
printVarSummary(bias_clm)


Month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
Num = (/"(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)","(j)","(k)","(l)"/)
;=========================================================
; plot 
;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/a20201012/"
figname = "Fig3"

pwks = "png"
pwks@wkHeight = 3048
pwks@wkWidth  = 6048
wks = gsn_open_wks(pwks, figpath+figname)
plot = new(12, graphic)


res                       = True     ; plot mods desired
res@gsnDraw               = False
res@gsnFrame              = False

res@gsnMaximize           = True     ; Make plot fit the frame
res@cnFillOn              = True     ; turn on color fill
res@cnLinesOn             = False    ; turn of contour lines
;res@cnLevelSpacingF       = 0.5      ; contour spacing
res@cnFillPalette         = "hotcold_18lev"
res@lbLabelBarOn          = False
res@cnInfoLabelOn         = False
res@cnLineLabelsOn        = False

res@mpMinLatF            = -30      ; range to zoom in on
res@mpMaxLatF            =  30.
;res@mpMinLonF            =  30.
;res@mpMaxLonF            = 120.
res@mpCenterLonF         = 180.

res@tmXBLabelFontHeightF = 0.02
res@tmYLLabelFontHeightF = 0.02
res@gsnStringFontHeightF = 0.025

res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -2
res@cnMaxLevelValF = 2
res@cnLevelSpacingF = 0.2

res@gsnRightString = "K"

do i = 0, 11
    res@gsnLeftString =  Num(i)+" "+Month(i)
    plot(i) = gsn_csm_contour_map(wks, bias_clm(i,:,:), res)
end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
pres@pmLabelBarOrthogonalPosF = -0.015
pres@lbLabelFontHeightF = 0.012
gsn_panel(wks, plot, (/4,3/), pres)


end