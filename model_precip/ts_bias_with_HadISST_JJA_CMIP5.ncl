
begin


yrStrt = 1950
yrLast = 1999
year = ispan(yrStrt, yrLast, 1)
nyear = dimsizes(year)
season = "JJA"

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

  obs_jjax = month_to_season(sstobs, season)
  obs_jja = dim_avg_n_Wrap(obs_jjax, 0)
  printVarSummary(obs_jja)

;===================================================================
; CMIP5 28 GCMs sst data
;===================================================================
dirigcm = "/home/yangsong3/data-model/wuyt/CMIP5/"

model1 = (/"ACCESS1-0","ACCESS1-3","bcc-csm1-1","CanESM2","CCSM4","CESM1-BGC","CESM1-CAM5",\
        "CMCC-CM","CMCC-CMS","CNRM-CM5","CSIRO-Mk3-6-0","GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M",\
        "GISS-E2-H","GISS-E2-R","HadGEM2-CC","HadGEM2-ES","inmcm4","IPSL-CM5A-LR","IPSL-CM5A-MR",\
        "IPSL-CM5B-LR","MIROC5","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-M","NorESM1-ME"/)

model2 = (/"ACCESS1-0","bcc-csm1-1","CanCM4","CanESM2","CCSM4","CESM1-BGC",\
          "CESM1-CAM5","CESM1-FASTCHEM","CESM1-WACCM","CNRM-CM5",\
          "FIO-ESM","GFDL-CM2p1","GFDL-ESM2M","GISS-E2-H","GISS-E2-R-CC",\
          "GISS-E2-R","HadCM3","inmcm4","IPSL-CM5A-LR","IPSL-CM5B-LR",\
          "MIROC4h","MIROC5","MIROC-ESM-CHEM","MIROC-ESM","MPI-ESM-LR","MPI-ESM-MR",\
          "MPI-ESM-P","MRI-CGCM3","NorESM1-M"/)

models = venn2_intersection(model1, model2)

nm = dimsizes(models)
diff_sst = new((/nm,nlat,nlon/), float)
diff_sst!0 = "model"
diff_sst&model = models
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
    SST_GCM = mask(TS_GCM, ismissing(sstobs), False)
    copy_VarCoords(TS_GCM, SST_GCM)
    SST_JJAx = month_to_season(SST_GCM, season)
    SST_JJA = dim_avg_n_Wrap(SST_JJAx, 0)
    diff_sst(i,:,:) = SST_JJA - obs_jja
    copy_VarCoords(obs_jja, diff_sst(0,:,:))
    print("end of "+models(i))
end do
printVarSummary(diff_sst)
diff_sstA = dim_avg_n_Wrap(diff_sst, 0)
;=========================================================
; plot 
;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/model_precip/"
figname = "ts_bias_with_ERSST_JJA_CMIP5"

pwks = "png"
pwks@wkHeight = 3048
pwks@wkWidth  = 6048
wks = gsn_open_wks("pdf", figpath+figname)
plot = new(1, graphic)
plot_prob = plot


res                       = True     ; plot mods desired
res@gsnDraw               = False
res@gsnFrame              = False

res@gsnMaximize           = True     ; Make plot fit the frame
res@cnFillOn              = True     ; turn on color fill
res@cnLinesOn             = False    ; turn of contour lines
;res@cnLevelSpacingF       = 0.5      ; contour spacing
res@cnFillPalette         = "hotcold_18lev"
; res@lbLabelBarOn          = False
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

res3@cnFillPattern = 17 ;(/17/)
res3@cnFillColor  =  "black" ;"green4"
res3@gsnLeftString = ""
res3@cnFillDotSizeF = 0.003
res3@cnFillScaleF  = 1 ;1.2

sig = 0.95
res@gsnRightString = "K" ;sig*100+"%"
res@gsnLeftString = "SST: CMIP5 MME - ERSST"


res@pmLabelBarOrthogonalPosF = 0.4

plot(0) = gsn_csm_contour_map(wks, diff_sstA, res)


gres = True
gres@gsLineColor = "red"
gres@gsLineThicknessF = 5.0
glat = (/-10,-10,0,0,-10/)
glon = (/90,110,110,90,90/)

dum = gsn_add_polyline(wks,plot(0),glon,glat,gres)

draw(plot(0))
end