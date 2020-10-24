;===========================================
; draw ensmble sst bias in CMIP5 models
; 12 month maps
;===========================================

begin

yrStrt = 1950
yrLast = 1999
year = ispan(yrStrt, yrLast, 1)
nyear = dimsizes(year)

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
SST_GCM = new((/nm,600,nlat,nlon/), float)
SST_GCM!0 = "model"
SST_GCM&model = models
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
    SST_GCM(i,:,:,:) = mask(TS_GCM, ismissing(sstobs), False)
    ; bias_sst5(i,:,:,:) = TS_GCM - sstobs
    print("end of "+models(i))
end do
TS_GCMA = dim_avg_n_Wrap(SST_GCM, 0)
copy_VarCoords(sstobs, TS_GCMA)
printVarSummary(TS_GCMA)

sstobs_12 = reshape(sstobs, (/nyear,12,nlat,nlon/))
sstmod_12 = reshape(TS_GCMA, (/nyear,12,nlat,nlon/))
copy_VarCoords(sstobs(0,:,:), sstobs_12(0,0,:,:))
copy_VarCoords(sstobs(0,:,:), sstmod_12(0,0,:,:))
obs_ave = dim_avg_n_Wrap(sstobs_12, 0)
obs_var = dim_variance_n_Wrap(sstobs_12, 0)
mod_ave = dim_avg_n_Wrap(sstmod_12, 0)
mod_var = dim_variance_n_Wrap(sstmod_12, 0)

prob = obs_ave
prob = -999.
do i = 0, 11
    prob(i,:,:) = ttest(obs_ave(i,:,:), obs_var(i,:,:), nyear, mod_ave(i,:,:), mod_var(i,:,:), nyear, False, False)
end do
printVarSummary(prob)

; copy_VarCoords(sstobs, bias_sst5(0,:,:,:))
; bias_sst5A = dim_avg_n_Wrap(bias_sst5, 0)
; printVarSummary(bias_sst5A)

bias_clm = obs_ave
bias_clm = mod_ave - obs_ave
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
wks = gsn_open_wks("pdf", figpath+figname)
plot = new(12, graphic)
plot_prob = plot


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


do i = 0, 11
    res@gsnLeftString =  Num(i)+" "+Month(i)
    plot(i) = gsn_csm_contour_map(wks, bias_clm(i,:,:), res)
    plot_prob(i) = gsn_csm_contour(wks, prob(i,:,:), res3)
    plot_prob(i) = ShadeLtContour(plot_prob(i),(1-sig),6)
    overlay(plot(i), plot_prob(i))
end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
pres@pmLabelBarOrthogonalPosF = -0.015
pres@lbLabelFontHeightF = 0.012
gsn_panel(wks, plot, (/6,2/), pres)


end