
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

models1 = (/"ACCESS1-0","bcc-csm1-1","CanCM4","CanESM2","CCSM4","CESM1-BGC",\
          "CESM1-CAM5","CESM1-FASTCHEM","CESM1-WACCM","CMCC-CM","CNRM-CM5",\
          "EC-EARTH","FIO-ESM","GFDL-CM2p1","GFDL-ESM2M","GISS-E2-H","GISS-E2-R-CC",\
          "GISS-E2-R","HadCM3","inmcm4","IPSL-CM5A-LR","IPSL-CM5A-MR","IPSL-CM5B-LR",\
          "MIROC4h","MIROC5","MIROC-ESM-CHEM","MIROC-ESM","MPI-ESM-LR","MPI-ESM-MR",\
          "MRI-CGCM3","NorESM1-M"/)

models2 = (/"ACCESS1-0","ACCESS1-3","bcc-csm1-1","CanESM2","CCSM4","CESM1-BGC","CESM1-CAM5",\
        "CMCC-CM","CMCC-CMS","CNRM-CM5","CSIRO-Mk3-6-0","GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M",\
        "GISS-E2-H","GISS-E2-R","HadGEM2-CC","HadGEM2-ES","inmcm4","IPSL-CM5A-LR","IPSL-CM5A-MR",\
        "IPSL-CM5B-LR","MIROC5","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-M","NorESM1-ME"/)

models = venn2_intersection(models1,models2)

nm = dimsizes(models)
SST_GCM = new((/nm,600,nlat,nlon/), float)
SST_GCM!0 = "model"
SST_GCM&model = models
copy_VarCoords(sstobs, SST_GCM(0,:,:,:))
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


;====================================================
; Pattern correlation coeiffcient
;====================================================
; annual mean, tropical ocean

; Pattern corr for each model

SST_MME_bias = TS_GCMA - sstobs
copy_VarCoords(sstobs, SST_MME_bias)

SST_MME_bias_med = month_to_annual(SST_MME_bias, 1)
SST_MME_bias_new = dim_avg_n_Wrap(SST_MME_bias_med, 0)
printVarSummary(SST_MME_bias_new)

o_1d = ndtooned(SST_MME_bias_new({-30:30},:))
o_1d_bar = avg(o_1d)
printVarSummary(o_1d_bar)

corr = new(nm, float)
do i = 0, nm-1
    SST_GCM_bias = SST_GCM(i,:,:,:) - sstobs
    copy_VarCoords(sstobs, SST_GCM_bias)

    SST_GCM_med = month_to_annual(SST_GCM_bias, 1)
    SST_GCM_new = dim_avg_n_Wrap(SST_GCM_med, 0)

    y_1d = ndtooned(SST_GCM_new({-30:30},:))
    y_1d_bar = avg(y_1d)

    upper_term = (y_1d - y_1d_bar)*(o_1d - o_1d_bar)
    lower_term_1 = (y_1d - y_1d_bar)^2
    lower_term_2 = (o_1d - o_1d_bar)^2
    corr(i) = sum(upper_term)/sqrt(sum(lower_term_1)*sum(lower_term_2))
end do
print(models+"  "+sprintf("%4.2f", corr))


filipath = "./"
filiname = "calc_pattern_corr_SST_bias_GCM_with_MME_CMIP5.txt"

; write_table(filipath+filiname, "w", [/"models","corr_with_MME"/], "%10s %10s")
; write_table(filipath+filiname, "a", [/models,sprintf("%4.2f", corr)/], "%10.2f %10.2f")

write_table(filipath+filiname, "w", [/"corr_with_MME"/], "%10s")
write_table(filipath+filiname, "a", [/corr/], "%10.2f")

end
