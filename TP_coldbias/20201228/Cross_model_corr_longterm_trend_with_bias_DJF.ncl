
begin

lonstr = 75 ; 72 ;75
lonend = 105 ;136 ;105
latstr = 25 ;18 ;25
latend = 40 ;54 ; 40

season = "DJF"
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
yrLast = 2005
year = ispan(yrStrt, yrLast, 1)
nyear = dimsizes(year)
month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
nm = dimsizes(month)

ts_obsx = new((/nyear*nm,72,128/), float,-99.)

diriobs = "/home/yangsong3/data-observation/stationChina2000/Grid/t2m/"

do i = 0, nyear-1
    yr = yrStrt+i
    do j = 0, 11
        mn = sprinti("%0.2i", j+1)
        filiobs = "SURF_CLI_CHN_TEM_MON_GRID_0.5-MEAN-"+yr+mn+".txt"
        x = asciiread(diriobs+filiobs, -1, "float")
        ts_obsx(i*12+j,:,:) = reshape(x(6:), (/72,128/))
    end do
end do
printVarSummary(ts_obsx)

ts_obs = ts_obsx(:,::-1,:)

time = ispan(1, nyear*12, 1)
lat = fspan(18.25, 53.75, 72)
lon = fspan(72.25, 135.75, 128)
;print(lat)
;print(lon)

ts_obs!0 = "time"
ts_obs!1 = "lat"
ts_obs!2 = "lon"
ts_obs&time = time
ts_obs&lat = lat
ts_obs&lon = lon
ts_obs&lat@units = "degrees_north"
ts_obs&lon@units = "degrees_east"
printVarSummary(ts_obs)

ts_obs_sea = month_to_season(ts_obs, season)
printVarSummary(ts_obs_sea)

topox1 = ftopo->topo({18:54},{72:136})
lonx = topox1&lon
latx = topox1&lat
; printVarSummary(topox1)
topox2 = linint2_Wrap(lonx, latx, topox1, False, lon, lat, 0)
; print(topox2({35},{90}))

ts_obs_sea = mask(ts_obs_sea, topox2.ge.1500, True)
TP_obs_sea = wgt_areaave_Wrap(ts_obs_sea(:,{latstr:latend},{lonstr:lonend}), 1, 1, 0)
printVarSummary(TP_obs_sea)

TP_obs_sea_avg = avg(TP_obs_sea)
print(TP_obs_sea_avg)

;===================================================================
; GCM 20 data
;===================================================================
dirigcm = "/home/yangsong3/data-model/wuyt/CMIP5/"

models = (/"ACCESS1-0","ACCESS1-3","bcc-csm1-1","CanESM2","CCSM4","CESM1-BGC","CESM1-CAM5",\
        "CMCC-CM","CMCC-CMS","CNRM-CM5","CSIRO-Mk3-6-0","GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M",\
        "GISS-E2-H","GISS-E2-R","HadGEM2-CC","HadGEM2-ES","inmcm4","IPSL-CM5A-LR","IPSL-CM5A-MR",\
        "IPSL-CM5B-LR","MIROC5","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-M","NorESM1-ME"/)

nm = dimsizes(models)
TP_GCM5 = new((/nm,nyear/), float)
TP_GCM5!0 = "model"
TP_GCM5&model = models
do i = 0, nm-1
    filigcm = systemfunc("ls "+dirigcm+models(i)+"/tas_Amon_"+models(i)+"_historical_r1i1p1_18*01-20*12.nc")
    print(filigcm)
    fgcm = addfile(filigcm, "r")
    TIME := fgcm->time
    YYYY := cd_calendar(TIME, -1)/100
    iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.2005)
    TS_GCMs := fgcm->tas(iYYYY,:,:)
    TS_GCMs!1 = "lat"
    TS_GCMs!2 = "lon"
    TS_GCMs = TS_GCMs - 273.15
    TS_GCM := month_to_season(TS_GCMs, season)
    TOPO := linint2_Wrap(topo&lon, topo&lat, topo, True, TS_GCM&lon, TS_GCM&lat, 0)
    TS_GCM = mask(TS_GCM, TOPO.ge.1500, True)
    TP_GCM5(i,:) = wgt_areaave_Wrap(TS_GCM(:,{latstr:latend},{lonstr:lonend}), 1, 1, 0)
    print("end of "+models(i))
end do
printVarSummary(TP_GCM5)

TP_GCM5_avg = dim_avg_n_Wrap(TP_GCM5, 1)

TP_GCM5_bias = TP_GCM5_avg
TP_GCM5_bias = TP_GCM5_avg - TP_obs_sea_avg
printVarSummary(TP_GCM5_bias)

TP_GCM5_trend = TP_GCM5(:,0)
do i = 0, nm-1
    TP_GCM5_trend(i) = regline(year, TP_GCM5(i,:))
end do
TP_GCM5_trend = TP_GCM5_trend*nyear
printVarSummary(TP_GCM5_trend)

print(TP_GCM5_trend)

corr = escorc(TP_GCM5_bias, TP_GCM5_trend)
print(corr)

lines = new((/2,nm/), float)
lines(0,:) = TP_GCM5_bias
lines(1,:) = TP_GCM5_trend
;=========================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/20201228/"
figname = "Corss_model_corr_longterm_trend_with bias_"+season
; pwks = "png"
; pwks@wkHeight              = 3048
; pwks@wkWidth               = 6048
wks = gsn_open_wks("pdf", figpath+figname)

plot = new(2, graphic)

res = True
res@gsnDraw  = False
res@gsnFrame = False
res@gsnMaximize = True

res@vpWidthF = 0.7
res@vpHeightF = 0.4

; res@trXMinF = 1
; res@trXMaxF = 12
res@tmXBMode = "Explicit"
res@tmXBValues = ispan(1, nm, 1)
res@tmXBLabels = models
res@tmXTOn = False
res@tmXBLabelAngleF = 45
res@tmXBLabelFontHeightF = 0.007
; res@tiXAxisString = "Month"

res@trYMinF = -8
res@trYMaxF = 4
; res@tmYLMode = "Explicit"
; res@tmYLValues = ispan(-5, 4, 1)
res@tmYLPrecision = 2
; res@tmYLLabels = ispan(-20, 20, 5)
; res@tmYROn = False
res@tiYAxisString = ""
res@gsnLeftString = season
res@gsnRightString = "Corr: "+sprintf("%4.2f", corr)
res@gsnStringFontHeightF = 0.02

res@gsnYRefLine = 0
res@gsnYRefLineDashPattern = 2
res@gsnYRefLineColor = "grey40"
res@gsnYRefLineThicknessF = 8

;Colors = (/"deeppink2", "blue", "purple","green4","orange2", "black" /);
Colors = (/"black","red","forestgreen"/)
Patterns = (/0,0,0/)
res@xyLineThicknessF = 15
res@xyLineColors = Colors
res@xyDashPatterns = Patterns

plot(0) = gsn_csm_xy(wks, ispan(1, nm, 1), lines, res)

;==================================================================================
    res_lines                   = True                  ; polyline mods desired
    res_lines@gsLineThicknessF  = 5                     ; line thicker
    

    res_text                    = True                  ; text mods desired
    res_text@txFontHeightF      = 0.018                 ; change text size
    res_text@txJust             = "CenterLeft"          ; text justification

    Texts  = (/"Tas bias","Tas trend"/)

;============================================================
    gx = (/1.5,2.5/)
    text_pos = 2.8
    gy1 = (/3,3/)
    gy2 = (/2,2/)


    res_lines@gsLineThicknessF = 10
    res_lines@gsLineColor = Colors(0)
    res_lines@gsLineDashPattern = Patterns(0)
    dum11 = gsn_add_polyline(wks, plot(0), gx, gy1, res_lines)
    dum12 = gsn_add_text(wks, plot(0), Texts(0), text_pos, gy1, res_text)
    res_lines@gsLineColor = Colors(1)
    res_lines@gsLineDashPattern = Patterns(1)
    dum21 = gsn_add_polyline(wks, plot(0), gx, gy2, res_lines)
    dum22 = gsn_add_text(wks, plot(0), Texts(1), text_pos, gy2, res_text)

draw(plot(0))
end