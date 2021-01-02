
begin

lonstr = 60
lonend = 135
latstr = 5
latend = 50

latstr1 = 22
latend1 = 30
lonstr1 = 110
lonend1 = 120

;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)
printVarSummary(topo)

;===================================================================
; OBS t2m data
;===================================================================
yrStrt = 1979
yrLast = 2005
year = ispan(yrStrt, yrLast, 1)
nyear = dimsizes(year)
month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
nmonth = dimsizes(month)
season = (/"DJF","MAM"/)

ts_obsx = new((/nyear*nmonth,72,128/), float,-99.)

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

ts_obsx = ts_obsx(:,::-1,:)
time = ispan(1, nyear*12, 1)
ts_obs = dtrend_msg_n(time, ts_obsx, False, False, 0)
; ts_obs = reshape(ts_obs_dtrend, (/nyr,nm,72,128/))

lato = fspan(18.25, 53.75, 72)
lono = fspan(72.25, 135.75, 128)
print(lato)
print(lono)

ts_obs!0 = "time"
ts_obs!1 = "lat"
ts_obs!2 = "lon"
ts_obs&time = ispan(1, nyear*nmonth, 1)
ts_obs&lat = lato
ts_obs&lon = lono
ts_obs&lat@units = "degrees_north"
ts_obs&lon@units = "degrees_east"
printVarSummary(ts_obs)

TS_obsx = month_to_season(ts_obs, season(0))
TS_obs = dim_avg_n_Wrap(TS_obsx, 0)

topox1 = ftopo->topo({18:54},{72:136})
lonx = topox1&lon
latx = topox1&lat
printVarSummary(topox1)
topox2 = linint2_Wrap(lonx, latx, topox1, False, lono, lato, 0)
print(topox2({35},{90}))

TS_obs = mask(TS_obs, topox2.ge.1500, True)
TP_obs = wgt_areaave_Wrap(TS_obs({latstr:latend},{lonstr:lonend}), 1, 1, 0)
printVarSummary(TP_obs)



;========================================================
; obs prec
;========================================================
diriobs = "/home/yangsong3/data-observation/GPCP/"
filiobs = "precip.mon.mean.nc"
fobs = addfile(diriobs+filiobs, "r")
TIME = fobs->time
YYYYo = cd_calendar(TIME, -1)/100
iYYYYo = ind(YYYYo.ge.yrStrt .and. YYYYo.le.yrLast)
PR_OBS = fobs->precip(iYYYYo,:,:)
timeo = fobs->time(iYYYYo)
PR_OBS = dtrend_msg_n(timeo, PR_OBS, False, False, 0)
PR_OBSx = month_to_season(PR_OBS, season(1))
pr_obs = dim_avg_n_Wrap(PR_OBSx, 0)
pr_obs_ea = wgt_areaave_Wrap(pr_obs({latstr1:latend1},{lonstr1:lonend1}), 1., 1., 0)
printVarSummary(pr_obs_ea)

lat = pr_obs&lat
lon = pr_obs&lon
;===================================================================
; GCM 19 data
;===================================================================


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
; print(models)
; print(dimsizes(models))

dirigcm = "/home/yangsong3/wuyt/sysu/scripts/model_precip/CMIP5/"
nm = dimsizes(models)
diff = new((/nm/), float)
diff!0 = "model"
diff&model = models
PR_GCM6_ea = diff
do i = 0, nm-1
    filigcm = systemfunc("ls "+dirigcm+"pr_Amon_"+models(i)+"_historical_r1i1p1_*.nc")
    fgcm = addfile(filigcm, "r")
    TIME := fgcm->time
    YYYY := cd_calendar(TIME, -1)/100
    iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
    PR_GCMs = fgcm->pr(iYYYY,:,:)
    PR_GCMs = PR_GCMs * 86400
    times = fgcm->time(iYYYY)
    PR_GCMs = dtrend_msg_n(times, PR_GCMs, False, False, 0)
    PR_GCMs_jja := month_to_season(PR_GCMs, season(1))
    ; printVarSummary(PR_GCMs_jja)
    PR_GCM6 := dim_avg_n_Wrap(PR_GCMs_jja, 0)
    PR_GCM6_ea(i) = wgt_areaave_Wrap(PR_GCM6({latstr1:latend1},{lonstr1:lonend1}), 1., 1., 0)
    diff(i) = PR_GCM6_ea(i) - pr_obs_ea
    delete([/PR_GCMs/])
    print("end of prec "+models(i))
end do
printVarSummary(diff)


;===================================================================
; GCM 20 data
;===================================================================
dirigcm := "/home/yangsong3/data-model/wuyt/CMIP5/"

TP_ts_diff = new((/nm/), float)
TP_ts_diff!0 = "model"
TP_ts_diff&model = models
TP_GCM5 = TP_ts_diff
do i = 0, nm-1
    filigcm = systemfunc("ls "+dirigcm+models(i)+"/tas_Amon_"+models(i)+"_historical_r1i1p1_18*01-20*12.nc")
    ; print(filigcm)
    fgcm = addfile(filigcm, "r")
    TIME := fgcm->time
    YYYY := cd_calendar(TIME, -1)/100
    iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
    TS_GCMs := fgcm->tas(iYYYY,:,:)
    TS_GCMs!1 = "lat"
    TS_GCMs!2 = "lon"
    TS_GCMs = TS_GCMs - 273.15
    TS_GCMs = dtrend_msg_n(TS_GCMs&time, TS_GCMs, False, False, 0)
    TS_GCMx := month_to_season(TS_GCMs, season(0))
    TS_GCM := dim_avg_n_Wrap(TS_GCMx, 0)
    TOPO := linint2_Wrap(topo&lon, topo&lat, topo, True, TS_GCM&lon, TS_GCM&lat, 0)
    TS_GCM = mask(TS_GCM, TOPO.ge.1500, True)
    TP_GCM5(i) = wgt_areaave_Wrap(TS_GCM({latstr:latend},{lonstr:lonend}), 1, 1, 0)
    TP_ts_diff(i) = TP_GCM5(i) - TP_obs
    print("end of ts "+models(i))
end do
printVarSummary(TP_ts_diff)


;==================================================================
; sort by bias
;==================================================================

ip_tbias = dim_pqsort(TP_ts_diff,-1)
print(models(ip_tbias)+"  "+sprintf("%4.2f",TP_ts_diff(ip_tbias)))

lines = new((/4,nm/), float)
lines(0,:) = TP_GCM5
lines(1,:) = PR_GCM6_ea
lines(2,:) = TP_ts_diff(ip_tbias)
lines(3,:) = diff(ip_tbias)
corr = escorc(lines(0,:), lines(1,:))
print("Corr: "+sprintf("%4.2f", corr))

;===================================================================
; Plot
;===================================================================
filipath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20201216/"
filiname = "CMIP5_time_series_DJF_cold_bias_MAM_EA_prec"
wks = gsn_open_wks("pdf", filipath+filiname)
plot = new(1, graphic)
plot_topo = plot

res = True
res@gsnDraw  = False
res@gsnFrame = False
res@gsnMaximize = True

res@vpWidthF = 0.7
res@vpHeightF = 0.4

; res@trXMinF = 1
; res@trXMaxF = 12
; res@tmXBMode = "Explicit"
; res@tmXBValues = Month
; res@tmXBLabels = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
res@tmXTOn = False
; res@tiXAxisString = "Month"

; res@trYMinF = -20
; res@trYMaxF = 20
; res@tmYLMode = "Explicit"
; res@tmYLValues = ispan(-20, 20, 5)
; res@tmYLLabels = ispan(-20, 20, 5)
; res@tmYROn = False
; res@tiYAxisString = "Temperature (~S~o~N~C)"
;res@gsnLeftString = "Tas: grid2station"

;Colors = (/"deeppink2", "blue", "purple","green4","orange2", "black" /);
Colors = (/"black","red","black","red"/)
Patterns = (/0,0,16,16/)
res@xyLineThicknessF = 15
res@xyLineColors = Colors
res@xyDashPatterns = Patterns

plot(0) = gsn_csm_y(wks, lines, res)

;==================================================================================
    res_lines                   = True                  ; polyline mods desired
    res_lines@gsLineThicknessF  = 5                     ; line thicker
    

    res_text                    = True                  ; text mods desired
    res_text@txFontHeightF      = 0.018                 ; change text size
    res_text@txJust             = "CenterLeft"          ; text justification

    Texts  = (/"DJF TP Ts bias","MAM EA prec bias"/)

;============================================================
    gx = (/1.5,2.2/)
    text_pos = 2.3
    gy1 = (/18,18/)
    gy2 = (/15,15/)
    gy3 = (/12,12/)
    gy4 = (/9,9/)
    ; gy5 = (/10,10/)
    ; gy6 = (/8,8/)

    res_lines@gsLineThicknessF = 15
    res_lines@gsLineColor = Colors(0)
    res_lines@gsLineDashPattern = Patterns(0)
    dum11 = gsn_add_polyline(wks, plot(0), gx, gy1, res_lines)
    dum12 = gsn_add_text(wks, plot(0), Texts(0), text_pos, gy1, res_text)
    res_lines@gsLineColor = Colors(1)
    res_lines@gsLineDashPattern = Patterns(1)
    dum21 = gsn_add_polyline(wks, plot(0), gx, gy2, res_lines)
    dum22 = gsn_add_text(wks, plot(0), Texts(1), text_pos, gy2, res_text)


draw(plot(0))
;draw(plot(1))
frame(wks)
end