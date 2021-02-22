
begin
lonstr = 75 ; 72 ;75
lonend = 105 ;136 ;105
latstr = 25 ;18 ;25
latend = 40 ;54 ; 40

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
ts_obs_dtrend = dtrend_msg_n(time, ts_obsx, False, False, 0)
ts_obs = reshape(ts_obs_dtrend, (/nyr,nm,72,128/))

lat = fspan(18.25, 53.75, 72)
lon = fspan(72.25, 135.75, 128)
print(lat)
print(lon)

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

TS_obs = dim_avg_n_Wrap(ts_obs, 0)

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
; GCM 20 data
;===================================================================
dirigcm = "/home/yangsong3/data-model/wuyt/CMIP5/"

models = (/"ACCESS1-0","ACCESS1-3","bcc-csm1-1","CanESM2","CCSM4","CESM1-BGC","CESM1-CAM5",\
        "CMCC-CM","CMCC-CMS","CNRM-CM5","CSIRO-Mk3-6-0","GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M",\
        "GISS-E2-H","GISS-E2-R","HadGEM2-CC","HadGEM2-ES","inmcm4","IPSL-CM5A-LR","IPSL-CM5A-MR",\
        "IPSL-CM5B-LR","MIROC5","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-M","NorESM1-ME"/)

nm = dimsizes(models)
TP_GCM5 = new((/nm,12/), float)
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
    TS_GCMs = dtrend_msg_n(TS_GCMs&time, TS_GCMs, False, False, 0)
    TS_GCM := clmMonTLL(TS_GCMs)
    TOPO := linint2_Wrap(topo&lon, topo&lat, topo, True, TS_GCM&lon, TS_GCM&lat, 0)
    TS_GCM = mask(TS_GCM, TOPO.ge.1500, True)
    TP_GCM5(i,:) = wgt_areaave_Wrap(TS_GCM(:,{latstr:latend},{lonstr:lonend}), 1, 1, 0)
    print("end of "+models(i))
end do
printVarSummary(TP_GCM5)

; TS_GCM = TS_GCM - 273.15
; TS_GCM = mask(TS_GCM, topo.ge.1500, True)
; TP_GCM = wgt_areaave_Wrap(TS_GCM(:,:,{latstr:latend},{lonstr:lonend}), 1, 1, 0)
TP_GCMA5 = dim_avg_n_Wrap(TP_GCM5, 0)
printVarSummary(TP_GCMA5)

print(TP_GCMA5)


;===================================================================
; GCM 20 data
;===================================================================
;dirigcm = "/home/yangsong3/data-model/wuyt/CMIP6/"
dirigcm = "/home/yangsong3/CMIP6/"

models := (/"BCC-CSM2-MR","BCC-ESM1","CanESM5","CAS-ESM2-0","CESM2","CESM2-FV2","CESM2-WACCM",\
          "CESM2-WACCM-FV2","CIESM","E3SM-1-0","E3SM-1-1","E3SM-1-1-ECA","FGOALS-f3-L",\
          "FGOALS-g3","KIOST-ESM","MCM-UA-1-0","MIROC6","MRI-ESM2-0","NESM3","SAM0-UNICON"/)

nm := dimsizes(models)
TP_GCM6 = new((/nm,12/), float)
TP_GCM6!0 = "model"
TP_GCM6&model = models
do i = 0, nm-1
    filigcm = systemfunc("ls "+dirigcm+models(i)+"/historical/tas_Amon_"+models(i)+"_historical_r1i1p1f1_g*_185001-201*12.nc")
    fgcm = addfile(filigcm, "r")
    TIME := fgcm->time
    YYYY := cd_calendar(TIME, -1)/100
    iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.2005)
    TS_GCMs := fgcm->tas(iYYYY,:,:)
    TS_GCMs!1 = "lat"
    TS_GCMs!2 = "lon"
    TS_GCMs = TS_GCMs - 273.15
    TS_GCMs = dtrend_msg_n(TS_GCMs&time, TS_GCMs, False, False, 0)
    TS_GCM := clmMonTLL(TS_GCMs)
    TOPO := linint2_Wrap(topo&lon, topo&lat, topo, True, TS_GCM&lon, TS_GCM&lat, 0)
    TS_GCM = mask(TS_GCM, TOPO.ge.1500, True)
    TP_GCM6(i,:) = wgt_areaave_Wrap(TS_GCM(:,{latstr:latend},{lonstr:lonend}), 1, 1, 0)
    print("end of "+models(i))
end do
printVarSummary(TP_GCM6)

; TS_GCM = TS_GCM - 273.15
; TS_GCM = mask(TS_GCM, topo.ge.1500, True)
; TP_GCM = wgt_areaave_Wrap(TS_GCM(:,:,{latstr:latend},{lonstr:lonend}), 1, 1, 0)
TP_GCMA6 = dim_avg_n_Wrap(TP_GCM6, 0)
printVarSummary(TP_GCMA6)

; print(TP_GCMA6)
;===================================================================
diff5 = TP_GCMA5 - TP_obs
diff6 = TP_GCMA6 - TP_obs
ave5  = (diff5(0)+diff5(1)+diff5(11))/3
ave6  = (diff6(0)+diff6(1)+diff6(11))/3

ave5_s  = (diff5(5)+diff5(6)+diff5(7))/3
ave6_s  = (diff6(5)+diff6(6)+diff6(7))/3
print(ave5)
print(ave6)
print(ave5_s)
print(ave6_s)

;===================================================================
; Plot
;===================================================================

Solid = new((/3,12/), float)
Solid(0,:) = (/TP_obs/)
Solid(1,:) = (/TP_GCMA5/)
Solid(2,:) = (/TP_GCMA6/)


Month = ispan(1, 12, 1)
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/a20201012/"
figname = "Fig1"
pwks = "png"
pwks@wkHeight              = 3048
pwks@wkWidth               = 6048
wks = gsn_open_wks("pdf", figpath+figname)

plot = new(2, graphic)

res = True
res@gsnDraw  = False
res@gsnFrame = False
res@gsnMaximize = True

res@vpWidthF = 0.7
res@vpHeightF = 0.4

res@trXMinF = 1
res@trXMaxF = 12
res@tmXBMode = "Explicit"
res@tmXBValues = Month
res@tmXBLabels = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
res@tmXTOn = False
; res@tiXAxisString = "Month"

res@trYMinF = -20
res@trYMaxF = 20
res@tmYLMode = "Explicit"
res@tmYLValues = ispan(-20, 20, 5)
res@tmYLLabels = ispan(-20, 20, 5)
res@tmYROn = False
res@tiYAxisString = "Temperature (~S~o~N~C)"
;res@gsnLeftString = "Tas: grid2station"

;Colors = (/"deeppink2", "blue", "purple","green4","orange2", "black" /);
Colors = (/"black","purple","forestgreen"/)
Patterns = (/0,0,0/)
res@xyLineThicknessF = 15
res@xyLineColors = Colors
res@xyDashPatterns = Patterns

plot(0) = gsn_csm_xy(wks, Month, Solid, res)

delete(res@xyLineColors)
delete(res@xyDashPatterns)
res@xyLineThicknessF = 7

gsn_define_colormap(wks, "WhiteYellowOrangeRed")
res@xyLineColors = ispan(25, 25+28*8, 8)
res@xyDashPattern = 16
plot(1) = gsn_csm_xy(wks, Month, TP_GCM5, res)
overlay(plot(1), plot(0))


;==================================================================================
    res_lines                   = True                  ; polyline mods desired
    res_lines@gsLineThicknessF  = 5                     ; line thicker
    

    res_text                    = True                  ; text mods desired
    res_text@txFontHeightF      = 0.018                 ; change text size
    res_text@txJust             = "CenterLeft"          ; text justification

    Texts  = (/"OBS","CMIP5-MME","CMIP6-MME"/)

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
    dum11 = gsn_add_polyline(wks, plot(1), gx, gy1, res_lines)
    dum12 = gsn_add_text(wks, plot(1), Texts(0), text_pos, gy1, res_text)
    res_lines@gsLineColor = Colors(1)
    res_lines@gsLineDashPattern = Patterns(1)
    dum21 = gsn_add_polyline(wks, plot(1), gx, gy2, res_lines)
    dum22 = gsn_add_text(wks, plot(1), Texts(1), text_pos, gy2, res_text)
    res_lines@gsLineColor = Colors(2)
    res_lines@gsLineDashPattern = Patterns(2)
    dum31 = gsn_add_polyline(wks, plot(1), gx, gy3, res_lines)
    dum32 = gsn_add_text(wks, plot(1), Texts(2), text_pos, gy3, res_text)

    res_lines@gsLineColor = "red3"
    res_lines@gsLineThicknessF = 7
    res_lines@gsLineDashPattern = 16
    dum41 = gsn_add_polyline(wks, plot(1), gx, gy4, res_lines)
    dum42 = gsn_add_text(wks, plot(1), "CMIP5 GCMs", text_pos, gy4, res_text)


draw(plot(1))
;draw(plot(1))
frame(wks)

end