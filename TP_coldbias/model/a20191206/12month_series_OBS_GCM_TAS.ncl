
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
; GCM 24 data
;===================================================================
dirigcm = "/home/yangsong3/data-model/wuyt/CMIP5/aaa-post/tas/"

models = (/"ACCESS1-0","ACCESS1-3","bcc-csm1-1","CanESM2","CCSM4","CESM1-BGC","CESM1-CAM5",\
        "CMCC-CM","CMCC-CMS","CNRM-CM5","CSIRO-Mk3-6-0","GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M",\
        "GISS-E2-H","GISS-E2-R","HadGEM2-CC","HadGEM2-ES","inmcm4","IPSL-CM5A-LR","IPSL-CM5A-MR",\
        "IPSL-CM5B-LR","MIROC5","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-M","NorESM1-ME"/)

nm = dimsizes(models)
TS_GCM = new((/nm,12,181,360/), float)
TS_GCM!0 = "model"
TS_GCM&model = models
do i = 0, nm-1
    filigcm = systemfunc("ls "+dirigcm+"tas_mon_"+models(i)+"_historical_r1i1p1_195101-200012.nc")
    fgcm = addfile(filigcm, "r")
    TS_GCMs = fgcm->tas
    TS_GCMs = dtrend_msg_n(TS_GCMs&time, TS_GCMs, False, False, 0)
    TS_GCM(i,:,:,:) = clmMonTLL(TS_GCMs)
end do
printVarSummary(TS_GCM)

TS_GCM = TS_GCM - 273.15
TS_GCM = mask(TS_GCM, topo.ge.1500, True)
TP_GCM = wgt_areaave_Wrap(TS_GCM(:,:,{latstr:latend},{lonstr:lonend}), 1, 1, 0)
TP_GCMA = dim_avg_n_Wrap(TP_GCM, 0)
printVarSummary(TP_GCMA)

print(TP_GCMA)

;===================================================================
; Plot
;===================================================================

Solid = new((/2,12/), float)
Solid(0,:) = (/TP_obs/)
Solid(1,:) = (/TP_GCMA/)


Month = ispan(1, 12, 1)
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/a20191206/"
figname = "12month_series_OBS_GCM_TAS"
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
res@tmXBLabels = Month
res@tmXTOn = False
res@tiXAxisString = "Month"

res@trYMinF = -20
res@trYMaxF = 20
res@tmYLMode = "Explicit"
res@tmYLValues = ispan(-20, 20, 5)
res@tmYLLabels = ispan(-20, 20, 5)
res@tmYROn = False
res@tiYAxisString = "Temperature (~S~o~N~C)"
;res@gsnLeftString = "Tas: grid2station"

;Colors = (/"deeppink2", "blue", "purple","green4","orange2", "black" /);
Colors = (/"black","forestgreen"/)
Patterns = (/0,0,0/)
res@xyLineThicknessF = 8
res@xyLineColors = Colors
res@xyDashPatterns = Patterns

plot(0) = gsn_csm_xy(wks, Month, Solid, res)

delete(res@xyLineColors)
delete(res@xyDashPatterns)
res@xyLineThicknessF = 3

gsn_define_colormap(wks, "WhiteYellowOrangeRed")
res@xyLineColors = ispan(25, 25+28*8, 8)
res@xyDashPattern = 16
plot(1) = gsn_csm_xy(wks, Month, TP_GCM, res)
overlay(plot(1), plot(0))


;==================================================================================
    res_lines                   = True                  ; polyline mods desired
    res_lines@gsLineThicknessF  = 5                     ; line thicker
    

    res_text                    = True                  ; text mods desired
    res_text@txFontHeightF      = 0.018                 ; change text size
    res_text@txJust             = "CenterLeft"          ; text justification

    Texts  = (/"OBS","GCM-MME"/)

;============================================================
    gx = (/1.5,2.2/)
    text_pos = 2.3
    gy1 = (/18,18/)
    gy2 = (/15,15/)
    ;gy3 = (/12,12/)
    gy4 = (/12,12/)
    ; gy5 = (/10,10/)
    ; gy6 = (/8,8/)

    res_lines@gsLineColor = Colors(0)
    res_lines@gsLineDashPattern = Patterns(0)
    dum11 = gsn_add_polyline(wks, plot(1), gx, gy1, res_lines)
    dum12 = gsn_add_text(wks, plot(1), Texts(0), text_pos, gy1, res_text)
    res_lines@gsLineColor = Colors(1)
    res_lines@gsLineDashPattern = Patterns(1)
    dum21 = gsn_add_polyline(wks, plot(1), gx, gy2, res_lines)
    dum22 = gsn_add_text(wks, plot(1), Texts(1), text_pos, gy2, res_text)
    ; res_lines@gsLineColor = Colors(2)
    ; res_lines@gsLineDashPattern = Patterns(2)
    ; dum31 = gsn_add_polyline(wks, plot(0), gx, gy3, res_lines)
    ; dum32 = gsn_add_text(wks, plot(0), Texts(2), text_pos, gy3, res_text)

    res_lines@gsLineColor = "red3"
    res_lines@gsLineDashPattern = 16
    dum41 = gsn_add_polyline(wks, plot(1), gx, gy4, res_lines)
    dum42 = gsn_add_text(wks, plot(1), "GCMs", text_pos, gy4, res_text)


draw(plot(1))
;draw(plot(1))
frame(wks)

end