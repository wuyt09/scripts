
begin
lonstr = 75
lonend = 105
latstr = 25
latend = 40

;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)
printVarSummary(topo)

;===================================================================
; OBS grid data
;===================================================================
yrStrt = 1980
yrLast = 2015
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

; ts_obs = dtrend_msg_n(time, ts_obs, False, False, 0)
ts_obs_sea = clmMonTLL(ts_obs)
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

; diriobs = "/home/yangsong3/data-observation/CHCN_CAMS_Gridded_2m_Temperature_Land/"
; filiobs = "air.mon.mean.nc"
; fobs = addfile(diriobs+filiobs, "r")
; Time := fobs->time
; YYYY := cd_calendar(Time, -1)/100
; iYYYY := ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
; ts_obs = fobs->air(iYYYY,::-1,:)

; ts_obs_sea = clmMonTLL(ts_obs)
; topo_obs = linint2_Wrap(topo&lon, topo&lat, topo, True, ts_obs&lon, ts_obs&lat, 0)
; ts_obs_sea = mask(ts_obs_sea, topo_obs.ge.1500, True)
; printVarSummary(ts_obs_sea)
; TP_obs_sea = wgt_areaave(ts_obs_sea(:,{latstr:latend},{lonstr:lonend}), 1, 1, 0) 
; printVarSummary(TP_obs_sea)
;===================================================================
; Reanalysis ERAI data
;===================================================================
diriera = "/home/yangsong3/data-observation/wuyt/ERAI/monthly/1x1/"
filiera = "t2m.mon.mean.7901.1804.nc"
; filiera = "skt.mon.mean.7917.nc"
fera = addfile(diriera+filiera, "r")
timer = fera->time
YYYY := cd_calendar(timer, -1)/100
iYYYY := ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
TS_era = short2flt(fera->t2m(iYYYY,::-1,:))
printVarSummary(TS_era)
;print(TS_era&latitude)
ts_era = clmMonTLL(TS_era)
printVarSummary(ts_era)
ts_era = mask(ts_era, topo.ge.1500, True)
TP_era = wgt_areaave(ts_era(:,{latstr:latend},{lonstr:lonend}), 1, 1, 0) 
printVarSummary(TP_era)

;===================================================================
; Reanalysis NCEP_DOE data
;===================================================================
dirin2 = "/home/yangsong3/data-observation/NCEP_DOE/monthly/"
filin2 = "air.2m.mon.mean.nc"
; filin2 = "skt.sfc.mon.mean.nc"
fn2 = addfile(dirin2+filin2, "r")
timen2 = fn2->time
YYYY := cd_calendar(timen2, -1)/100
iYYYY := ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
;TS_era = short2flt(fr->ts(iYYYY,:,:))
TS_n2 = short2flt(fn2->air(iYYYY,0,::-1,:))
printVarSummary(TS_n2)
ts_n2 = clmMonTLL(TS_n2)

topo_n2 = linint2_Wrap(topo&lon, topo&lat, topo, True, ts_n2&lon, ts_n2&lat, 0)
ts_n2 = mask(ts_n2, topo_n2.ge.1500, True)
printVarSummary(ts_n2)
TP_n2 = wgt_areaave(ts_n2(:,{latstr:latend},{lonstr:lonend}), 1, 1, 0) 
printVarSummary(TP_n2)

;===================================================================
; Reanalysis NCEP_NCAR data
;===================================================================
dirin1 = "/home/yangsong3/data-observation/NCEP_NCAR/monthly/"
filin1 = "air.sig995.mon.mean.nc"
fn1 = addfile(dirin1+filin1, "r")
timen1 = fn1->time
YYYY := cd_calendar(timen1, -1)/100
iYYYY := ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

TS_n1 = short2flt(fn1->air(iYYYY,::-1,:))
TS_n1 = TS_n1 + 273.15
ts_n1 = clmMonTLL(TS_n1)
printVarSummary(ts_n1)
topo_n1 = linint2_Wrap(topo&lon, topo&lat, topo, True, ts_n1&lon, ts_n1&lat, 0)
ts_n1 = mask(ts_n1, topo_n1.ge.1500, True)
TP_n1 = wgt_areaave(ts_n1(:,{latstr:latend},{lonstr:lonend}), 1, 1, 0) 
printVarSummary(TP_n1)

;===================================================================
; Reanalysis MERRA2 data
;===================================================================
dirime = "/home/yangsong3/data-observation/Yangsx/merra2/monthly/"
filime = "2m_t.198001-201811.mon.nc4"
fme = addfile(dirime+filime, "r")
timeme = fme->time
YYYYme = cd_calendar(timeme, -1)/100
iYYYYme = ind(YYYYme.ge.yrStrt .and. YYYYme.le.yrLast)
;TS_era = short2flt(fr->ts(iYYYY,:,:))
TS_me = short2flt(fme->T2M(iYYYYme,:,:))
ts_me = clmMonTLL(TS_me)
printVarSummary(ts_me)
topo_me = linint2_Wrap(topo&lon, topo&lat, topo, True, ts_me&lon, ts_me&lat, 0)
ts_me = mask(ts_me, topo_me.ge.1500, True)
TP_me = wgt_areaave(ts_me(:,{latstr:latend},{lonstr:lonend}), 1, 1, 0)
printVarSummary(TP_me)

;===================================================================
; Reanalysis ERA5 data
;===================================================================
diriera5 = "/home/yangsong3/data-observation/ERA5-monthly/surface/"
filiera5 = "t2m_mon_1x1_1979_2019.nc"
fera5 = addfile(diriera5+filiera5, "r")

timer5 = fera5->time
YYYY := cd_calendar(timer5, -1)/100
iYYYY := ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

TS_era5 = short2flt(fera->t2m(iYYYY,::-1,:))
printVarSummary(TS_era5)
;print(TS_era&latitude)
ts_era5 = clmMonTLL(TS_era5)
printVarSummary(ts_era5)
ts_era5 = mask(ts_era5, topo.ge.1500, True)
TP_era5 = wgt_areaave(ts_era5(:,{latstr:latend},{lonstr:lonend}), 1, 1, 0) 
printVarSummary(TP_era5)


;===================================================================
; Plot
;===================================================================

Solid = new((/6,12/), float)
Solid(0,:) = TP_obs_sea ;- 273.15
Solid(1,:) = TP_era5 - 273.15
Solid(2,:) = TP_era - 273.15
Solid(3,:) = TP_n2 - 273.15
Solid(4,:) = TP_n1 - 273.15
Solid(5,:) = TP_me - 273.15
; Solid = Solid - 273.15

; print(TP_era5)
; print(TP_era)
Month = ispan(1, 12, 1)
filipath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/20201228/"
filiname = "draw_12month_OBS_REANA"
wks = gsn_open_wks("pdf", filipath+filiname)

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

res@trYMinF = -16
res@trYMaxF = 20
res@tmYLMode = "Explicit"
res@tmYLValues = ispan(-15, 20, 5)
res@tmYLLabels = ispan(-15, 20, 5)
res@tmYROn = False
res@tiYAxisString = "Temperature (~S~o~N~C)"

Colors = (/"black","red", "blue", "purple","green4","orange2"/);
Thick = (/20,15,10,5,5,5/)
res@xyLineThicknesses = Thick
res@xyLineColors = Colors
res@xyDashPattern = 0 ;(/0,0,2/)


plot(0) = gsn_csm_xy(wks, Month, Solid, res)

;==================================================================================
    res_lines                   = True                  ; polyline mods desired
    res_lines@gsLineThicknessF  = 5                     ; line thicker
    

    res_text                    = True                  ; text mods desired
    res_text@txFontHeightF      = 0.018                 ; change text size
    res_text@txJust             = "CenterLeft"          ; text justification

    Texts  = (/"OBS","ERA5","ERAI","NCEP_DOE","NCEP_NCAR","MERRA2"/)

;============================================================
    gx = (/1.5,2.2/)
    text_pos = 2.3
    gy1 = (/18,18/)
    gy2 = (/16,16/)
    gy3 = (/14,14/)
    gy4 = (/12,12/)
    gy5 = (/10,10/)
    gy6 = (/8,8/)

    res_lines@gsLineColor = Colors(0)
    res_lines@gsLineDashPattern = 0
    dum11 = gsn_add_polyline(wks, plot(0), gx, gy1, res_lines)
    dum12 = gsn_add_text(wks, plot(0), Texts(0), text_pos, gy1, res_text)
    res_lines@gsLineColor = Colors(1)
    ;res_lines@gsLineDashPattern = 2
    dum21 = gsn_add_polyline(wks, plot(0), gx, gy2, res_lines)
    dum22 = gsn_add_text(wks, plot(0), Texts(1), text_pos, gy2, res_text)
    res_lines@gsLineColor = Colors(2)
    ;res_lines@gsLineDashPattern = 2
    dum31 = gsn_add_polyline(wks, plot(0), gx, gy3, res_lines)
    dum32 = gsn_add_text(wks, plot(0), Texts(2), text_pos, gy3, res_text)
    res_lines@gsLineColor = Colors(3)
    ;res_lines@gsLineDashPattern = 2
    dum41 = gsn_add_polyline(wks, plot(0), gx, gy4, res_lines)
    dum42 = gsn_add_text(wks, plot(0), Texts(3), text_pos, gy4, res_text)
    res_lines@gsLineColor = Colors(4)
    ;res_lines@gsLineDashPattern = 2
    dum51 = gsn_add_polyline(wks, plot(0), gx, gy5, res_lines)
    dum52 = gsn_add_text(wks, plot(0), Texts(4), text_pos, gy5, res_text)
    res_lines@gsLineColor = Colors(5)
    ;res_lines@gsLineDashPattern = 2
    dum61 = gsn_add_polyline(wks, plot(0), gx, gy6, res_lines)
    dum62 = gsn_add_text(wks, plot(0), Texts(5), text_pos, gy6, res_text)


draw(plot(0))
;draw(plot(1))
frame(wks)





end