
begin
lonstr = 75
lonend = 105
latstr = 25
latend = 40

;===================================================================
; read topo data
;===================================================================
filitp = "~/wuyt/data/topo.sfc.1x1.nc"
ftp = addfile(filitp, "r")
topo = ftp->topo(::-1,:)
printVarSummary(topo)

;===================================================================
; OBS station data
;===================================================================
diriobs = "/home/yangsong3/wuyt/data/Station/TEMP/scripts/"
filiobs = "4-remove-some-stations.nc"
fobs = addfile(diriobs+filiobs, "r")
TS_obsx = fobs->TAS
printVarSummary(TS_obsx)
LAT = fobs->LAT
LON = fobs->LON
HGT = fobs->HGT
hgt = dim_avg_n_Wrap(HGT, (/0,1/))
printVarSummary(LAT)

nstations = dimsizes(TS_obsx(:,0,0))
TS_obs = reshape(TS_obsx, (/nstations, 540/))
TS_obs!0 = "station"
TS_obs!1 = "time"
TS_obs&time = ispan(1, 540, 1)
TS_obs = dtrend_msg_n(TS_obs&time, TS_obs, False, False, 1)
TS_obsX = reshape(TS_obs, dimsizes(TS_obsx))

index = ind(hgt.ge.1500)
TS_OBS = TS_obsX(index,:,:)
ts_obs = dim_avg_n_Wrap(TS_OBS, (/0,1/))
ts_obs = ts_obs + 273.15
printVarSummary(ts_obs)

;===================================================================
; OBS station data -2
;===================================================================
 fx = "All_stations_1500m_12mon_clm.txt"
 ts_obs2 = asciiread(fx, 12, "float")
 print(ts_obs2)
 ts_obs2 = ts_obs2+273.15

;===================================================================
; Reanalysis ERAI data
;===================================================================
diriera = "/home/yangsong3/data-observation/wuyt/ERAI/monthly/1x1/"
filiera = "t2m.mon.mean.7901.1804.nc"
fera = addfile(diriera+filiera, "r")
TIMER = fera->time
YYYYera  = cd_calendar(TIMER, -1)/100
iYYYYera = ind(YYYYera.ge.1979 .and. YYYYera.le.2016)
TS_era = short2flt(fera->t2m(iYYYYera,::-1,:))
printVarSummary(TS_era)
TS_era = dtrend_msg_n(TS_era&time, TS_era, False, False, 0)
TS_era = mask(TS_era, topo.ge.1500, True)

ts_era = clmMonTLL(TS_era)
printVarSummary(ts_era)
TP_era = wgt_areaave(ts_era(:,{latstr:latend},{lonstr:lonend}), 1, 1, 0) 
printVarSummary(TP_era)

;===================================================================
; Reanalysis NCEP_DOE data
;===================================================================
dirin2 = "/home/yangsong3/data-observation/NCEP_DOE/monthly/"
filin2 = "air.2m.mon.mean.nc"
fn2 = addfile(dirin2+filin2, "r")
timen2 = fn2->time
YYYY = cd_calendar(timen2, -1)/100
iYYYY = ind(YYYY.ge.1979 .and. YYYY.le.2016)
;TS_era = short2flt(fr->ts(iYYYY,:,:))
TS_n2 = short2flt(fn2->air(iYYYY,0,::-1,:))
TS_n2 = dtrend_msg_n(TS_n2&time, TS_n2, False, False, 0)
topon2 = linint2_Wrap(topo&lon, topo&lat, topo, True, TS_n2&lon, TS_n2&lat, 0)
TS_n2 = mask(TS_n2, topon2.ge.1500, True)

ts_n2 = clmMonTLL(TS_n2)
printVarSummary(ts_n2)
TP_n2 = wgt_areaave(ts_n2(:,{latstr:latend},{lonstr:lonend}), 1, 1, 0) 
printVarSummary(TP_n2)

;===================================================================
; Reanalysis NCEP_NCAR data
;===================================================================
diriN1 = "/home/yangsong3/data-observation/NCEP_NCAR/monthly/"
filiN1 = "air.sig995.mon.mean.nc"
fN1 = addfile(diriN1+filiN1, "r")
timeN1 = fN1->time
YYYYN1 = cd_calendar(timeN1, -1)/100
iYYYYN1 = ind(YYYYN1.ge.1948 .and. YYYYN1.le.2016)
;TS_era = short2flt(fr->ts(iYYYY,:,:))
TS_N1 = short2flt(fN1->air(iYYYYN1,::-1,:))
TS_N1 = dtrend_msg_n(TS_N1&time, TS_N1, False, False, 0)
topoN1 = linint2_Wrap(topo&lon, topo&lat, topo, True, TS_N1&lon, TS_N1&lat, 0)
TS_N1 = mask(TS_N1, topoN1.ge.1500, True)

ts_N1 = clmMonTLL(TS_N1)
ts_N1 = ts_N1 + 273.15
printVarSummary(ts_N1)
TP_N1 = wgt_areaave(ts_N1(:,{latstr:latend},{lonstr:lonend}), 1, 1, 0) 
printVarSummary(TP_N1)

;===================================================================
; Reanalysis MERRA2 data
;===================================================================
dirime = "/home/yangsong3/data-observation/Yangsx/merra2/monthly/"
filime = "2m_t..198001-201811.mon.nc4"
fme = addfile(dirime+filime, "r")
timeme = fme->time
YYYYme = cd_calendar(timeme, -1)/100
iYYYYme = ind(YYYYme.ge.1980 .and. YYYYme.le.2017)
;TS_era = short2flt(fr->ts(iYYYY,:,:))
TS_me = short2flt(fme->T2M(iYYYYme,:,:))
TS_me = lonFlip(TS_me)
printVarSummary(TS_me)

TS_me = dtrend_msg_n(TS_me&time, TS_me, False, False, 0)
topome = linint2_Wrap(topo&lon, topo&lat, topo, True, TS_me&lon, TS_me&lat, 0)
TS_me = mask(TS_me, topome.ge.1500, True)

ts_me = clmMonTLL(TS_me)
printVarSummary(ts_me)
TP_me = wgt_areaave(ts_me(:,{latstr:latend},{lonstr:lonend}), 1, 1, 0) 
printVarSummary(TP_me)


;===================================================================
; Post process
;===================================================================
Month = ispan(1, 12, 1)
Solid = new((/6,12/), float)
Solid(0,:) = ts_obs
Solid(1,:) = TP_era
Solid(2,:) = TP_n2
Solid(3,:) = TP_N1
Solid(4,:) = TP_me
Solid(5,:) = ts_obs2
Solid = Solid - 273.15 ; degC

dt = Solid
do i = 1, 11
    dt(:,i) = Solid(:,i) - Solid(:,i-1)
end do
dt(:,0) = Solid(:,0) - Solid(:,11)

rc1 = Solid(:,0)
rc2 = rc1
do i = 0, dimsizes(rc1)-1
    rc1(i) = regline(Month(2:5), Solid(i,2:5))
    rc2(i) = regline(Month(7:10), Solid(i,7:10))
end do

filename = "./OBS_REANA_dt.txt"
system("rm -rf "+filename)
write_table(filename, "w", ([/"delta-temp", "OBS-remove", "ERAI", "NCEP2", "NCEP1", "MERRA2", "OBS-all"/]), "%10s %10s %10s %10s %10s %10s  %10s")
write_table(filename, "a", [/Month, dt(0,:), dt(1,:), dt(2,:), dt(3,:), dt(4,:), dt(5,:)/], "%10i %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f")
write_table(filename, "a", [/"rc_(3-6)", rc1(0), rc1(1), rc1(2), rc1(3), rc1(4), rc1(5)/], "%10s %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f")
write_table(filename, "a", [/"rc_(8-11)", rc2(0), rc2(1), rc2(2), rc2(3), rc2(4), rc2(5)/], "%10s %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f")



; filename = "./OBS_REANA_rc.txt"
; system("rm -rf "+filename)
; write_table(filename, "w", ([/"OBS-remove", "ERAI", "NCEP2", "NCEP1", "MERRA2", "OBS-all"/]), "%10s %10s %10s %10s %10s %10s")
; write_table(filename, "a", [/rc1(0), rc1(1), rc1(2), rc1(3), rc1(4), rc1(5)/], "%10.2f %10.2f %10.2f %10.2f %10.2f %10.2f")
; write_table(filename, "a", [/rc2(0), rc2(1), rc2(2), rc2(3), rc2(4), rc2(5)/], "%10.2f %10.2f %10.2f %10.2f %10.2f %10.2f")

;===================================================================
; Plot
;===================================================================

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_Spr_Aut_heating/obs/subseasonal/"
figname = "draw_monthly_surf_temp_OBS_REANA"

wks = gsn_open_wks("pdf", figpath+figname)

plot = new(1, graphic)

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

Colors = (/"black", "blue", "purple", "green4", "orange2", "red"/);
res@xyLineThicknessF = 5
res@xyLineColors = Colors
res@xyDashPattern = 0 ;(/0,0,2/)
res@gsnLeftString = "Surf air temp"
res@gsnStringFontHeightF = 0.03

plot(0) = gsn_csm_xy(wks, Month, Solid(0:4,:), res)
;overlay(plot(1), plot(0))


;==================================================================================
    res_lines                   = True                  ; polyline mods desired
    res_lines@gsLineThicknessF  = 5                     ; line thicker
    

    res_text                    = True                  ; text mods desired
    res_text@txFontHeightF      = 0.018                 ; change text size
    res_text@txJust             = "CenterLeft"          ; text justification

    Texts  = (/"OBS-station","ERAI","NCEP2","NCEP1","MERRA2","OBS-all stations"/)

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
    ; res_lines@gsLineColor = Colors(5)
    ; res_lines@gsLineDashPattern = 2
    ; dum61 = gsn_add_polyline(wks, plot(0), gx, gy6, res_lines)
    ; dum62 = gsn_add_text(wks, plot(0), Texts(5), text_pos, gy6, res_text)


draw(plot(0))
;draw(plot(1))
frame(wks)





end