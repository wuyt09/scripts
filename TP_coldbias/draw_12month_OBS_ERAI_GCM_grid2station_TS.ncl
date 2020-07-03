
begin
; lonstr = 75
; lonend = 105
; latstr = 25
; latend = 40

;===================================================================
; OBS station data
;===================================================================
diriobs = "/home/yangsong3/wuyt/data/Station/GST/scripts/"
filiobs = "4-remove-some-stations.nc"
fobs = addfile(diriobs+filiobs, "r")
TS_obs = fobs->TS
TS_obs = TS_obs + 273.15
printVarSummary(TS_obs)
ts_obs = dim_avg_n_Wrap(TS_obs, 1)
TP_obs = dim_avg_n_Wrap(ts_obs, 0)
printVarSummary(TP_obs)
;print(TP_obs(0))
;staID = fobs->Station
;print(staID+"   "+ts_obs(:,0))
LAT = fobs->LAT
LON = fobs->LON
printVarSummary(LAT)
nstations = dimsizes(TS_obs(:,0,0))

;===================================================================
; Reanalysis ERAI data
;===================================================================
diriera = "/home/yangsong3/data-observation/wuyt/ERAI/monthly/1x1/"
filiera = "skt.mon.mean.7917.nc"
fera = addfile(diriera+filiera, "r")
TS_era = short2flt(fera->skt(:,::-1,:))
printVarSummary(TS_era)
;print(TS_era&latitude)
ts_era = clmMonTLL(TS_era)
printVarSummary(ts_era)

latera = fera->latitude(::-1)
lonera = fera->longitude
printVarSummary(latera)
nlatera = dimsizes(latera)
nlonera = dimsizes(lonera)
tp_era = new((/nstations,12/), float)
do ns = 0, nstations-1
    a = ind(abs(latera(:) - LAT(ns)) .le. 2)
    b = ind(abs(lonera(:) - LON(ns)) .le. 2)
    tp_era(ns,:) = dim_avg_n_Wrap(ts_era(:,a,b), (/1,2/))
    delete([/a,b/])
end do
printVarSummary(tp_era)
TP_era = dim_avg_n_Wrap(tp_era, 0)

;TP_era = wgt_areaave_Wrap(ts_era(:,{latstr:latend},{lonstr:lonend}), 1, 1, 0)
printVarSummary(TP_era)

;===================================================================
; GCM 24 data
;===================================================================
dirigcm = "/home/yangsong3/data-model/wuyt/CMIP5/aaa-post/ts/"

models = (/"ACCESS1-0","ACCESS1-3","bcc-csm1-1","CanESM2","CCSM4","CESM1-BGC","CESM1-CAM5",\
        "CMCC-CM","CMCC-CMS","CNRM-CM5","CSIRO-Mk3-6-0","GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M",\
        "GISS-E2-H","GISS-E2-R","HadGEM2-CC","HadGEM2-ES","inmcm4","IPSL-CM5A-LR","IPSL-CM5A-MR",\
        "IPSL-CM5B-LR","MIROC5","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-M","NorESM1-ME"/)

nm = dimsizes(models)
TS_GCM = new((/nm,12,181,360/), float)
TS_GCM!0 = "model"
TS_GCM&model = models
do i = 0, nm-1
    filigcm = systemfunc("ls "+dirigcm+"ts_mon_"+models(i)+"_historical_r1i1p1_195101-200012.nc")
    fgcm = addfile(filigcm, "r")
    TS_GCMs = fgcm->ts
    TS_GCMs = dtrend_msg_n(TS_GCMs&time, TS_GCMs, False, False, 0)
    TS_GCM(i,:,:,:) = clmMonTLL(TS_GCMs)
end do
printVarSummary(TS_GCM)

latgcm = fgcm->lat
longcm = fgcm->lon
printVarSummary(latgcm)
nlatgcm = dimsizes(latgcm)
nlongcm = dimsizes(longcm)

tp_gcm = new((/nm,nstations,12/), float)
do ns = 0, nstations-1
    a = ind(abs(latgcm(:) - LAT(ns)) .le. 2)
    b = ind(abs(longcm(:) - LON(ns)) .le. 2)
    tp_gcm(:,ns,:) = dim_avg_n_Wrap(TS_GCM(:,:,a,b), (/2,3/))
    delete([/a,b/])
end do
printVarSummary(tp_gcm)
TP_GCM = dim_avg_n_Wrap(tp_gcm, 1)
printVarSummary(TP_GCM)
TP_GCMA = dim_avg_n_Wrap(TP_GCM, 0)
printVarSummary(TP_GCMA)

;===================================================================
; Plot
;===================================================================

Solid = new((/3,12/), float)
Solid(0,:) = (/TP_obs/)
Solid(1,:) = (/TP_era/)
Solid(2,:) = (/TP_GCMA/)
Solid = Solid - 273.15 ; degC
TP_GCM = TP_GCM - 273.15


Month = ispan(1, 12, 1)
filipath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/"
filiname = "draw_12month_OBS_ERAI_GCM_grid2station_TS"
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
res@gsnLeftString = "Ts: grid2station"

;Colors = (/"deeppink2", "blue", "purple","green4","orange2", "black" /);
Colors = (/"black","blue","forestgreen"/)
Patterns = (/0,0,0/)
res@xyLineThicknessF = 5
res@xyLineColors = Colors
res@xyDashPatterns = Patterns


plot(0) = gsn_csm_xy(wks, Month, Solid, res)

delete(res@xyLineColors)
delete(res@xyDashPatterns)
res@xyLineThicknessF = 2
res@xyLineColor = "green2"
res@xyDashPattern = 2
plot(1) = gsn_csm_xy(wks, Month, TP_GCM, res)
overlay(plot(1), plot(0))


;==================================================================================
    res_lines                   = True                  ; polyline mods desired
    res_lines@gsLineThicknessF  = 5                     ; line thicker
    

    res_text                    = True                  ; text mods desired
    res_text@txFontHeightF      = 0.018                 ; change text size
    res_text@txJust             = "CenterLeft"          ; text justification

    Texts  = (/"OBS","ERAI","GCM-MME"/)

;============================================================
    gx = (/1.5,2.2/)
    text_pos = 2.3
    gy1 = (/18,18/)
    gy2 = (/16,16/)
    gy3 = (/14,14/)
    gy4 = (/12,12/)
    ; gy5 = (/10,10/)
    ; gy6 = (/8,8/)

    res_lines@gsLineColor = Colors(0)
    res_lines@gsLineDashPattern = Patterns(0)
    dum11 = gsn_add_polyline(wks, plot(0), gx, gy1, res_lines)
    dum12 = gsn_add_text(wks, plot(0), Texts(0), text_pos, gy1, res_text)
    res_lines@gsLineColor = Colors(1)
    res_lines@gsLineDashPattern = Patterns(1)
    dum21 = gsn_add_polyline(wks, plot(0), gx, gy2, res_lines)
    dum22 = gsn_add_text(wks, plot(0), Texts(1), text_pos, gy2, res_text)
    res_lines@gsLineColor = Colors(2)
    res_lines@gsLineDashPattern = Patterns(2)
    dum31 = gsn_add_polyline(wks, plot(0), gx, gy3, res_lines)
    dum32 = gsn_add_text(wks, plot(0), Texts(2), text_pos, gy3, res_text)

    res_lines@gsLineColor = "green2"
    res_lines@gsLineDashPattern = 2
    dum41 = gsn_add_polyline(wks, plot(0), gx, gy4, res_lines)
    dum42 = gsn_add_text(wks, plot(0), "GCMs", text_pos, gy4, res_text)


;draw(plot(0))
draw(plot(1))
frame(wks)





end