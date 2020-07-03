
begin
lonstr = 75
lonend = 105
latstr = 25
latend = 40

;===================================================================
; Reanalysis ERAI data
;===================================================================
diriera = "/home/yangsong3/data-observation/wuyt/ERAI/monthly/1x1/"
;filiera = "t2m.mon.mean.7901.1804.nc"
filiera = "skt.mon.mean.7917.nc"
fera = addfile(diriera+filiera, "r")
; timer = fera->time
; TIME = cd_calendar(timer, -1)
; print(TIME)
;TS_era = short2flt(fera->t2m(0:467,:,:))
TS_era = short2flt(fera->skt(:,::-1,:))
TS_era = TS_era - 273.15
printVarSummary(TS_era)
TS_era = dtrend_msg_n(TS_era&time, TS_era, False, False, 0)
ts_era = clmMonTLL(TS_era)
printVarSummary(ts_era)

latera = fera->latitude(::-1)
lonera = fera->longitude
printVarSummary(latera)
TP_era = wgt_areaave_Wrap(ts_era(:,{latstr:latend},{lonstr:lonend}), 1, 1, 1)
printVarSummary(TP_era)

;===================================================================
; OBS station data
;===================================================================
diriobs = "/home/yangsong3/wuyt/data/Station/GST/scripts/"
filiobs = "4-remove-some-stations.nc"
fobs = addfile(diriobs+filiobs, "r")
TS_obsx = fobs->TS
printVarSummary(TS_obsx)
LAT = fobs->LAT
LON = fobs->LON
printVarSummary(LAT)

nstations = dimsizes(TS_obsx(:,0,0))
TS_obs = reshape(TS_obsx, (/nstations, 540/))
TS_obs!0 = "station"
TS_obs!1 = "time"
TS_obs&time = ispan(1, 540, 1)
TS_obs = dtrend_msg_n(TS_obs&time, TS_obs, False, False, 1)
printVarSummary(TS_obs)

;===================================================================
; station to grid
;===================================================================
TS_OBS = TS_obs(time|:, station|:)
rscan = (/10,5,3/)
ts_obs = obj_anal_ic_Wrap(LON,LAT,TS_OBS,lonera,latera,rscan,False)
ts_obs := reshape(ts_obs, (/45,12,181,360/))
copy_VarCoords(TS_obsx(0,:,:), ts_obs(:,:,0,0))
copy_VarCoords(TS_era(0,:,:), ts_obs(0,0,:,:))
printVarSummary(ts_obs)

file_path = "./"
file_name = "station2grid_TP_TS.nc"
system("rm -rf "+file_path+file_name)
fout = addfile(file_path+file_name, "c")
fout->TS = ts_obs
; fout->LAT = lat
; fout->LON = lon

tp_obs = wgt_areaave_Wrap(ts_obs(:,:,{latstr:latend},{lonstr:lonend}), 1, 1, 1)
TP_obs = dim_avg_n_Wrap(tp_obs, 0) ; year
printVarSummary(TP_obs)
;print(TP_obs(0))
;staID = fobs->Station
;print(staID+"   "+ts_obs(:,0))
print(TP_obs)
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
TS_GCM = TS_GCM - 273.15
printVarSummary(TS_GCM)

latgcm = fgcm->lat
longcm = fgcm->lon
printVarSummary(latgcm)

TP_GCM = wgt_areaave_Wrap(TS_GCM(:,:,{latstr:latend},{lonstr:lonend}), 1, 1, 1)
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

Month = ispan(1, 12, 1)
filipath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/"
filiname = "draw_12month_OBS_ERAI_GCM_x"
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