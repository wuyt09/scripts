
begin
lonstr = 75
lonend = 105
latstr = 25
latend = 40

;===================================================================
; OBS station data
;===================================================================



;===================================================================
; Reanalysis ERAI data
;===================================================================
diriera = "/home/yangsong3/data-observation/wuyt/ERAI/monthly/1x1/"
;filir = "t2m.mon.mean.7901.1804.nc"
filiera = "skt.mon.mean.7917.nc"
fera = addfile(diriera+filiera, "r")
; timer = fera->time
; TIME = cd_calendar(timer, -1)
; print(TIME)
;TS_era = short2flt(fr->ts(0:467,:,:))
TS_era = short2flt(fera->skt(:,::-1,:))
printVarSummary(TS_era)
;print(TS_era&latitude)
ts_era = clmMonTLL(TS_era)
printVarSummary(ts_era)
TP_era = wgt_areaave(ts_era(:,{latstr:latend},{lonstr:lonend}), 1, 1, 0) 
printVarSummary(TP_era)

;===================================================================
; Reanalysis NCEP_DOE data
;===================================================================
dirin2 = "/home/yangsong3/data-observation/NCEP_DOE/monthly/"
;filir = "t2m.mon.mean.7901.1804.nc"
filin2 = "skt.sfc.mon.mean.nc"
fn2 = addfile(dirin2+filin2, "r")
timen2 = fn2->time
YYYY = cd_calendar(timen2, -1)/100
iYYYY = ind(YYYY.ge.1979 .and. YYYY.le.2018)
;TS_era = short2flt(fr->ts(iYYYY,:,:))
TS_n2 = short2flt(fn2->skt(iYYYY,::-1,:))
printVarSummary(TS_n2)
ts_n2 = clmMonTLL(TS_n2)
printVarSummary(ts_n2)
TP_n2 = wgt_areaave(ts_n2(:,{latstr:latend},{lonstr:lonend}), 1, 1, 0) 
printVarSummary(TP_n2)

;===================================================================
; Reanalysis NCEP_NCAR data
;===================================================================
dirin1 = "/home/yangsong3/data-observation/NCEP_NCAR/ltmmean/"
filin1 = "skt.sfc.mon.ltm.nc"
fn1 = addfile(dirin1+filin1, "r")
timen1 = fn1->time
ts_n1 = short2flt(fn1->skt(:,::-1,:))
ts_n1 = ts_n1 + 273.15
printVarSummary(ts_n1)
TP_n1 = wgt_areaave(ts_n1(:,{latstr:latend},{lonstr:lonend}), 1, 1, 0) 
printVarSummary(TP_n1)

;===================================================================
; Reanalysis NCEP_NCAR data
;===================================================================
diriN1 = "/home/yangsong3/data-observation/NCEP_NCAR/monthly/"
filiN1 = "skt.sfc.mon.mean.nc"
fN1 = addfile(diriN1+filiN1, "r")
timeN1 = fN1->time
YYYYN1 = cd_calendar(timeN1, -1)/100
iYYYYN1 = ind(YYYYN1.ge.1948 .and. YYYYN1.le.2018)
;TS_era = short2flt(fr->ts(iYYYY,:,:))
TS_N1 = short2flt(fN1->skt(iYYYYN1,::-1,:))
ts_N1 = clmMonTLL(TS_N1)
ts_N1 = ts_N1 + 273.15
printVarSummary(ts_N1)
TP_N1 = wgt_areaave(ts_N1(:,{latstr:latend},{lonstr:lonend}), 1, 1, 0) 
printVarSummary(TP_N1)

;===================================================================
; Reanalysis MERRA2 data
;===================================================================
dirime = "/home/yangsong3/data-observation/Yangsx/merra2/monthly/"
filime = "skint.1980-201808.mon.nc4"
fme = addfile(dirime+filime, "r")
timeme = fme->time
YYYYme = cd_calendar(timeme, -1)/100
iYYYYme = ind(YYYYme.ge.1980 .and. YYYYme.le.2017)
;TS_era = short2flt(fr->ts(iYYYY,:,:))
TS_me = short2flt(fme->TS(iYYYYme,:,:))
ts_me = clmMonTLL(TS_me)
printVarSummary(ts_me)
TP_me = wgt_areaave(ts_me(:,{latstr:latend},{lonstr:lonend}), 1, 1, 0) 
printVarSummary(TP_me)

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
    filigcm = systemfunc("ls "+dirigcm+"ts_12mon_"+models(i)+"_historical_r1i1p1_195101-200012.nc")
    fgcm = addfile(filigcm, "r")
    TS_GCM(i,:,:,:) = fgcm->ts
end do
printVarSummary(TS_GCM)
TP_GCM = wgt_areaave(TS_GCM(:,:,{latstr:latend},{lonstr:lonend}), 1, 1, 0)
TP_GCMA = dim_avg_n_Wrap(TP_GCM, 0)
printVarSummary(TP_GCMA)


;===================================================================
; Plot
;===================================================================

Solid = new((/6,12/), float)
Solid(0,:) = TP_era
Solid(1,:) = TP_n2
Solid(2,:) = TP_n1
Solid(3,:) = TP_N1
Solid(4,:) = TP_me
Solid(5,:) = TP_GCMA
Solid = Solid - 273.15 ; degC
TP_GCM = TP_GCM - 273.15


Month = ispan(1, 12, 1)
filipath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/"
filiname = "draw_12month_OBS_REANA_GCM"
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

Colors = (/"deeppink2", "blue", "purple","green4","orange2", "black" /);
res@xyLineThicknessF = 3
res@xyLineColors = Colors
res@xyDashPattern = 0 ;(/0,0,2/)


plot(0) = gsn_csm_xy(wks, Month, Solid, res)

delete(res@xyLineColors)
delete(res@xyDashPatterns)
res@xyLineThicknessF = 2
res@xyLineColor = "green4"
res@xyDashPattern = 2
plot(1) = gsn_csm_xy(wks, Month, TP_GCM, res)
;overlay(plot(1), plot(0))


;==================================================================================
    res_lines                   = True                  ; polyline mods desired
    res_lines@gsLineThicknessF  = 5                     ; line thicker
    

    res_text                    = True                  ; text mods desired
    res_text@txFontHeightF      = 0.018                 ; change text size
    res_text@txJust             = "CenterLeft"          ; text justification

    Texts  = (/"ERAI","NCEP2","NCEP1-ltm","NCEP1","MERRA2","GCM-A"/)

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