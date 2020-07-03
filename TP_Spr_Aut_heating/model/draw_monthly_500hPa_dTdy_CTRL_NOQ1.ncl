
begin
;=========================================================
filitp = "~/wuyt/data/topo.sfc.1x1.nc"
ftp = addfile(filitp, "r")
TOPO = ftp->topo(::-1,:)
printVarSummary(TOPO)
;=========================================================
Month = ispan(1, 12, 1)
casename = (/"CTRL","NOQ1"/)
dt = new((/dimsizes(casename),12/), float)

do j = 0, dimsizes(casename)-1
    diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(j)+"_F2000_F19G16_CAM4/pro/"
    fili = "F2000_"+casename(j)+".plevel.T.mon.0601-3012.nc"
    fx = addfile(diri+fili, "r")
    TS = fx->T(:,{500},:,:)
    TS = dtrend_msg_n(TS&time, TS, False, False, 0)
    t  = clmMonTLL(TS)
    printVarSummary(t)

    lon = tofloat(t&lon)
    lat = tofloat(t&lat)
    dlat = (lat(2)-lat(1))*0.0174533
    dy = 6378388.*dlat
    dtdx = t
    dtdx = center_finite_diff_n(t, dy, True, 0, 1)

    topo = linint2_Wrap(TOPO&lon, TOPO&lat, TOPO, True, dtdx&lon, dtdx&lat, 0)
    dtdx  = mask(dtdx, topo.ge.1500, True)
    dt(j,:) = wgt_areaave_Wrap(dtdx(:,{10:40},{70:110}), 1, 1, 0)
delete([/diri,fili,fx,TS,t,lon,lat,dtdx,dlat,dy/])
end do
printVarSummary(dt)
;============================================================

dk = dt
do i = 1, 11
    dk(:,i) = dt(:,i) - dt(:,i-1)
end do
dk(:,0) = dt(:,0) - dt(:,11)

rc1 = regline(Month(2:5), dt(0,2:5))
rc2 = regline(Month(7:10), dt(0,7:10))
print("CTRL: rc1 "+rc1+" rc2 "+rc2)

rc1 = regline(Month(2:5), dt(1,2:5))
rc2 = regline(Month(7:10), dt(1,7:10))
print("NOQ1: rc1 "+rc1+" rc2 "+rc2)

;===============================================================
; plot
;===============================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_Spr_Aut_heating/model/"
figname = "draw_monthly_500hPa_dTdy_CTRL_NOQ1"

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

; res@trYMinF = -16
; res@trYMaxF = 20
; res@tmYLMode = "Explicit"
; res@tmYLValues = ispan(-15, 20, 5)
; res@tmYLLabels = ispan(-15, 20, 5)
res@tmYROn = False
res@tiYAxisString = "Temperature (~S~o~N~C)"
res@gsnYRefLine = 0
res@gsnYRefLineColor = "red"
res@gsnYRefLineDashPattern = 2
res@gsnYRefLineThicknessF = 5

Colors = (/"black", "blue", "purple", "green4", "orange2", "red"/);
res@xyLineThicknessF = 5
res@xyLineColors = Colors
res@xyDashPattern = 0 ;(/0,0,2/)
res@gsnLeftString = "500hPa Temp"
res@gsnStringFontHeightF = 0.02

plot(0) = gsn_csm_xy(wks, Month, dt, res)
;overlay(plot(1), plot(0))


;==================================================================================
    res_lines                   = True                  ; polyline mods desired
    res_lines@gsLineThicknessF  = 5                     ; line thicker
    

    res_text                    = True                  ; text mods desired
    res_text@txFontHeightF      = 0.018                 ; change text size
    res_text@txJust             = "CenterLeft"          ; text justification

    Texts  = (/"CTRL","NOQ1"/)

;============================================================
    gx = (/1.5,2.2/)
    text_pos = 2.3
    gy1 = (/15,15/)
    gy2 = (/13,13/)
    gy3 = (/11,11/)
    gy4 = (/9,9/)
    ; gy5 = (/300,300/)
    ; gy6 = (/280,280/)

    res_lines@gsLineColor = Colors(0)
    res_lines@gsLineDashPattern = 0
    dum11 = gsn_add_polyline(wks, plot(0), gx, gy1, res_lines)
    dum12 = gsn_add_text(wks, plot(0), Texts(0), text_pos, gy1, res_text)
    res_lines@gsLineColor = Colors(1)
    ;res_lines@gsLineDashPattern = 2
    dum21 = gsn_add_polyline(wks, plot(0), gx, gy2, res_lines)
    dum22 = gsn_add_text(wks, plot(0), Texts(1), text_pos, gy2, res_text)
    ; res_lines@gsLineColor = Colors(2)
    ; ;res_lines@gsLineDashPattern = 2
    ; dum31 = gsn_add_polyline(wks, plot(0), gx, gy3, res_lines)
    ; dum32 = gsn_add_text(wks, plot(0), Texts(2), text_pos, gy3, res_text)
    ; res_lines@gsLineColor = Colors(3)
    ; ;res_lines@gsLineDashPattern = 2
    ; dum41 = gsn_add_polyline(wks, plot(0), gx, gy4, res_lines)
    ; dum42 = gsn_add_text(wks, plot(0), Texts(3), text_pos, gy4, res_text)
    ; res_lines@gsLineColor = Colors(4)
    ; ;res_lines@gsLineDashPattern = 2
    ; dum51 = gsn_add_polyline(wks, plot(0), gx, gy5, res_lines)
    ; dum52 = gsn_add_text(wks, plot(0), Texts(4), text_pos, gy5, res_text)
    ; res_lines@gsLineColor = Colors(5)
    ; ;res_lines@gsLineDashPattern = 2
    ; dum61 = gsn_add_polyline(wks, plot(0), gx, gy6, res_lines)
    ; dum62 = gsn_add_text(wks, plot(0), Texts(5), text_pos, gy6, res_text)


draw(plot(0))
;draw(plot(1))
frame(wks)



end