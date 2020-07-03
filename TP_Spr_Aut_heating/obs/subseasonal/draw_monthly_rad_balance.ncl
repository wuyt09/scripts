
begin

lonstr = 75
lonend = 105
latstr = 25
latend = 40

yrStrt = 1979
yrLast = 2016
year = ispan(yrStrt, yrLast, 1)
Month = ispan(1, 12, 1)

;================================================================
filitp = "~/wuyt/data/topo.sfc.1x1.nc"
ftp = addfile(filitp, "r")
topo = ftp->topo(::-1,:)
printVarSummary(topo)
;================================================================

diri = "/home/yangsong3/data-observation/wuyt/ERAI/monthly/1x1/"
filisr = "net_solar_sur_mon_7916_1.nc"
fililr = "net_thermal_sur_mon_7916_1.nc"
filild = "surface_thermal_downwards_monthly_1979_2018.nc"
filisd = "ssrd.mon.mean.7916.r1.nc"

fsh = addfile(diri+filild, "r")
flh = addfile(diri+filisd, "r")
fsr = addfile(diri+filisr, "r")
flr = addfile(diri+fililr, "r")

TIME = fsr->time
YYYY = cd_calendar(TIME, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

STRD = short2flt(fsh->strd(iYYYY,::-1,:))
SSRD = short2flt(flh->ssrd(iYYYY,::-1,:))
SSR  = short2flt(fsr->ssr(iYYYY,::-1,:))
STR  = short2flt(flr->str(iYYYY,::-1,:))
STRU = STRD - STR
SSRU = SSRD - SSR
copy_VarCoords(STRD, STRU)
copy_VarCoords(SSRD, SSRU)

STRD = dtrend_msg_n(STRD&time, STRD, False, False, 0)
SSRD = dtrend_msg_n(SSRD&time, SSRD, False, False, 0)
SSRU = dtrend_msg_n(SSRU&time, SSRU, False, False, 0)
STRU = dtrend_msg_n(STRU&time, STRU, False, False, 0)
printVarSummary(SSRU)

strd = clmMonTLL(STRD)
ssrd = clmMonTLL(SSRD)
ssru = clmMonTLL(SSRU)
stru = clmMonTLL(STRU)
strd = mask(strd, topo.ge.1500, True)
ssrd = mask(ssrd, topo.ge.1500, True)
ssru = mask(ssru, topo.ge.1500, True)
stru = mask(stru, topo.ge.1500, True)
printVarSummary(stru)

td = wgt_areaave_Wrap(strd(:,{latstr:latend},{lonstr:lonend}), 1, 1, 0)
sd = wgt_areaave_Wrap(ssrd(:,{latstr:latend},{lonstr:lonend}), 1, 1, 0)
su = wgt_areaave_Wrap(ssru(:,{latstr:latend},{lonstr:lonend}), 1, 1, 0)
tu = wgt_areaave_Wrap(stru(:,{latstr:latend},{lonstr:lonend}), 1, 1, 0)
printVarSummary(sd)

pp = new((/4,12/), float)
pp(0,:) = sd
pp(1,:) = su
pp(2,:) = tu
pp(3,:) = td


;===============================================================
; plot
;===============================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_Spr_Aut_heating/obs/"
figname = "draw_monthly_rad_balance"

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
res@tiYAxisString = "Flux (W/m~S~2~N~)"
res@gsnYRefLine = 0
res@gsnYRefLineColor = "red"
res@gsnYRefLineDashPattern = 2
res@gsnYRefLineThicknessF = 5

Colors = (/"black", "blue", "purple", "green4", "orange2", "red"/);
res@xyLineThicknessF = 5
res@xyLineColors = Colors
res@xyDashPattern = 0 ;(/0,0,2/)
res@gsnLeftString = "Surface radiation balance"
res@gsnStringFontHeightF = 0.02

plot(0) = gsn_csm_xy(wks, Month, pp, res)
;overlay(plot(1), plot(0))


;==================================================================================
    res_lines                   = True                  ; polyline mods desired
    res_lines@gsLineThicknessF  = 5                     ; line thicker
    

    res_text                    = True                  ; text mods desired
    res_text@txFontHeightF      = 0.018                 ; change text size
    res_text@txJust             = "CenterLeft"          ; text justification

    Texts  = (/"downward solar","upward solar","upward thermal","downward thermal"/)

;============================================================
    gx = (/1.5,2.2/)
    text_pos = 2.3
    gy1 = (/380,380/)
    gy2 = (/360,360/)
    gy3 = (/340,340/)
    gy4 = (/320,320/)
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
    res_lines@gsLineColor = Colors(2)
    ;res_lines@gsLineDashPattern = 2
    dum31 = gsn_add_polyline(wks, plot(0), gx, gy3, res_lines)
    dum32 = gsn_add_text(wks, plot(0), Texts(2), text_pos, gy3, res_text)
    res_lines@gsLineColor = Colors(3)
    ;res_lines@gsLineDashPattern = 2
    dum41 = gsn_add_polyline(wks, plot(0), gx, gy4, res_lines)
    dum42 = gsn_add_text(wks, plot(0), Texts(3), text_pos, gy4, res_text)
    res_lines@gsLineColor = Colors(4)
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