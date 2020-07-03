
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
filish = "sshf.mon.mean.7916.r1.nc"
fililh = "slhf.mon.mean.7916.r1.nc"
filisr = "net_solar_sur_mon_7916_1.nc"
fililr = "net_thermal_sur_mon_7916_1.nc"

fsh = addfile(diri+filish, "r")
flh = addfile(diri+fililh, "r")
fsr = addfile(diri+filisr, "r")
flr = addfile(diri+fililr, "r")

TIME = fsh->time
YYYY = cd_calendar(TIME, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

SSHF = short2flt(fsh->sshf(iYYYY,::-1,:))
SLHF = short2flt(flh->slhf(iYYYY,::-1,:))
SSR  = short2flt(fsr->ssr(iYYYY,::-1,:))
STR  = short2flt(flr->str(iYYYY,::-1,:))
SSHF = SSHF *(-1.)
SLHF = SLHF *(-1.)
SSHF = dtrend_msg_n(SSHF&time, SSHF, False, False, 0)
SLHF = dtrend_msg_n(SLHF&time, SLHF, False, False, 0)
SSR  = dtrend_msg_n(SSR&time, SSR, False, False, 0)
STR  = dtrend_msg_n(STR&time, STR, False, False, 0)
printVarSummary(SSR)

sshf = clmMonTLL(SSHF)
slhf = clmMonTLL(SLHF)
ssr  = clmMonTLL(SSR)
str  = clmMonTLL(STR)
sshf = mask(sshf, topo.ge.1500, True)
slhf = mask(slhf, topo.ge.1500, True)
ssr  = mask(ssr, topo.ge.1500, True)
str  = mask(str, topo.ge.1500, True)
printVarSummary(str)

sh = wgt_areaave_Wrap(sshf(:,{latstr:latend},{lonstr:lonend}), 1, 1, 0)
lh = wgt_areaave_Wrap(slhf(:,{latstr:latend},{lonstr:lonend}), 1, 1, 0)
sr = wgt_areaave_Wrap(ssr(:,{latstr:latend},{lonstr:lonend}), 1, 1, 0)
tr = wgt_areaave_Wrap(str(:,{latstr:latend},{lonstr:lonend}), 1, 1, 0)
printVarSummary(sr)

netR = sr
g0 = sr
netR = sr+tr
g0 = netR - sh - lh
printVarSummary(netR)

pp = new((/6,12/), float)
pp(0,:) = sh
pp(1,:) = lh
pp(2,:) = netR
pp(3,:) = g0
pp(4,:) = sr
pp(5,:) = tr

;===============================================================
; plot
;===============================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_Spr_Aut_heating/obs/"
figname = "draw_monthly_rad_SH_LH_G0"

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
res@gsnLeftString = "Surface energy balance"
res@gsnStringFontHeightF = 0.02

plot(0) = gsn_csm_xy(wks, Month, pp(0:3,:), res)
;overlay(plot(1), plot(0))


;==================================================================================
    res_lines                   = True                  ; polyline mods desired
    res_lines@gsLineThicknessF  = 5                     ; line thicker
    

    res_text                    = True                  ; text mods desired
    res_text@txFontHeightF      = 0.018                 ; change text size
    res_text@txJust             = "CenterLeft"          ; text justification

    Texts  = (/"SH","LH","netR","G (netR-SH-LH)","SR","LR"/)

;============================================================
    gx = (/1.5,2.2/)
    text_pos = 2.3
    ; gy1 = (/18,18/)
    ; gy2 = (/16,16/)
    ; gy3 = (/14,14/)
    ; gy4 = (/12,12/)
    ; gy5 = (/10,10/)
    ; gy6 = (/8,8/)
    gy1 = (/140,140/)
    gy2 = (/130,130/)
    gy3 = (/120,120/)
    gy4 = (/110,110/)
    gy5 = (/100,100/)
    gy6 = (/90,90/)

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