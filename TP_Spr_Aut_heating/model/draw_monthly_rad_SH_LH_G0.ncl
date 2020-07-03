
begin

lonstr = 75
lonend = 105
latstr = 25
latend = 40
;================================================================
filitp = "~/wuyt/data/topo.sfc.1x1.nc"
ftp = addfile(filitp, "r")
TOPO = ftp->topo(::-1,:)
printVarSummary(TOPO)
;================================================================
Month = ispan(1, 12, 1)
casename = (/"CTRL","NOQ1"/)
pp = new((/dimsizes(casename),4,12/), float)

do i = 0, dimsizes(casename)-1

    diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(i)+"_F2000_F19G16_CAM4/pro/"
    filild = "F2000_"+casename(i)+".SHFLX.mon.0601-3012.nc"
    filisd = "F2000_"+casename(i)+".LHFLX.mon.0601-3012.nc"
    filisr = "F2000_"+casename(i)+".FSNS.mon.0601-3012.nc"
    fililr = "F2000_"+casename(i)+".FLNS.mon.0601-3012.nc"

    fsh = addfile(diri+filild, "r")
    flh = addfile(diri+filisd, "r")
    fsr = addfile(diri+filisr, "r")
    flr = addfile(diri+fililr, "r")

    STRD = short2flt(fsh->SHFLX)
    SSRD = short2flt(flh->LHFLX)
    SSRU = short2flt(fsr->FSNS)
    STRU = short2flt(flr->FLNS)

    STRD = dtrend_msg_n(STRD&time, STRD, False, False, 0)
    SSRD = dtrend_msg_n(SSRD&time, SSRD, False, False, 0)
    SSRU = dtrend_msg_n(SSRU&time, SSRU, False, False, 0)
    STRU = dtrend_msg_n(STRU&time, STRU, False, False, 0)
    printVarSummary(SSRU)

    strd = clmMonTLL(STRD)
    ssrd = clmMonTLL(SSRD)
    ssru = clmMonTLL(SSRU)
    stru = clmMonTLL(STRU)

    topo = linint2_Wrap(TOPO&lon, TOPO&lat, TOPO, True, strd&lon, strd&lat, 0)
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

    netR = su+tu
    G0 = netR - td - sd
    pp(i,0,:) = td ; SH
    pp(i,1,:) = sd ; LH
    pp(i,2,:) = (/netR/) ;
    pp(i,3,:) = (/G0/)
end do

ppp = new((/8,12/), float)
ppp(0:3,:) = pp(0,:,:)
ppp(4:7,:) = pp(1,:,:)

;===============================================================
; plot
;===============================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_Spr_Aut_heating/model/"
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

Colors = (/"black", "blue", "purple", "green4", "black", "blue", "purple", "green4"/);
res@xyLineThicknessF = 5
res@xyLineColors = Colors
res@xyDashPatterns = (/0,0,0,0,2,2,2,2/)
res@gsnLeftString = "Surface energy balance"
res@gsnStringFontHeightF = 0.02

plot(0) = gsn_csm_xy(wks, Month, ppp, res)
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
    gy1 = (/330,330/)
    gy2 = (/310,310/)
    gy3 = (/290,290/)
    gy4 = (/270,270/)
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