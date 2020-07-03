
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
    filisr = "F2000_"+casename(i)+".FSNS.mon.0601-3012.nc"
    fililr = "F2000_"+casename(i)+".FLNS.mon.0601-3012.nc"
    filild = "F2000_"+casename(i)+".FLDS.mon.0601-3012.nc"
    filisd = "F2000_"+casename(i)+".FSDS.mon.0601-3012.nc"

    fsh = addfile(diri+filild, "r")
    flh = addfile(diri+filisd, "r")
    fsr = addfile(diri+filisr, "r")
    flr = addfile(diri+fililr, "r")

    STRD = short2flt(fsh->FLDS)
    SSRD = short2flt(flh->FSDS)
    SSR  = short2flt(fsr->FSNS)
    STR  = short2flt(flr->FLNS)
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

    pp(i,0,:) = sd
    pp(i,1,:) = su
    pp(i,2,:) = tu
    pp(i,3,:) = td
end do
printVarSummary(pp)

ppp = new((/8,12/), float)
ppp(0:3,:) = pp(0,:,:)
ppp(4:7,:) = pp(1,:,:)

;===============================================================
; plot
;===============================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_Spr_Aut_heating/model/"
figname = "draw_monthly_rad_balance_CTRL_NOQ1"

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

res@trYMinF = 0
res@trYMaxF = 350
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
res@gsnLeftString = "Surface radiation balance"
res@gsnStringFontHeightF = 0.02

plot(0) = gsn_csm_xy(wks, Month, ppp, res)

;==================================================================================
    res_lines                   = True                  ; polyline mods desired
    res_lines@gsLineThicknessF  = 5                     ; line thicker
    

    res_text                    = True                  ; text mods desired
    res_text@txFontHeightF      = 0.018                 ; change text size
    res_text@txJust             = "CenterLeft"          ; text justification

    Texts  = (/"downward solar","upward solar","upward thermal","downward thermal","CTRL","NOQ1"/)

;============================================================
    gx = (/1.5,2.2/)
    text_pos = 2.3
    gy1 = (/340,340/)
    gy2 = (/320,320/)
    gy3 = (/300,300/)
    gy4 = (/280,280/)
    gy5 = (/260,260/)
    gy6 = (/240,240/)

    res_lines@gsLineColor = Colors(0)
    res_lines@gsLineDashPattern = 0
    dum11 = gsn_add_polyline(wks, plot(0), gx, gy1, res_lines)
    dum12 = gsn_add_text(wks, plot(0), Texts(0), text_pos, gy1, res_text)
    res_lines@gsLineColor = Colors(1)
    ;res_lines@gsLineDashPattern = 2
    dum21 = gsn_add_polyline(wks, plot(0), gx, gy2, res_lines)
    dum22 = gsn_add_text(wks, plot(0), Texts(1), text_pos, gy2, res_text)
    res_lines@gsLineColor = Colors(2)
    ;res_lines@gsLineDashPattern = 0
    dum31 = gsn_add_polyline(wks, plot(0), gx, gy3, res_lines)
    dum32 = gsn_add_text(wks, plot(0), Texts(2), text_pos, gy3, res_text)
    res_lines@gsLineColor = Colors(3)
    ;res_lines@gsLineDashPattern = 2
    dum41 = gsn_add_polyline(wks, plot(0), gx, gy4, res_lines)
    dum42 = gsn_add_text(wks, plot(0), Texts(3), text_pos, gy4, res_text)
    res_lines@gsLineColor = "black"
    ;res_lines@gsLineDashPattern = 2
    dum51 = gsn_add_polyline(wks, plot(0), gx, gy5, res_lines)
    dum52 = gsn_add_text(wks, plot(0), Texts(4), text_pos, gy5, res_text)
    res_lines@gsLineColor = "black"
    res_lines@gsLineDashPattern = 2
    dum61 = gsn_add_polyline(wks, plot(0), gx, gy6, res_lines)
    dum62 = gsn_add_text(wks, plot(0), Texts(5), text_pos, gy6, res_text)


; pres = True
; pres@gsnMaximize = True
; gsn_panel(wks, plot, (/2,1/), pres)
draw(plot)
frame(wks)

end