
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
fili = "cloudcover.mon.mean.7916.r1.nc"
fx = addfile(diri+fili, "r")

TIME = fx->time
YYYY = cd_calendar(TIME, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

CC = short2flt(fx->cc(iYYYY,:,::-1,:))
CC = dtrend_msg_n(CC&time, CC, False, False, 0)
printVarSummary(CC)

filisp = "surf_pressure.mon.mean.7916.r1.nc"
fsp = addfile(diri+filisp, "r")
SP = short2flt(fsp->sp(iYYYY,::-1,:))

;===========integration========================================
lev = CC&level
plev = lev*100.
psfc = SP
ptop = min(plev)

dP = dpres_plevel_Wrap(plev, psfc, ptop, 0)
dP!1 = "lev"
dP&lev = lev
printVarSummary(dP)

CCDP = CC*dP
copy_VarCoords(CC, CCDP)

levsfc = 400
levtop = 200

cc_lay = dim_sum_n_Wrap(CCDP, 1)/dim_sum_n_Wrap(dP, 1)
copy_VarCoords(CC(:,0,:,:), cc_lay)
printVarSummary(cc_lay)

cc = clmMonTLL(cc_lay)
cc = mask(cc, topo.ge.1500, True)
printVarSummary(cc)

CX = wgt_areaave_Wrap(cc(:,{latstr:latend},{lonstr:lonend}), 1, 1, 0)
printVarSummary(CX)

; pp = new((/4,12/), float)
; pp(0,:) = td
; pp(1,:) = sd
; pp(2,:) = su
; pp(3,:) = tu
; pp = pp - 273.15

;===============================================================
; plot
;===============================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_Spr_Aut_heating/obs/"
figname = "draw_monthly_cloudcover"

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
res@tiYAxisString = " "
res@gsnYRefLine = 0
res@gsnYRefLineColor = "red"
res@gsnYRefLineDashPattern = 2
res@gsnYRefLineThicknessF = 5

Colors = (/"black", "blue", "purple", "green4", "orange2", "red"/);
res@xyLineThicknessF = 5
res@xyLineColors = Colors
res@xyDashPattern = 0 ;(/0,0,2/)
res@gsnLeftString = "Cloud cover"
res@gsnStringFontHeightF = 0.02

plot(0) = gsn_csm_xy(wks, Month, CX, res)
;overlay(plot(1), plot(0))


; ;==================================================================================
;     res_lines                   = True                  ; polyline mods desired
;     res_lines@gsLineThicknessF  = 5                     ; line thicker
    

;     res_text                    = True                  ; text mods desired
;     res_text@txFontHeightF      = 0.018                 ; change text size
;     res_text@txJust             = "CenterLeft"          ; text justification

;     Texts  = (/"layer 1","layer 2","layer 3","layer 4"/)

; ;============================================================
;     gx = (/1.5,2.2/)
;     text_pos = 2.3
;     gy1 = (/15,15/)
;     gy2 = (/13,13/)
;     gy3 = (/11,11/)
;     gy4 = (/9,9/)
;     ; gy5 = (/300,300/)
;     ; gy6 = (/280,280/)

;     res_lines@gsLineColor = Colors(0)
;     res_lines@gsLineDashPattern = 0
;     dum11 = gsn_add_polyline(wks, plot(0), gx, gy1, res_lines)
;     dum12 = gsn_add_text(wks, plot(0), Texts(0), text_pos, gy1, res_text)
;     res_lines@gsLineColor = Colors(1)
;     ;res_lines@gsLineDashPattern = 2
;     dum21 = gsn_add_polyline(wks, plot(0), gx, gy2, res_lines)
;     dum22 = gsn_add_text(wks, plot(0), Texts(1), text_pos, gy2, res_text)
;     res_lines@gsLineColor = Colors(2)
;     ;res_lines@gsLineDashPattern = 2
;     dum31 = gsn_add_polyline(wks, plot(0), gx, gy3, res_lines)
;     dum32 = gsn_add_text(wks, plot(0), Texts(2), text_pos, gy3, res_text)
;     res_lines@gsLineColor = Colors(3)
;     ;res_lines@gsLineDashPattern = 2
;     dum41 = gsn_add_polyline(wks, plot(0), gx, gy4, res_lines)
;     dum42 = gsn_add_text(wks, plot(0), Texts(3), text_pos, gy4, res_text)
;     res_lines@gsLineColor = Colors(4)
;     ; ;res_lines@gsLineDashPattern = 2
;     ; dum51 = gsn_add_polyline(wks, plot(0), gx, gy5, res_lines)
;     ; dum52 = gsn_add_text(wks, plot(0), Texts(4), text_pos, gy5, res_text)
;     ; res_lines@gsLineColor = Colors(5)
;     ; ;res_lines@gsLineDashPattern = 2
;     ; dum61 = gsn_add_polyline(wks, plot(0), gx, gy6, res_lines)
;     ; dum62 = gsn_add_text(wks, plot(0), Texts(5), text_pos, gy6, res_text)


draw(plot(0))
;draw(plot(1))
frame(wks)

end