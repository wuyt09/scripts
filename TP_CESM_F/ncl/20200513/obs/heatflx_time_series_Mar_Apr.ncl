
begin

yrStrt = 1979
yrLast = 2015
year = ispan(yrStrt, yrLast, 1)
nyear = dimsizes(year)

;================================================================
filitp = "~/wuyt/data/topo.sfc.1x1.nc"
ftp = addfile(filitp, "r")
topo = ftp->topo(::-1,:)
printVarSummary(topo)
;================================================================

diri = "/home/yangsong3/data-observation/wuyt/ERAI/monthly/1x1/"
filish = "sshf.mon.mean.7916.r1.nc"
fililh = "slhf.mon.mean.7916.r1.nc"

fsh = addfile(diri+filish, "r")
TIME = fsh->time
YYYY = cd_calendar(TIME, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

flh = addfile(diri+fililh, "r")

SH = short2flt(fsh->sshf(iYYYY,::-1,:))
LH = short2flt(flh->slhf(iYYYY,::-1,:))
SH = -1.*SH
LH = -1.*LH
heat = SH
heat = SH+LH
heat@long_name = "total heat flux"
printVarSummary(heat)

nlon = dimsizes(heat&longitude)
nlat = dimsizes(heat&latitude)
;===============================================================

; Tclm = clmMonTLL(TTT)
; Tano = calcMonAnomTLL(TTT, Tclm)
; ;print(Tano(20,{25:40},{75:105}))
; printVarSummary(Tano)

t = new((/3,nyear,nlat,nlon/), float)
t_Mar = heat(2::12,:,:)
t_Apr = heat(3::12,:,:)
t_MarApr = (t_Mar+t_Apr)/2
t(0,:,:,:) = t_Mar
t(1,:,:,:) = t_Apr
t(2,:,:,:) = t_MarApr
copy_VarCoords(heat(0,:,:), t(0,0,:,:))
printVarSummary(t)
t = mask(t, topo.ge.1500, True)

ttt = wgt_areaave_Wrap(t(:,:,{25:40},{75:105}), 1, 1, 0)
printVarSummary(ttt)

t_dtr = dtrend_msg_n(year, ttt, False, False, 1)
t_std = dim_standardize_n_Wrap(t_dtr, 1, 1)

printVarSummary(t_std)
write_table("heatflx_TP_dtrend_standardize_Mar_Apr_ERAI.txt", "w", [/"Mar","Apr","MarApr"/], "%10s %10s %10s")
write_table("heatflx_TP_dtrend_standardize_Mar_Apr_ERAI.txt", "a", [/t_std(0,:),t_std(1,:),t_std(2,:)/], "%10.2f %10.2f %10.2f")
;=================================================================
; pp = new((/2,dimsizes(t_std(0,:))/), float)
; pp(0,:) = t_std(0,:)
; pp(1,:) = t_std(1,:)
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20200513/obs/"
figname = "heatflx_time_series_Mar_Apr"
Colors = (/"green4","goldenrod2","black"/)

wks = gsn_open_wks("pdf", figpath+figname)
nplot = 1
plot = new(nplot, graphic)

res = True
res@gsnFrame = False
res@gsnDraw  = False
res@gsnMaximize = True

res@vpWidthF = 0.7
res@vpHeightF = 0.45

res@xyLineThicknessF = 5
res@xyLineColors = Colors
res@xyDashPattern = 0

res@trXMinF = min(year)
res@trXMaxF = max(year)

res@gsnRightString = ""
res@gsnStringFontHeightF = 0.025

res@gsnYRefLine = (/0.8,-0.8/)
res@gsnYRefLineColor = "gray54"
res@gsnYRefLineThicknesses = 3
res@gsnYRefLineDashPattern = 2
res@gsnLeftString = "HeatFlux"

plot(0) = gsn_csm_xy(wks, year, t_std, res)

;==================================================================================
    res_lines                   = True                  ; polyline mods desired
    res_lines@gsLineThicknessF  = 5                     ; line thicker
    

    res_text                    = True                  ; text mods desired
    res_text@txFontHeightF      = 0.018                 ; change text size
    res_text@txJust             = "CenterLeft"          ; text justification

    Texts  = (/"March","April","Mar+Apr"/)

; ;============================================================
    interval = 0.4
    gx = (/1980,1983/)
    text_pos = 1984
    gy1 = (/-2.4,-2.4/)
    gy2 = gy1-interval ;(/15,15/)
    gy3 = gy2-interval
    gy4 = gy3-interval
    print(gy1)

    res_lines@gsLineColor = Colors(0)
    ;res_lines@gsLineDashPattern = Patterns(0)
    dum11 = gsn_add_polyline(wks, plot(0), gx, gy1, res_lines)
    dum12 = gsn_add_text(wks, plot(0), Texts(0), text_pos, gy1, res_text)
    res_lines@gsLineColor = Colors(1)
    ;res_lines@gsLineDashPattern = Patterns(1)
    dum21 = gsn_add_polyline(wks, plot(0), gx, gy2, res_lines)
    dum22 = gsn_add_text(wks, plot(0), Texts(1), text_pos, gy2, res_text)
    res_lines@gsLineColor = Colors(2)
    ;res_lines@gsLineDashPattern = Patterns(2)
    dum31 = gsn_add_polyline(wks, plot(0), gx, gy3, res_lines)
    dum32 = gsn_add_text(wks, plot(0), Texts(2), text_pos, gy3, res_text)

draw(plot(0))
frame(wks)
end