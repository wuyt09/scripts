
begin

yrStrt = 1979
yrLast = 2015
year = ispan(yrStrt, yrLast, 1)
nyear = dimsizes(year)
season = (/"MAM","SON"/)
; season = (/"MAM","JJA","SON"/)

;================================================================
filitp = "~/wuyt/data/topo.sfc.1x1.nc"
ftp = addfile(filitp, "r")
TOPO = ftp->topo(::-1,:)
printVarSummary(TOPO)
;================================================================

diri = "/home/yangsong3/wuyt/data/ERA-Interim/monthly/"
fili = "q1.mon.mean.1979.2016.ERAI.nc"

f = addfile(diri+fili, "r")
TIME = f->time
YYYY = cd_calendar(TIME, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

T = short2flt(f->q1(iYYYY,:,::-1,:))
printVarSummary(T)
printMinMax(T, True)

diri2 = "/home/yangsong3/data-observation/ERA-interim/monthly/"
filip = "Surface_pressure.mon.mean.nc"
fp = addfile(diri2+filip, "r")
PS0 = short2flt(fp->sp(iYYYY,::-1,:))
PS = linint2_Wrap(PS0&longitude, PS0&latitude, PS0, True, T&lon, T&lat, 0)
print(T&lat)
nlat = dimsizes(T&lat)
nlon = dimsizes(T&lon)
;===========integration========================================
levs = 600
levt = 300

lev = T&lev
plev = lev*100.
psfc = PS
ptop = min(plev)

dP = dpres_plevel_Wrap(plev, psfc, ptop, 0)
dP!1 = "lev"
dP&lev = lev
printVarSummary(dP)

Tdp = T*dP
copy_VarCoords(T, Tdp)
TTT = dim_sum_n_Wrap(Tdp(:,{levs:levt},:,:), 1)/dim_sum_n_Wrap(dP(:,{levs:levt},:,:),1)
copy_VarCoords(T(:,0,:,:), TTT)
printVarSummary(TTT)
printMinMax(TTT, True)

;===============================================================

; Tclm = clmMonTLL(TTT)
; Tano = calcMonAnomTLL(TTT, Tclm)
; ;print(Tano(20,{25:40},{75:105}))
; printVarSummary(Tano)


; t = new((/3,nyear,nlat,nlon/), float)
; t_Mar = TTT(2::12,:,:)
; t_Apr = TTT(3::12,:,:)
; t_MarApr = (t_Mar+t_Apr)/2
; t(0,:,:,:) = t_Mar
; t(1,:,:,:) = t_Apr
; t(2,:,:,:) = t_MarApr
; copy_VarCoords(TTT(0,:,:), t(0,0,:,:))

topo = linint2_Wrap(TOPO&lon, TOPO&lat, TOPO, True, T&lon, T&lat, 0)

t = month_to_seasonN(TTT, season)
printVarSummary(t)
t = mask(t, topo.ge.1500, True)

ttt = wgt_areaave_Wrap(t(:,:,{25:40},{75:105}), 1, 1, 0)
printVarSummary(ttt)

t_dtr = dtrend_msg_n(year, ttt, False, False, 1)
t_std = dim_standardize_n_Wrap(t_dtr, 1, 1)

printVarSummary(t_std)
; write_table("TTT_TP_dtrend_standardize_Mar_Apr_ERAI.txt", "w", [/"Mar","Apr","MarApr"/], "%10s %10s %10s")
; write_table("TTT_TP_dtrend_standardize_Mar_Apr_ERAI.txt", "a", [/t_std(0,:),t_std(1,:),t_std(2,:)/], "%10.2f %10.2f %10.2f")
; write_table("TTT_TP_dtrend_standardize_MAM_SON_ERAI.txt", "w", [/"MAM","SON"/], "%10s %10s")
; write_table("TTT_TP_dtrend_standardize_MAM_SON_ERAI.txt", "a", [/t_std(0,:),t_std(1,:)/], "%10.2f %10.2f")
write_table("Q1_TP_dtrend_standardize_MAM_SON_ERAI.txt", "w", [/"MAM","SON"/], "%10s %10s")
write_table("Q1_TP_dtrend_standardize_MAM_SON_ERAI.txt", "a", [/t_std(0,:),t_std(1,:)/], "%10.2f %10.2f")

;=================================================================
; pp = new((/2,dimsizes(t_std(0,:))/), float)
; pp(0,:) = t_std(0,:)
; pp(1,:) = t_std(1,:)
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20200704/"
figname = "Q1_time_series_MAM_SON" ;Mar_Apr
; Colors = (/"green4","goldenrod2","black"/)
Colors = (/"red","blue"/)

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

; res@gsnYRefLine = (/0.8,-0.8/)
res@gsnYRefLine = (/1,-1/)
res@gsnYRefLineColor = "gray54"
res@gsnYRefLineThicknesses = 3
res@gsnYRefLineDashPattern = 2
res@gsnLeftString = "Q1 ("+levs+"-"+levt+" hPa)"


plot(0) = gsn_csm_xy(wks, year, t_std, res)

;==================================================================================
    res_lines                   = True                  ; polyline mods desired
    res_lines@gsLineThicknessF  = 5                     ; line thicker
    

    res_text                    = True                  ; text mods desired
    res_text@txFontHeightF      = 0.018                 ; change text size
    res_text@txJust             = "CenterLeft"          ; text justification

    Texts  = (/"MAM","SON"/)

; ;============================================================
    interval = 0.3
    gx = (/1985,1988/)
    text_pos = 1989
    gy1 = (/2.7,2.7/)
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
    ; res_lines@gsLineColor = Colors(2)
    ; ;res_lines@gsLineDashPattern = Patterns(2)
    ; dum31 = gsn_add_polyline(wks, plot(0), gx, gy3, res_lines)
    ; dum32 = gsn_add_text(wks, plot(0), Texts(2), text_pos, gy3, res_text)

draw(plot(0))
frame(wks)
end