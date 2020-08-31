
begin

yrStrt = 1979
yrLast = 2015
year = ispan(yrStrt, yrLast, 1)
nyear = dimsizes(year)
;season = (/"MAM","SON"/)
; season = (/"MAM","JJA","SON"/)

;================================================================
filitp = "~/wuyt/data/topo.sfc.1x1.nc"
ftp = addfile(filitp, "r")
topo = ftp->topo(::-1,:)
printVarSummary(topo)
;================================================================

diri = "/home/yangsong3/data-observation/ERA-interim/monthly/"
fili = "air.mon.mean.nc"

f = addfile(diri+fili, "r")
TIME = f->time
YYYY = cd_calendar(TIME, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

T = short2flt(f->t(iYYYY,:,::-1,:))
printVarSummary(T)
printMinMax(T, True)

filip = "Surface_pressure.mon.mean.nc"
fp = addfile(diri+filip, "r")
PS = short2flt(fp->sp(iYYYY,::-1,:))

; print(T&latitude)
nlat = dimsizes(T&latitude)
nlon = dimsizes(T&longitude)
;===========integration========================================
levs = 500
levt = 200

lev = T&level
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


t = new((/3,nyear,nlat,nlon/), float)
t_Mar = TTT(2::12,:,:)
t_Apr = TTT(3::12,:,:)
; t_MarApr = (t_Mar+t_Apr)/2
t_May = TTT(4::12,:,:)
t(0,:,:,:) = t_Mar
t(1,:,:,:) = t_Apr
t(2,:,:,:) = t_May ;t_MarApr
copy_VarCoords(TTT(0,:,:), t(0,0,:,:))
; t = month_to_seasonN(TTT, season)
; printVarSummary(t)
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
; write_table("TTT_TP_dtrend_standardize_Mar_Apr_May_ERAI.txt", "w", [/"MAM","JJA","SON"/], "%10s %10s %10s")
; write_table("TTT_TP_dtrend_standardize_Mar_Apr_May_ERAI.txt", "a", [/t_std(0,:),t_std(1,:),t_std(2,:)/], "%10.2f %10.2f %10.2f")
corr_34 = escorc(t_std(0,:), t_std(1,:))
corr_35 = escorc(t_std(0,:), t_std(2,:))
corr_45 = escorc(t_std(1,:), t_std(2,:))
ave3 = dim_avg_n_Wrap(t_std(0,:), 0)
ave4 = dim_avg_n_Wrap(t_std(1,:), 0)
ave5 = dim_avg_n_Wrap(t_std(2,:), 0)
var3 = dim_variance_n_Wrap(t_std(0,:), 0)
var4 = dim_variance_n_Wrap(t_std(1,:), 0)
var5 = dim_variance_n_Wrap(t_std(2,:), 0)
; print(corr_35)
prob34 = ttest(ave3, var3, nyear, ave4, var4, nyear, False, False)
prob35 = ttest(ave3, var3, nyear, ave5, var5, nyear, False, False)
prob45 = ttest(ave5, var5, nyear, ave4, var4, nyear, False, False)
; print(prob35)
; exit
print(sprintf("%4.2f", prob34))
print(sprintf("%4.2f", prob35))
print(sprintf("%4.2f", prob45))
; exit
;=================================================================
; pp = new((/2,dimsizes(t_std(0,:))/), float)
; pp(0,:) = t_std(0,:)
; pp(1,:) = t_std(1,:)
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20200704/"
figname = "TTT_time_series_Mar_Apr_May" ;Mar_Apr
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
res@gsnStringFontHeightF = 0.02


res@gsnRightString = "Corr_34="+sprintf("%4.2f", corr_34)+",  Corr_35="+sprintf("%4.2f", corr_35)+",  Corr_45="+sprintf("%4.2f", corr_45)
res@gsnRightStringFontHeightF = 0.015

res@gsnYRefLine = (/-1,1/) ; (/0.8,-0.8/)
res@gsnYRefLineColor = "gray54"
res@gsnYRefLineThicknesses = 3
res@gsnYRefLineDashPattern = 2
res@gsnLeftString = "TP TT ("+levs+"-"+levt+" hPa)"


plot(0) = gsn_csm_xy(wks, year, t_std, res)

;==================================================================================
    res_lines                   = True                  ; polyline mods desired
    res_lines@gsLineThicknessF  = 5                     ; line thicker
    

    res_text                    = True                  ; text mods desired
    res_text@txFontHeightF      = 0.018                 ; change text size
    res_text@txJust             = "CenterLeft"          ; text justification

    Texts  = (/"March","April","May"/)

; ;============================================================
    interval = 0.3
    gx = (/1985,1988/)
    text_pos = 1989
    gy1 = (/-2.0,-2.0/)
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