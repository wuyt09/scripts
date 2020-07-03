
begin

yrStrt = 1979
yrLast = 2015
year = ispan(yrStrt, yrLast, 1)
season = (/"MAM","SON"/)

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
;===========integration========================================
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
TTT = dim_sum_n_Wrap(Tdp(:,{500:200},:,:), 1)/dim_sum_n_Wrap(dP(:,{500:200},:,:),1)
copy_VarCoords(T(:,0,:,:), TTT)
printVarSummary(TTT)
printMinMax(TTT, True)

;===============================================================

Tclm = clmMonTLL(TTT)
Tano = calcMonAnomTLL(TTT, Tclm)
;print(Tano(20,{25:40},{75:105}))
printVarSummary(Tano)

t = month_to_seasonN(Tano, season)
printVarSummary(t)
t = mask(t, topo.ge.1500, True)

ttt = wgt_areaave_Wrap(t(:,:,{25:40},{75:105}), 1, 1, 0)
printVarSummary(ttt)

t_dtr = dtrend_msg_n(year, ttt, False, False, 1)
t_std = dim_standardize_n_Wrap(t_dtr, 1, 1)

printVarSummary(t_std)
write_table("TTT_TP_dtrend_standardize_MAM_SON_ERAI.txt", "w", [/"MAM","SON"/], "%10s %10s")
write_table("TTT_TP_dtrend_standardize_MAM_SON_ERAI.txt", "a", [/t_std(0,:),t_std(1,:)/], "%10.2f %10.2f")
;=================================================================
pp = new((/2,dimsizes(t_std(0,:))/), float)
pp(0,:) = t_std(0,:)
pp(1,:) = t_std(1,:)
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_Spr_Aut_heating/obs/seasonal/"
figname = "draw_time_series_TTT_MAM_SON"

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
res@xyLineColors = (/"green4","goldenrod2"/)
res@xyDashPattern = 0

res@trXMinF = min(year)
res@trXMaxF = max(year)

res@gsnRightString = "K"
res@gsnStringFontHeightF = 0.025

res@gsnYRefLine = (/0.8,-0.8/)
res@gsnYRefLineColor = "gray54"
res@gsnYRefLineThicknesses = 3
res@gsnYRefLineDashPattern = 2


plot(0) = gsn_csm_xy(wks, year, pp, res)
draw(plot(0))
frame(wks)
end