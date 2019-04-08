
begin

yrStrt = 1979
yrLast = 2016
year = ispan(yrStrt, yrLast, 1)

diri = "/home/yangsong3/wangzq/DATA/jra/"
fili = "air.jra55.mon.mean.nc"

f = addfile(diri+fili, "r")
TIME = f->time
YYYY = cd_calendar(TIME, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

T = short2flt(f->air(iYYYY,:,:,:))
printVarSummary(T)

;===========integration========================================
lev = T&level
plev = lev*100.
psfc = max(plev)
ptop = min(plev)

dP = dpres_plevel_Wrap(plev, psfc, ptop, 0)
dP!0 = "lev"
dP&lev = lev
printVarSummary(dP)

dp = conform_dims(dimsizes(T), dP, 1)
Tdp = T*dp
copy_VarCoords(T, Tdp)
TTT = dim_sum_n_Wrap(Tdp(:,{400:200},:,:), 1)/sum(dP({400:200}))
copy_VarCoords(T(:,0,:,:), TTT)
printVarSummary(TTT)
;===============================================================

Tclm = clmMonTLL(TTT)
Tano = calcMonAnomTLL(TTT, Tclm)
;print(Tano(20,{25:40},{75:105}))

t = month_to_season(Tano, "MAM")
printVarSummary(t)
;print(t(2,{25:40},{75:105}))

ttt = wgt_areaave_Wrap(t(:,{25:40},{75:105}), 1, 1, 0)
printVarSummary(ttt)

t_dtr = dtrend_msg(year, ttt, False, False)
t_std = dim_standardize_n_Wrap(t_dtr, 1, 0)
printVarSummary(t_std)
write_table("TTT_TP_dtrend_standardize_JRA55.txt", "w", [/"MAM"/], "%10s")
write_table("TTT_TP_dtrend_standardize_JRA55.txt", "a", [/t_std/], "%10.2f")
;=================================================================

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_upstream/observation/TTT/"
figname = "draw_MAM_time_series_TTT_JRA55"

wks = gsn_open_wks("pdf", figpath+figname)
nplot = 1
plot = new(nplot, graphic)

res = True
res@gsnFrame = False
res@gsnDraw  = False
res@gsnMaximize = True

res@vpWidthF = 0.7
res@vpHeightF = 0.45

res@xyLineThicknesses = 3
res@xyLineColors = "black"

res@trXMinF = min(year)
res@trXMaxF = max(year)

res@gsnRightString = "K"
res@gsnStringFontHeightF = 0.025

res@gsnYRefLine = (/0.8,-0.8/)
res@gsnYRefLineColor = "green4"
res@gsnYRefLineThicknesses = 3

res@gsnLeftString = "JRA55"
plot(0) = gsn_csm_xy(wks, year, t_std, res)
draw(plot(0))

end