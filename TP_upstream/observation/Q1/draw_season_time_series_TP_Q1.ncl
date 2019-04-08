
begin

year = ispan(1979, 2013, 1)

diri = "/home/yangsong3/wuyt/data/"
fili = "q1.mon.mean.1979.2013.nc" ; NCEP_DOE 

f = addfile(diri+fili, "r")
q1 = f->q1(:,:,::-1,:)
printVarSummary(q1)

Cp = 1004.07 ; J/(Kg*K)
q1 = q1 * Cp ; units W/Kg

;----------------integration-----------------------------
plev = q1&lev*100.
psfc = max(plev)
ptop = min(plev)
dP = dpres_plevel_Wrap(plev, psfc, ptop, 0)
printVarSummary(dP)

dp = conform_dims(dimsizes(q1), dP, 1)
q1dp = q1 * dp
copy_VarCoords(q1, q1dp)
printVarSummary(q1dp)
Q1 = dim_sum_n_Wrap(q1dp(:,{1000:1},:,:), 1)
copy_VarCoords(q1(:,0,:,:), Q1)
printVarSummary(Q1)
Q1 = Q1/9.8  ; units W/m2
;---------------------------------------------------------
Season = (/"MAM","JJA","SON","DJF"/)

q1Sea = month_to_seasonN(Q1, Season)
printVarSummary(q1Sea)
q1Sea(2,2,:,:) = q1@_FillValue
q1Sea(3,3,:,:) = q1@_FillValue
q1Sea(3,22,:,:) = q1@_FillValue

;------------------------------------
 f2 = "~/wuyt/data/TOPO_1500.nc"
 topo_in = addfile(f2,"r")
 topo = topo_in->TP

 xlat = topo_in->LATIXY(:,0)
 xlon = topo_in->LONGXY(0,:)
 topo!0 = "lsmlat"
 topo!1 = "lsmlon"
 topo&lsmlat = xlat
 topo&lsmlon = xlon
;------------------------------------

TP = linint2_Wrap(xlon, xlat, topo, True, Q1&lon, Q1&lat, 0)
tp = conform_dims(dimsizes(q1Sea), TP, (/2,3/))

q1Sea = mask(q1Sea, tp.ge.1500, True)
print(q1Sea(0,0,{35},{90}))

q1sea = wgt_areaave_Wrap(q1Sea(:,:,{20:45},{80:110}), 1, 1, 0)
printVarSummary(q1sea)

q1season = dtrend_msg(year, q1sea, False, False)
;write_table("Q1_TP80E_dtrend.txt", "w", [/"MAM","JJA","SON","DJF"/], "%10s %10s %10s %10s")
;write_table("Q1_TP80E_dtrend.txt", "a", [/q1sea(0,:),q1sea(1,:),q1sea(2,:),q1sea(3,:)/], "%10.2f %10.2f %10.2f %10.2f")

q1season = dim_standardize_n_Wrap(q1season, 1, 1)
printVarSummary(q1season)

write_table("Q1_TP80E_dtrend_standardize.txt", "w", [/"MAM","JJA","SON","DJF"/], "%10s %10s %10s %10s")
write_table("Q1_TP80E_dtrend_standardize.txt", "a", [/q1season(0,:),q1season(1,:),q1season(2,:),q1season(3,:)/], "%10.2f %10.2f %10.2f %10.2f")

;========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_upstream/climatology/"
figname = "draw_season_time_series_TP_Q1"

wks = gsn_open_wks("pdf", figpath+figname)
nplot = 4
plot = new(nplot, graphic)

res = True
res@gsnFrame = False
res@gsnDraw  = False

res@vpWidthF = 0.7
res@vpHeightF = 0.45

res@xyLineThicknesses = 3
res@xyLineColors = "black"

res@trXMinF = min(year)
res@trXMaxF = max(year)

res@gsnRightString = "W/m~S~2~N~"

res@gsnYRefLine = (/-0.8,0.8/)
res@gsnYRefLineColor = "green4"
res@gsnYRefLineThicknesses = 3

do i = 0, dimsizes(Season)-1
    res@gsnLeftString = Season(i)
    plot(i) = gsn_csm_xy(wks, year, q1season(i,:), res)
    ;gsn_csm_y(wks, q1sea(i,:), res)
end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelMainString = "Q1 time evolution"


gsn_panel(wks, plot, (/2,2/), pres)
end