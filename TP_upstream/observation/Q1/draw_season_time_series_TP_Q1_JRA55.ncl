
begin

yrStrt = 1958
yrLast = 2014
year = ispan(yrStrt, yrLast, 1)

diri = "/home/yangsong3/wangzq/DATA/jra/"
fili = "JRA55.mon.Q1.1958-2014.Yanai.nc" ; JRA55

f = addfile(diri+fili, "r")
TIME = f->time
YYYY = cd_calendar(TIME, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

Q1 = f->Q1(iYYYY,:,:)    ; units W/m^2

;Q1 = f->Q1I
;Q1&year = year
printVarSummary(Q1)
;---------------------------------------------------------
Season = (/"MAM","JJA","SON","DJF"/)

q1Sea = month_to_seasonN(Q1, Season)
printVarSummary(q1Sea)

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

q1sea = wgt_areaave_Wrap(q1Sea(:,:,{20:45},{60:110}), 1, 1, 0)
printVarSummary(q1sea)

q1season = dtrend_msg(year, q1sea, False, False)
;write_table("Q1_TP80E_dtrend.txt", "w", [/"MAM","JJA","SON","DJF"/], "%10s %10s %10s %10s")
;write_table("Q1_TP80E_dtrend.txt", "a", [/q1sea(0,:),q1sea(1,:),q1sea(2,:),q1sea(3,:)/], "%10.2f %10.2f %10.2f %10.2f")

q1season = dim_standardize_n_Wrap(q1season, 1, 1)
printVarSummary(q1season)

write_table("Q1_TP_dtrend_standardize_JRA55.txt", "w", [/"MAM","JJA","SON","DJF"/], "%10s %10s %10s %10s")
write_table("Q1_TP_dtrend_standardize_JRA55.txt", "a", [/q1season(0,:),q1season(1,:),q1season(2,:),q1season(3,:)/], "%10.2f %10.2f %10.2f %10.2f")

;========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_upstream/climatology/"
figname = "draw_season_time_series_TP_Q1_JRA55"

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

;res@gsnYRefLine = (/-0.8,0.8/)
res@gsnYRefLine = (/-1,1/)
res@gsnYRefLineColor = "green4"
res@gsnYRefLineThicknesses = 3

do i = 0, dimsizes(Season)-1
    res@gsnLeftString = Season(i)
    plot(i) = gsn_csm_xy(wks, year, q1season(i,:), res)
    ;gsn_csm_y(wks, q1sea(i,:), res)
end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelMainString = "Q1 time evolution JRA55"


gsn_panel(wks, plot, (/2,2/), pres)
end