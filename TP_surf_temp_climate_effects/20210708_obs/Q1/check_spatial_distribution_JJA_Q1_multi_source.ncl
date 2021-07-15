begin

yrStrt = 1979
yrLast = 2013
year = ispan(yrStrt, yrLast, 1)
nyear = dimsizes(year)
season = "JJA"

LEVs = 600
LEVe = 300

latstr = 25
latend = 40
lonstr = 75
lonend = 105

;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)
printVarSummary(topo)


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
lev1 = T&lev
lev1!0 = "level"
lev1&level = lev1

t_med = month_to_season(T, season)
t = dim_avg_n_Wrap(t_med, 0)

ERAI_profile = wgt_areaave_Wrap(t(:,{latstr:latend},{lonstr:lonend}), 1., 1., 0)
printVarSummary(ERAI_profile)


fili4 = "/home/yangsong3/data-observation/ERA5-monthly/pressure/thermal_equation/q1.monmean.1979.2019.nc"
f4 = addfile(fili4, "r")
TIME4 = f4->time
YYYY4 = cd_calendar(TIME4, -1)/100
iYYYY4 = ind(YYYY4.ge.yrStrt .and. YYYY4.le.yrLast)

X4 = short2flt(f4->q1(iYYYY4,:,:,:))
printVarSummary(X4)
lev4 = X4&lev
lev4!0 = "level"
lev4&level = lev4

x4_med = month_to_season(X4, season)
x4 = dim_avg_n_Wrap(x4_med, 0)

ERA5_profile = wgt_areaave_Wrap(x4(:,{latstr:latend},{lonstr:lonend}), 1., 1., 0)
printVarSummary(ERA5_profile)



fili2 = "/home/yangsong3/wuyt/data/q1.mon.mean.1979.2013.nc"
f2 = addfile(fili2, "r")
TIME2 = f2->time
YYYY2 = cd_calendar(TIME2, -1)/100
iYYYY2 = ind(YYYY2.ge.yrStrt .and. YYYY2.le.yrLast)

X2 = short2flt(f2->q1(iYYYY2,:,:,:))
printVarSummary(X2)
lev2 = X2&lev
lev2!0 = "level"
lev2&level = lev2

x2_med = month_to_season(X2, season)
x2 = dim_avg_n_Wrap(x2_med, 0)

NCEP_profile = wgt_areaave_Wrap(x2(:,{latstr:latend},{lonstr:lonend}), 1., 1., 0)
printVarSummary(NCEP_profile)

fili3 = "/home/yangsong3/wuyt/data/Q1.mon.mean.1979.2013.nc"
f3 = addfile(fili3, "r")
TIME3 = f3->time
YYYY3 = cd_calendar(TIME3, -1)/100
iYYYY3 = ind(YYYY3.ge.yrStrt .and. YYYY3.le.yrLast)

X3 = short2flt(f3->q1(iYYYY3,:,:,:))
printVarSummary(X3)
lev3 = X3&lev
lev3!0 = "level"
lev3&level = lev3

x3_med = month_to_season(X3, season)
x3 = dim_avg_n_Wrap(x3_med, 0)

NCEP2_profile = wgt_areaave_Wrap(x3(:,{latstr:latend},{lonstr:lonend}), 1., 1., 0)
printVarSummary(NCEP2_profile)


;====================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_surf_temp_climate_effects/20210708/Q1/"
figname = "check_Q1_multi_source"

wks = gsn_open_wks("pdf", figpath+figname)

plot = new(4, graphic)


 res = True
 res@gsnDraw  = False
 res@gsnFrame = False

 res@gsnXRefLine = 0
 res@gsnXRefLineDashPattern = 2
 res@gsnXRefLineColor = "grey"
 res@gsnXRefLineThicknessF = 3
 res@trXMinF = -6.
 res@trXMaxF = 6
 res@trYMinF   = 100
 res@trYMaxF   = 1000
 res@trYLog = True
 
 res@tmYLMode      = "Explicit"
 res@tmYLValues    = (/"1000","850","750","500","300","100"/)
 res@tmYLLabels    = (/"1000","850","750","500","300","100"/)
 res@xyLineColors = "black"
 res@xyLineThicknessF = 3
 res@xyDashPattern = 0

 res@trYReverse = True
 
 res@gsnRightString = "K"
 res@gsnStringFontHeightF = 0.025


 res@gsnLeftString = "ERAI"
 plot(0) = gsn_csm_xy(wks, ERAI_profile({1000:50})*86400, lev1({1000:50}), res)  ;

 res@gsnLeftString = "ERA5"
 plot(1) = gsn_csm_xy(wks, ERA5_profile({1000:50})*86400, lev4({1000:50}), res)  ;

 res@gsnLeftString = "NCEP"
 plot(2) = gsn_csm_xy(wks, NCEP_profile({1000:50})*86400, lev2({1000:50}), res)  ;

 res@gsnLeftString = "NCEP2"
 plot(3) = gsn_csm_xy(wks, NCEP2_profile({1000:50})*86400, lev3({1000:50}), res)  ;


pres = True
pres@gsnMaximize = True
gsn_panel(wks, plot, (/2,2/), pres)

end