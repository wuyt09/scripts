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


;=================================================================================================

diriq1 = "/home/yangsong3/wuyt/data/ERA-Interim/monthly/"
filiq1 = "q1.mon.mean.1979.2016.ERAI.nc"
dirip1 = "/home/yangsong3/data-observation/wuyt/ERAI/monthly/1.5x1.5/"
filip1 = "sp_era_mn_197901_201507.nc"

fq1 = addfile(diriq1+filiq1, "r")
TIME = fq1->time
YYYY = cd_calendar(TIME, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
fp1 = addfile(dirip1+filip1, "r")

Q1 = short2flt(fq1->q1(iYYYY,:,::-1,:))
Q1 = Q1*1004.7/9.8
P1 = short2flt(fp1->sp(iYYYY,::-1,:))
printVarSummary(Q1)
lev1 = Q1&lev
lev1 = lev1*100
lev1!0 = "level"
lev1&level = lev1
dp1 = dpres_plevel_Wrap(lev1, P1, min(lev1), 0)
printVarSummary(dp1)
Q1dp = Q1*dp1
copy_VarCoords(Q1, Q1dp)
T1 = dim_sum_n_Wrap(Q1dp, 1)
t1_med = month_to_season(T1, season)
t1 = dim_avg_n_Wrap(t1_med, 0)
printVarSummary(t1)

;==================================================================================================
diriq4 = "/home/yangsong3/data-observation/ERA5-monthly/pressure/thermal_equation/"
filiq4 = "q1.monmean.1979.2019.nc"
dirip4 = "/home/yangsong3/data-observation/ERA5-monthly/surface/"
filip4 = "sur_pressure_mon_1x1_1979_2019.nc"
fq4 = addfile(diriq4+filiq4, "r")
fp4 = addfile(dirip4+filip4, "r")
TIME4 = fq4->time
YYYY4 = cd_calendar(TIME4, -1)/100
iYYYY4 = ind(YYYY4.ge.yrStrt .and. YYYY4.le.yrLast)

Q4 = short2flt(fq4->q1(iYYYY4,:,::-1,:))
Q4 = Q4*1004.7/9.8
P4 = short2flt(fp4->sp(iYYYY4,::-1,:))

lev4 = Q4&lev
lev4 = lev4*100
lev4!0 = "level"
lev4&level = lev4
dp4 = dpres_plevel_Wrap(lev4, P4, min(lev4), 0)
Q4dp = Q4*dp4
copy_VarCoords(Q4, Q4dp)
T4 = dim_sum_n_Wrap(Q4dp, 1)

t4_med = month_to_season(T4, season)
t4 = dim_avg_n_Wrap(t4_med, 0)

;==================================================================================================
diriq2 = "/home/yangsong3/wuyt/data/"
filiq2 = "q1.mon.mean.1979.2013.nc"
dirip2 = "/home/yangsong3/data-observation/NCEP_DOE/monthly/"
filip2 = "pres.mon.mean.nc"

fq2 = addfile(diriq2+filiq2, "r")
fp2 = addfile(dirip2+filip2, "r")
TIME2 = fq2->time
YYYY2 = cd_calendar(TIME2, -1)/100
iYYYY2 = ind(YYYY2.ge.yrStrt .and. YYYY2.le.yrLast)

Q2 = short2flt(fq2->q1(iYYYY2,:,::-1,:))
Q2 = Q2*1004.7/9.8
P2 = short2flt(fp2->pres(iYYYY2,::-1,:))
P2 = P2*100

lev2 = Q2&lev
lev2 = lev2*100
lev2!0 = "level"
lev2&level = lev2

dp2 = dpres_plevel_Wrap(lev2, P2, min(lev2), 0)
Q2dp = Q2*dp2
copy_VarCoords(Q2, Q2dp)
T2 = dim_sum_n_Wrap(Q2dp, 1)

t2_med = month_to_season(T2, season)
t2 = dim_avg_n_Wrap(t2_med, 0)
printVarSummary(t2)

;==================================================================================================
diriq3 = "/home/yangsong3/wangzq/DATA/jra/"
filiq3 = "Q1I.jra55.mon.mean.nc"
fq3 = addfile(diriq3+filiq3, "r")
TIME3 = fq3->time
YYYY3 = cd_calendar(TIME3, -1)/100
iYYYY3 = ind(YYYY3.ge.yrStrt .and. YYYY3.le.yrLast)

T3 = short2flt(fq3->Q1(iYYYY3,:,:))

t3_med = month_to_season(T3, season)
t3 = dim_avg_n_Wrap(t3_med, 0)
printVarSummary(t3)


;==================================================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_surf_temp_climate_effects/20210708_obs/Q1/"
figname = "check_spatial_distribution_JJA_Q1_multi_source"

wks = gsn_open_wks("pdf", figpath+figname)

plot = new(4, graphic)
plot_topo = plot

;----------------------------------------------------
 res_topo            = True
 res_topo@cnFillOn   = False
 res_topo@cnLinesOn  = True
 res_topo@gsnFrame   = False
 res_topo@gsnDraw    = False
 res_topo@cnLineLabelsOn         = False
 res_topo@cnLineThicknessF       = 4
 res_topo@cnLineDashPattern      = 0
 res_topo@cnLineColor            = "green4"
 res_topo@cnInfoLabelOn          = False
 res_topo@cnLevelSelectionMode   ="ExplicitLevels"
 res_topo@cnLevels               = 2000 ;(/0.9/)
 res_topo@gsnLeftString = ""
 res_topo@gsnRightString = ""
;----------------------------------------------------

res = True
res@gsnDraw  = False
res@gsnFrame = False
res@gsnMaximize = True
res@gsnAddCyclic = False

res@mpOutlineOn = True
res@mpGeophysicalLineThicknessF = 2
res@mpGeophysicalLineColor = "grey25" ;"red" 
res@mpFillDrawOrder = "PreDraw"
res@mpFillOn = True

res@mpMinLatF = latstr
res@mpMaxLatF = latend
res@mpMinLonF = lonstr
res@mpMaxLonF = lonend

res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False

res@gsnStringFontHeightF = 0.035
res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025

; res@gsnRightString = "mm/day"

res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = 0
res@cnMaxLevelValF = 300
res@cnLevelSpacingF = 30
res@lbLabelBarOn = False

res@gsnLeftString = "ERAI"
plot(0) = gsn_csm_contour_map(wks, t1, res)
plot_topo(0) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(0), plot_topo(0))

res@gsnLeftString = "ERA5"
plot(1) = gsn_csm_contour_map(wks, t4, res)
plot_topo(1) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(1), plot_topo(1))

res@gsnLeftString = "NCEP2"
plot(2) = gsn_csm_contour_map(wks, t2, res)
plot_topo(2) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(2), plot_topo(2))

res@gsnLeftString = "JRA55"
plot(3) = gsn_csm_contour_map(wks, t3, res)
plot_topo(3) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(3), plot_topo(3))



pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
gsn_panel(wks, plot, (/2,2/), pres)

end