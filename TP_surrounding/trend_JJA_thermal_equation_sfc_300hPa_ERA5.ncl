
begin

lonstr = 0
lonend = 130 ;360
latstr = 10 ;10
latend = 70

yrStrt = 1979
yrLast = 2019
year = ispan(yrStrt, yrLast, 1)
nyear = dimsizes(year)
season = "JJA"

LEVs = 1000
LEVe = 300

latstr_east = 35
latend_east = 55
lonstr_east = 84 ;80
lonend_east = 115 ;120

latstr_west = 40
latend_west = 55
lonstr_west = 10
lonend_west = 50

latstr_south = 15
latend_south = 25
lonstr_south = 60
lonend_south = 110

latstr_north = 45
latend_north = 60
lonstr_north = 70
lonend_north = 110

latstr_tp = 25 ;28
latend_tp = 40 ;38
lonstr_tp = 75 ;80
lonend_tp = 105

region = (/"West","East","South","North","TP"/)
;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)
printVarSummary(topo)


;===================================================================
diri = "/home/yangsong3/data-observation/ERA5-monthly/pressure/thermal_equation/"
fili1 = "adv_T.monmean.1979.2019.nc"
fili2 = "dTdp.monmean.1979.2019.nc"
fili3 = "dTdt.monmean.1979.2019.nc"
fili4 = "q1.monmean.1979.2019.nc"

f1 = addfile(diri+fili1, "r")
f2 = addfile(diri+fili2, "r")
f3 = addfile(diri+fili3, "r")
f4 = addfile(diri+fili4, "r")

Time = f1->time
YYYY = cd_calendar(Time, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

x1 = short2flt(f1->adv_T(iYYYY,::-1,::-1,:))
x2 = short2flt(f2->dTdp(iYYYY,::-1,::-1,:))
x3 = short2flt(f3->dTdt(iYYYY,::-1,::-1,:))
x4 = short2flt(f4->q1(iYYYY,::-1,::-1,:))
printVarSummary(x1)

lat = x1&lat
lon = x1&lon
nlat = dimsizes(lat)
nlon = dimsizes(lon)

filips = "/home/yangsong3/data-observation/ERA5-monthly/surface/sur_pressure_mon_1x1_1979_2019.nc"
fps = addfile(filips, "r")
ps = short2flt(fps->sp(iYYYY,::-1,:))
printVarSummary(ps)


plev = tofloat(x1&lev)
plev = plev*100
plev!0 = "lev"
plev&lev = x1&lev
plev@units = "Pa"
; print(plev)

dp = dpres_plevel_Wrap(plev, ps, min(plev), 0) ;LEVs(i)*100
; dP = dpres_plevel_Wrap(plev, max(plev), min(plev), 0)
; dp = conform_dims(dimsizes(t), dP, 1)
; copy_VarCoords(t, dp)
printVarSummary(dp)
; print(dp(90,:,{35},{90}))
; exit

x1dp = x1*dp
x2dp = x2*dp
x3dp = x3*dp
x4dp = x4*dp
copy_VarCoords(x1, x1dp)
copy_VarCoords(x2, x2dp)
copy_VarCoords(x3, x3dp)
copy_VarCoords(x4, x4dp)
printVarSummary(x1dp)


x1_tropo = dim_sum_n_Wrap(x1dp(:,{LEVs:LEVe},:,:), 1)/dim_sum_n_Wrap(dp(:,{LEVs:LEVe},:,:), 1)
copy_VarCoords(x1(:,0,:,:), x1_tropo)
x2_tropo = dim_sum_n_Wrap(x2dp(:,{LEVs:LEVe},:,:), 1)/dim_sum_n_Wrap(dp(:,{LEVs:LEVe},:,:), 1)
copy_VarCoords(x2(:,0,:,:), x2_tropo)
x3_tropo = dim_sum_n_Wrap(x3dp(:,{LEVs:LEVe},:,:), 1)/dim_sum_n_Wrap(dp(:,{LEVs:LEVe},:,:), 1)
copy_VarCoords(x3(:,0,:,:), x3_tropo)
x4_tropo = dim_sum_n_Wrap(x4dp(:,{LEVs:LEVe},:,:), 1)/dim_sum_n_Wrap(dp(:,{LEVs:LEVe},:,:), 1)
copy_VarCoords(x4(:,0,:,:), x4_tropo)

x1_clm = month_to_season(x1_tropo, season)
x2_clm = month_to_season(x2_tropo, season)
x3_clm = month_to_season(x3_tropo, season)
x4_clm = month_to_season(x4_tropo, season)
printVarSummary(x1_clm)

x = new((/4,nyear,nlat,nlon/), float)
x(0,:,:,:) = x1_clm
x(1,:,:,:) = x2_clm
x(2,:,:,:) = x3_clm
x(3,:,:,:) = x4_clm
copy_VarCoords(x4_clm, x(0,:,:,:))
term = (/"adv_T","dTdp","dTdt","Q1"/)
x!0 = "term"
x&term = term
x = x*86400
x@units = "k/day"

rc = x(:,0,:,:)
rc = 0
prob = rc
do i = 0, 3
    rc_med = regCoef_n(year, x(i,:,:,:), 0, 0)
    rc(i,:,:) = rc_med

    tval = onedtond(rc_med@tval, dimsizes(rc_med))
    df   = onedtond(rc_med@nptxy, dimsizes(rc_med))
    b    = tval
    b    = 0.5
    prob(i,:,:) = betainc(df/(df+tval^2), df/2.0, b)
end do
copy_VarCoords(x1_clm(0,:,:), rc(0,:,:))
copy_VarCoords(rc, prob)

rc@long_name = "regression coefficient"
prob@long_name = "probability"


trend = rc*10
copy_VarCoords(rc, trend)
trend@units = "K/day/decade"
printVarSummary(trend)

sig = 0.05
;===============================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_surrounding/"
figname = "trend_"+season+"_thermal_equation_sfc_300hPa_ERA5"
wks = gsn_open_wks("pdf", figpath+figname)

plot = new(4, graphic)
plot_topo = plot
plot_prob = new(4, graphic)
dum1 = plot
dum2 = plot
dum3 = plot
dum4 = plot

 res_topo            = True
 res_topo@cnFillOn   = False
 res_topo@cnLinesOn  = True
 res_topo@gsnFrame   = False
 res_topo@gsnDraw    = False
 res_topo@cnLineLabelsOn         = False
 res_topo@cnLineThicknessF       = 5
 res_topo@cnLineColor            = "forestgreen"
 res_topo@cnInfoLabelOn          = False
 res_topo@cnLevelSelectionMode   ="ExplicitLevels"
 res_topo@cnLevels               = (/2000/)
 res_topo@gsnLeftString = ""
 res_topo@gsnRightString = ""


;----------------------------------------------
; draw prob
;----------------------------------------------
 res3  = True
 res3@gsnDraw   = False
 res3@gsnFrame  = False

 ;res3@cnFillOn  = False 
 res3@cnLinesOn = False
 res3@cnLineLabelsOn = False
 res3@cnInfoLabelOn  = False
 ;res3@lbLabelBarOn   = False

 res3@gsnRightString = ""

 res3@cnLevelSelectionMode = "ManualLevels"
 res3@cnMinLevelValF = 0.00
 res3@cnMaxLevelValF = 0.10
 res3@cnLevelSpacingF = 0.01

 ; res3@cnFillPattern = 17 ;(/17/)
 res3@cnFillColor  =  "black"
 res3@gsnLeftString = ""
 res3@cnFillDotSizeF = 0.003
 res3@cnFillScaleF  = 1.2
 ;----------------------------------------------
 ;----------------------------------------------

res = True
res@gsnDraw = False
res@gsnFrame = False
res@gsnMaximize = True

res@mpOutlineOn = True
res@mpFillOn = True
res@mpGeophysicalLineThicknessF = 5
res@mpGeophysicalLineColor = "grey25" ;"red" 
; res@mpFillDrawOrder = "PostDraw"
; res@mpCenterLonF = 180
res@mpMinLatF = latstr
res@mpMaxLatF = latend
res@mpMinLonF = lonstr
res@mpMaxLonF = lonend
res@mpCenterLonF = (lonstr+lonend)/2


res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False
; res@cnFillMode = "CellFill"

res@gsnStringFontHeightF = 0.025
res@tmXBLabelFontHeightF = 0.02
res@tmYLLabelFontHeightF = 0.02

res@cnFillPalette = "temp_diff_18lev"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -0.5
res@cnMaxLevelValF = 0.5
res@cnLevelSpacingF = 0.05
res@lbLabelBarOn = False

res@pmLabelBarHeightF = 0.1
res@pmLabelBarWidthF = 0.7
res@pmLabelBarOrthogonalPosF = 0.1


;-----------------------------------
gres = True
gres@gsLineColor = "blue"
gres@gsLineThicknessF = 5.0
gres@gsLineDashPattern = 16

 glon1 = (/lonend_north,lonend_north,lonstr_north,lonstr_north,lonend_north/)
 glat1 = (/latstr_north,latend_north,latend_north,latstr_north,latstr_north/) ; north
 glon2 = (/lonend_east,lonend_east,lonstr_east,lonstr_east,lonend_east/)
 glat2 = (/latstr_east,latend_east,latend_east,latstr_east,latstr_east/) ; east
 glon3 = (/lonend_south,lonend_south,lonstr_south,lonstr_south,lonend_south/)
 glat3 = (/latstr_south,latend_south,latend_south,latstr_south,latstr_south/) ; south
 glon4 = (/lonend_west,lonend_west,lonstr_west,lonstr_west,lonend_west/)
 glat4 = (/latstr_west,latend_west,latend_west,latstr_west,latstr_west/) ; west
;-----------------------------------

do i = 0, 3
    res@gsnLeftString = term(i)
    plot(i) = gsn_csm_contour_map(wks, trend(i,:,:), res)
    plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(i), plot_topo(i))

    plot_prob(i) = gsn_csm_contour(wks,prob(i,:,:),res3)
    plot_prob(i) = ShadeLtContour(plot_prob(i),sig,17)
    overlay(plot(i), plot_prob(i))

    ; dum1 = gsn_add_polyline(wks, plot(0), glon1, glat1, gres)
    dum2(i) = gsn_add_polyline(wks, plot(i), glon2, glat2, gres)
    dum3(i) = gsn_add_polyline(wks, plot(i), glon3, glat3, gres)
    dum4(i) = gsn_add_polyline(wks, plot(i), glon4, glat4, gres)
end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
pres@lbLabelFontHeightF = 0.015
gsn_panel(wks, plot, (/2,2/), pres)

end