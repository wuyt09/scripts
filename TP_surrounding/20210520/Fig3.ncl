
begin


yrStrt = 1979
yrLast = 2019
year = ispan(yrStrt, yrLast, 1)
nyear = dimsizes(year)
season = "JJA"

latstr_east = 30
latend_east = 45
lonstr_east = 105
lonend_east = 140

latstr_west = 30
latend_west = 55
lonstr_west = -10
lonend_west = 60

latstr_south = 10
latend_south = 28
lonstr_south = 70
lonend_south = 90

latstr_north = 45
latend_north = 60
lonstr_north = 85
lonend_north = 120

latstr_tp = 25 ;28
latend_tp = 40 ;38
lonstr_tp = 75 ;80
lonend_tp = 105

;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo ;(::-1,:)
printVarSummary(topo)

fa      = addfile ("/home/yangsong3/wuyt/data/lsm_1279l4_0.1x0.1.grb_v4_unpack.nc","r")
; lsdata = doubletofloat(a->lsm(0,:,:))
lsdata = fa->lsm(0,:,:)
printVarSummary(lsdata)
;=============================================================================================================
diri = "/home/yangsong3/data-observation/ERA5-monthly/surface/"
fili1 = "surface_thermal_downwards_mon_1x1_1979_2019.nc"
fili2 = "surface_thermal_downwards_clear_sky_mon_1x1_1979_2019.nc"
fili3 = "surface_net_thermal_mon_1x1_1979_2019.nc"
fili4 = "surface_net_thermal_clear_sky_mon_1x1_1979_2019.nc"
fili5 = "surface_solar_downwards_mon_1x1_1979_2019.nc"
fili6 = "surface_solar_downwards_clear_sky_mon_1x1_1979_2019.nc"
fili7 = "surface_net_solar_mon_1x1_1979_2019.nc"
fili8 = "surface_net_solar_clear_sky_mon_1x1_1979_2019.nc"

f1 = addfile(diri+fili1, "r")
Time = f1->time
YYYY = cd_calendar(Time, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

f2 = addfile(diri+fili2, "r")
f3 = addfile(diri+fili3, "r")
f4 = addfile(diri+fili4, "r")
f5 = addfile(diri+fili5, "r")
f6 = addfile(diri+fili6, "r")
f7 = addfile(diri+fili7, "r")
f8 = addfile(diri+fili8, "r")

STRD = short2flt(f1->strd(iYYYY,:,:))
STRD = STRD/86400.
STRDC = short2flt(f2->strdc(iYYYY,:,:))
STRDC = STRDC/86400.
STR  = short2flt(f3->str(iYYYY,:,:))
STR  = STR/86400.
STRC  = short2flt(f4->strc(iYYYY,:,:))
STRC  = STRC/86400.

SSRD = short2flt(f5->ssrd(iYYYY,:,:))
SSRD = SSRD/86400.
SSRDC = short2flt(f6->ssrdc(iYYYY,:,:))
SSRDC = SSRDC/86400.
SSR  = short2flt(f7->ssr(iYYYY,:,:))
SSR  = SSR/86400.
SSRC  = short2flt(f8->ssrc(iYYYY,:,:))
SSRC  = SSRC/86400.


STRU = STRD - STR
copy_VarCoords(STRD, STRU)

SSRU = SSRD - SSR
copy_VarCoords(SSRD, SSRU)
SSRUC = SSRDC - SSRC
copy_VarCoords(SSRDC, SSRUC)


STRD_CRE = STRD - STRDC
copy_VarCoords(STRD, STRD_CRE)
printVarSummary(STRD_CRE)
SSRD_CRE = SSRD - SSRDC
copy_VarCoords(SSRD, SSRD_CRE)
printVarSummary(SSRD_CRE)
SSRU_CRE = SSRU - SSRUC
copy_VarCoords(SSRU, SSRU_CRE)
printVarSummary(SSRU_CRE)

stru = month_to_season(STRU, season)
stru&time = year

strd = month_to_season(STRD, season)
strd&time = year
strd_cre = month_to_season(STRD_CRE, season)
strd_cre&time = year
printVarSummary(strd_cre)

ssrd = month_to_season(SSRD, season)
ssrd&time = year
ssrd_cre = month_to_season(SSRD_CRE, season)
ssrd_cre&time = year
printVarSummary(ssrd_cre)

ssru = month_to_season(SSRU, season)
ssru&time = year
ssru_cre = month_to_season(SSRU_CRE, season)
ssru_cre&time = year


lat = STRD&latitude
lon = STRD&longitude

lsm = linint2_Wrap(lsdata&longitude, lsdata&latitude, lsdata, True, lon, lat, 0)
printVarSummary(lsm)
delete(lsdata)
lsm1 = conform_dims(dimsizes(ssrd), lsm, (/1,2/))
stru = where(lsm1.ge.0.8, stru, stru@_FillValue)
strd = where(lsm1.ge.0.8, strd, strd@_FillValue)
strd_cre = where(lsm1.ge.0.8, strd_cre, strd_cre@_FillValue)
ssrd = where(lsm1.ge.0.8, ssrd, ssrd@_FillValue)
ssrd_cre = where(lsm1.ge.0.8, ssrd_cre, ssrd_cre@_FillValue)
ssru = where(lsm1.ge.0.8, ssru, ssru@_FillValue)
ssru_cre = where(lsm1.ge.0.8, ssru_cre, ssru_cre@_FillValue)
delete(lsm1)


obs_stru = runave_n_Wrap(stru, 10, 0, 0)
obs_strd = runave_n_Wrap(strd, 10, 0, 0)
obs_strd_cre = runave_n_Wrap(strd_cre, 10, 0, 0)
obs_ssrd = runave_n_Wrap(ssrd, 10, 0, 0)
obs_ssrd_cre = runave_n_Wrap(ssrd_cre, 10, 0, 0)
obs_ssru = runave_n_Wrap(ssru, 10, 0, 0)
obs_ssru_cre = runave_n_Wrap(ssru_cre, 10, 0, 0)
;================================================================================
nn = 31
nlat = dimsizes(lat)
nlon = dimsizes(lon)

diri2 = "/home/yangsong3/wuyt/sysu/scripts/TP_surrounding/cfram20210513_clearsky/"
nvname = (/"rad_conv","lwu","lwd","swu","swd","lwu_toa","lwd_toa","swu_toa","swd_toa"/)  ; ,"error"
nvars = dimsizes(nvname)


UNDEF = -999.
x = new((/nn,8,nvars,nlat,nlon/),float,UNDEF)
x!0 = "time"
x!1 = "process"
x!2 = "var"
x!3 = "lat"
x!4 = "lon"
x&var = nvname
x&lat = lat
x&lon = lon
x&lat@units = "degrees_north"
x&lon@units = "degrees_east"
x@long_name = "forcing"
x@units     = "W/m2"


do nnn = 0, nn-1
    n = nnn+1
    fili1 = diri2+"baseline_radsfc_ranc_1.grd"

    fili2 = diri2+"co2_radsfc_ranc_"+n+".grd"
    fili3 = diri2+"ta_radsfc_ranc_"+n+".grd"

    fili4 = diri2+"o3_radsfc_ranc_"+n+".grd"
    fili5 = diri2+"wv_radsfc_ranc_"+n+".grd"

    fili6 = diri2+"solar_radsfc_ranc_"+n+".grd"
    fili7 = diri2+"albedo_radsfc_ranc_"+n+".grd"

    fili8 = diri2+"warm_radsfc_ranc_"+n+".grd"

    do nv = 0, nvars-1
        x(nnn,0,nv,:,:) = fbindirread(fili1,nv,(/nlat,nlon/),"float")
        x(nnn,1,nv,:,:) = fbindirread(fili2,nv,(/nlat,nlon/),"float")
        x(nnn,2,nv,:,:) = fbindirread(fili3,nv,(/nlat,nlon/),"float")
        x(nnn,3,nv,:,:) = fbindirread(fili4,nv,(/nlat,nlon/),"float")
        x(nnn,4,nv,:,:) = fbindirread(fili5,nv,(/nlat,nlon/),"float")
        x(nnn,5,nv,:,:) = fbindirread(fili6,nv,(/nlat,nlon/),"float")
        x(nnn,6,nv,:,:) = fbindirread(fili7,nv,(/nlat,nlon/),"float")
        x(nnn,7,nv,:,:) = fbindirread(fili8,nv,(/nlat,nlon/),"float")
    end do
end do
printVarSummary(x)

; remove the error points
error = x(:,:,7,:,:) - dim_sum_n_Wrap(x(:,:,1:6,:,:), 2) + 3*x(:,:,0,:,:)
copy_VarCoords(x(:,:,0,:,:), error)
x(:,:,7,:,:) = where(abs(error).ge.5, x(:,:,7,:,:) - error, x(:,:,7,:,:))

; mask the sea area
lsm1 = conform_dims(dimsizes(x), lsm, (/3,4/))
x = where(lsm1.ge.0.8, x, x@_FillValue)

x_dR = x(:,:,2,:,:)
x_dS = x(:,:,4,:,:)  ;*(-1.) ; reverse the direction
x_uS = x(:,:,3,:,:)
printVarSummary(x_dS)

yrs = 1984
yre = 2014
year_new = ispan(yrs, yre, 1)
nyearn = dimsizes(year_new)

R_title = (/"CO2","Aerosol","WV","Air Temp","CRE"/)
S_title = (/"Solar","Aerosol","WV","Albedo","CRE"/)
;=================================================================
; (a) South region, downR and upS
;=================================================================
S_uR_obs = wgt_areaave_Wrap(obs_stru(:,{latstr_south:latend_south},{lonstr_south:lonend_south}), 1., 1., 0)
S_dR_obs = wgt_areaave_Wrap(obs_strd(:,{latstr_south:latend_south},{lonstr_south:lonend_south}), 1., 1., 0)
S_dR_obs_cre = wgt_areaave_Wrap(obs_strd_cre(:,{latstr_south:latend_south},{lonstr_south:lonend_south}), 1., 1., 0)
S_uS_obs = wgt_areaave_Wrap(obs_ssru(:,{latstr_south:latend_south},{lonstr_south:lonend_south}), 1., 1., 0)
S_uS_obs_cre = wgt_areaave_Wrap(obs_ssru_cre(:,{latstr_south:latend_south},{lonstr_south:lonend_south}), 1., 1., 0)


S_dR = wgt_areaave_Wrap(x_dR(:,:,{latstr_south:latend_south},{lonstr_south:lonend_south}), 1., 1., 0)
S_uS = wgt_areaave_Wrap(x_uS(:,:,{latstr_south:latend_south},{lonstr_south:lonend_south}), 1., 1., 0)

S_dR_co2 = S_dR(:,1) - S_dR(:,0)
S_dR_ta  = S_dR(:,2) - S_dR(:,0)
S_dR_o3  = S_dR(:,3) - S_dR(:,0)
S_dR_wv  = S_dR(:,4) - S_dR(:,0)
S_dR_cld = S_dR_obs_cre(5:35) - S_dR_obs_cre(4)

S_uS_o3  = S_uS(:,3) - S_uS(:,0)
S_uS_wv  = S_uS(:,4) - S_uS(:,0)
S_uS_sol = S_uS(:,5) - S_uS(:,0)
S_uS_al  = S_uS(:,6) - S_uS(:,0)
S_uS_cld = S_uS_obs_cre(5:35) - S_uS_obs_cre(4)

S_dR_total = S_dR_co2+S_dR_ta+S_dR_o3+S_dR_wv+S_dR_cld
copy_VarCoords(S_dR(:,0), S_dR_total)
S_uS_total = S_uS_o3+S_uS_wv+S_uS_sol+S_uS_al+S_uS_cld
copy_VarCoords(S_uS(:,0), S_uS_total)

S_dR_total_obs = S_dR_obs(5:35) - S_dR_obs(4)
S_uS_total_obs = S_uS_obs(5:35) - S_uS_obs(4)

S_dR_areo = S_dR_total_obs - S_dR_total
S_uS_areo = S_uS_total_obs - S_uS_total


title = (/"CRE","WV","Aerosol","Air Temp.","Albedo","CO2","Solar"/)
lines = new((/7,nyearn/), float)
lines(0,:) = S_dR_cld - S_uS_cld
lines(1,:) = S_dR_wv - S_uS_wv
lines(2,:) = S_dR_areo - S_uS_areo
lines(3,:) = S_dR_ta
lines(4,:) = (-1.) * S_uS_al
lines(5,:) = S_dR_co2
lines(6,:) = (-1.) * S_uS_sol

;==================================================================================================
zero = new(nyearn, float)
zero = 0.

y = new((/2,nyearn/), float)
y(0,:) = S_uR_obs(5:35) - S_uR_obs(4) ;S_dR_total_obs - S_uS_total_obs
y(1,:) = dim_sum_n_Wrap(lines, 0)

d1 = lines(0,:)
d2 = lines(1,:)
d3 = lines(2,:)
d4 = lines(3,:)
d5 = lines(4,:)
d6 = lines(5,:)
d7 = lines(6,:)
; d8 = x_final(:,7)
; d9 = x_final(:,8)
; d0 = x_final(:,9)

;==================================================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_surrounding/20210520/"
figname = "Fig3"

delete(year)
year = year_new

wks = gsn_open_wks("pdf", figpath+figname)
; nplot = 1
; plot = new(nplot, graphic)

res = True
res@gsnMaximize = True
res@gsnDraw  = False
res@gsnFrame = False

res@vpWidthF = 0.8
res@vpHeightF = 0.4

;---Set axes limits. Add extra space for X max.
res@trYMinF                  = -4
res@trYMaxF                  = 8;(max(d2) + max(d1)+max(d3) + max(d4))*1.2  
res@trXMinF                  = yrStrt+4
res@trXMaxF                  = yrLast-4

;---Put city labels on X axis
; res@tmXBMode                 = "Explicit"
; res@tmXBValues               = ispan(1,ncities,1)
; res@tmXBLabels               = cities
res@tmXBLabelJust            = "BottomCenter"
res@tmXBLabelFontHeightF     = 0.016


res@gsnScale                 = True    ; Make labels on axes same size

;--Turn off X axis tickmarks
res@tmXBMajorLengthF         = 0.01
res@tmXBMajorOutwardLengthF  = 0.  
res@tmYLMajorLengthF         = 0.01
res@tmYLMajorOutwardLengthF  = 0.01

res@tmYROn                   = False
res@tiYAxisOn                = True

res@tmXTOn                   = True   ; Turn off top tickmarks
res@tmXTMajorLengthF         = 0.01
res@tmXTMajorOutwardLengthF  = 0.  

;---Create blank plot without Y grid lines
; plot_without_xgrid = gsn_csm_blank_plot(wks,res)
res@gsnLeftString = "(e) South region: ("+lonstr_south+"-"+lonend_south+"E,"+latstr_south+"-"+latend_south+"N)"
res@gsnRightString = " "
res@tiXAxisString = "Year"
res@tiYAxisFontHeightF = 0.015
res@tiYAxisString = "W/m~S~2~N~"

res@gsnYRefLine = 0.0
res@gsnYRefLineThicknessF = 3.0
res@xyLineColor = "black"
res@xyLineThicknessF = 3.0
plot_without_xgrid = gsn_csm_xy(wks, year, zero, res)


res@xyMonoDashPattern = False
res@xyDashPatterns = (/0,16/)
res@xyMonoLineColor = False
res@xyLineColors = (/"red","black"/)
res@xyMonoLineThickness = False
res@xyLineThicknesses = (/8.0,5.0/)

plot_without_xgrid_xy = gsn_csm_xy(wks, year, y, res)


;---------------------------------------------------------------
;---Arrays to hold the bar values and the bars themselves.
ncities = nn

dum1 = new(ncities,graphic)
dum2 = new(ncities,graphic)
dum3 = new(ncities,graphic)
dum4 = new(ncities,graphic)
dum5 = new(ncities,graphic)
dum6 = new(ncities,graphic)
dum7 = new(ncities,graphic)
dum8 = new(ncities,graphic)
; dum9 = new(ncities,graphic)
; dum0 = new(ncities,graphic)

;---Set some resources for the bars.
bres           = True
bres@gsEdgesOn = True      ; Outline the polygons (bars)
;
; Loop through each city and attach bar to plot that
; doesn't have grid lines.

do i=1,ncities
  ybar_p_min = 0.
  ybar_p_max = 0.

  ybar_n_min = 0.
  ybar_n_max = 0.

  ybar_positive = (/0.,0.,0.,0.,0./)
  ybar_negative = (/0.,0.,0.,0.,0./)

;*************************************************
  if (d1(i-1).ge.0) then
    ybar_p_max = ybar_p_max + d1(i-1)
    ybar_positive = (/ybar_p_min,ybar_p_max,ybar_p_max,ybar_p_min,ybar_p_min/)
    ybar = ybar_positive
    ybar_p_min = ybar_p_max
  end if

  if (d1(i-1).lt.0) then
    ybar_n_max = ybar_n_max + d1(i-1)
    ybar_negative = (/ybar_n_min,ybar_n_max,ybar_n_max,ybar_n_min,ybar_n_min/)
    ybar = ybar_negative
    ybar_n_min = ybar_n_max
  end if
  ; bres@gsFillColor           = "darkorange2"
  ; bres@gsEdgeColor           = "darkorange2"
  bres@gsFillColor           = "coral2"
  bres@gsEdgeColor           = "coral2"
  bres@gsFillIndex           = 0
  bres@gsEdgeThicknessF     = 2.5

  ; ybar      = (/0,d1(i-1),d1(i-1),0,0/)
  j = 1983+i  
  xbar      = (/j-0.25,j-0.25,j+0.25,j+0.25,j-0.25/)
  dum1(i-1) = gsn_add_polygon(wks,plot_without_xgrid,xbar,ybar,bres)

;*************************************************
  if (d2(i-1).ge.0) then
    ybar_p_max = ybar_p_max + d2(i-1)
    ybar_positive = (/ybar_p_min,ybar_p_max,ybar_p_max,ybar_p_min,ybar_p_min/)
    ybar = ybar_positive
    ybar_p_min = ybar_p_max
  end if

  if (d2(i-1).lt.0) then
    ybar_n_max = ybar_n_max + d2(i-1)
    ybar_negative = (/ybar_n_min,ybar_n_max,ybar_n_max,ybar_n_min,ybar_n_min/)
    ybar = ybar_negative
    ybar_n_min = ybar_n_max
  end if

  ; bres@gsFillColor          = "gold1"
  ; bres@gsEdgeColor          = "gold1"
  bres@gsFillColor          = "goldenrod2"
  bres@gsEdgeColor          = "goldenrod2"
  bres@gsFillIndex           = 0
  bres@gsEdgeThicknessF     = 2.5
  
  ; ybar      = (/d1(i-1),d2(i-1)+d1(i-1),d2(i-1)+d1(i-1),d1(i-1),d1(i-1)/)
  dum2(i-1) = gsn_add_polygon(wks,plot_without_xgrid,xbar,ybar,bres)
  
;*************************************************
  if (d3(i-1).ge.0) then
    ybar_p_max = ybar_p_max + d3(i-1)
    ybar_positive = (/ybar_p_min,ybar_p_max,ybar_p_max,ybar_p_min,ybar_p_min/)
    ybar = ybar_positive
    ybar_p_min = ybar_p_max
  end if

  if (d3(i-1).lt.0) then
    ybar_n_max = ybar_n_max + d3(i-1)
    ybar_negative = (/ybar_n_min,ybar_n_max,ybar_n_max,ybar_n_min,ybar_n_min/)
    ybar = ybar_negative
    ybar_n_min = ybar_n_max
  end if

  bres@gsFillColor          = "darkviolet"
  bres@gsEdgeColor          = "darkviolet"
  bres@gsFillIndex          = 16
  bres@gsFillScaleF         = 0.8
  bres@gsEdgeThicknessF     = 2.5

  ; ybar      = (/d1(i-1)+d2(i-1),d3(i-1)+d2(i-1)+d1(i-1),\
  ;               d3(i-1)+d2(i-1)+d1(i-1),d1(i-1)+d2(i-1),\
  ;               d1(i-1)+d2(i-1)/)
  dum3(i-1) = gsn_add_polygon(wks,plot_without_xgrid,xbar,ybar,bres)

;*************************************************
  if (d4(i-1).ge.0) then
    ybar_p_max = ybar_p_max + d4(i-1)
    ybar_positive = (/ybar_p_min,ybar_p_max,ybar_p_max,ybar_p_min,ybar_p_min/)
    ybar = ybar_positive
    ybar_p_min = ybar_p_max
  end if

  if (d4(i-1).lt.0) then
    ybar_n_max = ybar_n_max + d4(i-1)
    ybar_negative = (/ybar_n_min,ybar_n_max,ybar_n_max,ybar_n_min,ybar_n_min/)
    ybar = ybar_negative
    ybar_n_min = ybar_n_max
  end if

  bres@gsFillColor          = "lightskyblue3"
  bres@gsEdgeColor          = "lightskyblue3"
  bres@gsFillIndex          = 1
  bres@gsEdgeThicknessF     = 2.5

  ; ybar      = (/d1(i-1)+d2(i-1)+d3(i-1),d4(i-1)+d3(i-1)+d2(i-1)+d1(i-1),\
  ;               d4(i-1)+d3(i-1)+d2(i-1)+d1(i-1),d1(i-1)+d2(i-1)+d3(i-1),\
  ;               d1(i-1)+d2(i-1)+d3(i-1)/)
  dum4(i-1) = gsn_add_polygon(wks,plot_without_xgrid,xbar,ybar,bres)

;*************************************************
  if (d5(i-1).ge.0) then
    ybar_p_max = ybar_p_max + d5(i-1)
    ybar_positive = (/ybar_p_min,ybar_p_max,ybar_p_max,ybar_p_min,ybar_p_min/)
    ybar = ybar_positive
    ybar_p_min = ybar_p_max
  end if

  if (d5(i-1).lt.0) then
    ybar_n_max = ybar_n_max + d5(i-1)
    ybar_negative = (/ybar_n_min,ybar_n_max,ybar_n_max,ybar_n_min,ybar_n_min/)
    ybar = ybar_negative
    ybar_n_min = ybar_n_max
  end if

  bres@gsFillColor          = "green3"
  bres@gsEdgeColor          = "green3"
  bres@gsFillIndex          = 0
  bres@gsEdgeThicknessF     = 2.5

  dum5(i-1) = gsn_add_polygon(wks,plot_without_xgrid,xbar,ybar,bres)

;*************************************************
  if (d6(i-1).ge.0) then
    ybar_p_max = ybar_p_max + d6(i-1)
    ybar_positive = (/ybar_p_min,ybar_p_max,ybar_p_max,ybar_p_min,ybar_p_min/)
    ybar = ybar_positive
    ybar_p_min = ybar_p_max
  end if

  if (d6(i-1).lt.0) then
    ybar_n_max = ybar_n_max + d6(i-1)
    ybar_negative = (/ybar_n_min,ybar_n_max,ybar_n_max,ybar_n_min,ybar_n_min/)
    ybar = ybar_negative
    ybar_n_min = ybar_n_max
  end if

  bres@gsFillColor          = "dodgerblue"
  bres@gsEdgeColor          = "dodgerblue"
  bres@gsFillIndex          = 0
  bres@gsEdgeThicknessF     = 2.5

  dum6(i-1) = gsn_add_polygon(wks,plot_without_xgrid,xbar,ybar,bres)

;*************************************************
  if (d7(i-1).ge.0) then
    ybar_p_max = ybar_p_max + d7(i-1)
    ybar_positive = (/ybar_p_min,ybar_p_max,ybar_p_max,ybar_p_min,ybar_p_min/)
    ybar = ybar_positive
    ybar_p_min = ybar_p_max
  end if

  if (d7(i-1).lt.0) then
    ybar_n_max = ybar_n_max + d7(i-1)
    ybar_negative = (/ybar_n_min,ybar_n_max,ybar_n_max,ybar_n_min,ybar_n_min/)
    ybar = ybar_negative
    ybar_n_min = ybar_n_max
  end if

  bres@gsFillColor          = "wheat3"
  bres@gsEdgeColor          = "wheat3"
  bres@gsFillIndex          = 17
  bres@gsEdgeThicknessF     = 2.5

  dum7(i-1) = gsn_add_polygon(wks,plot_without_xgrid,xbar,ybar,bres)

; ;*************************************************
;   if (d8(i-1).ge.0) then
;     ybar_p_max = ybar_p_max + d8(i-1)
;     ybar_positive = (/ybar_p_min,ybar_p_max,ybar_p_max,ybar_p_min,ybar_p_min/)
;     ybar = ybar_positive
;     ybar_p_min = ybar_p_max
;   end if

;   if (d8(i-1).lt.0) then
;     ybar_n_max = ybar_n_max + d8(i-1)
;     ybar_negative = (/ybar_n_min,ybar_n_max,ybar_n_max,ybar_n_min,ybar_n_min/)
;     ybar = ybar_negative
;     ybar_n_min = ybar_n_max
;   end if

;   bres@gsFillColor          = "plum1"
;   bres@gsEdgeColor          = "plum1"
;   bres@gsFillIndex          = 8
;   bres@gsEdgeThicknessF     = 2.5

;   dum8(i-1) = gsn_add_polygon(wks,plot_without_xgrid,xbar,ybar,bres)


;*************************************************
  delete([/ybar_p_min, ybar_p_max, ybar_n_min, ybar_n_max,\
           ybar_positive, ybar_negative, ybar/])

end do



; Create some labelbars for a legend
lb_colors = (/"coral2","goldenrod2","darkviolet","lightskyblue3","green3","dodgerblue","wheat3","plum1"/)
lb_dash = (/0,0,16,1,0,0,17,8/)

lb_x = new((/10,5/), "float")
; lb_x(0,:) = (/0.5,0.5,2.,2.,0.5/)
lb_x(0,:) = (/1984.5,1984.5,1986.,1986.,1984.5/)
do i = 1, 9
  lb_x(i,:) = lb_x(i-1,:) + 3.0
end do

lb_y = (/-2.8, -2.6, -2.6, -2.8, -2.8/)

dum_lb = new(10, graphic)
lb_res = True
lb_res@gsEdgesOn = True
lb_res@gsEdgeThicknessF = 2.5
do i = 0, 6
  lb_res@gsFillColor = lb_colors(i)
  lb_res@gsEdgeColor = lb_colors(i)
  lb_res@gsFillIndex = lb_dash(i)
  dum_lb(i) = gsn_add_polygon(wks, plot_without_xgrid, lb_x(i,:), lb_y, lb_res)
end do

text = (/"CRE","WV","Aero","Ta","Alo","CO2","Solar","ERA5 upR","cfram total"/) ;(/"AL","WV","CLD","CO2","ATM","HS","SH","LH","Sum","ERA5"/)
text_res = True
text_res@txFontHeightF = 0.012
text_res@txJust             = "CenterLeft"
text_plot = new(10, graphic)
text_x = new(10, "float")
text_x(0) = 1986.1
do i = 1, 6
  text_x(i) = text_x(i-1) + 3.0
end do
do i = 7, 8
  text_x(i) = text_x(i-1) + 2.8
end do
text_y = -2.7
do i = 0, 8
  text_plot(i) = gsn_add_text(wks, plot_without_xgrid, text(i), text_x(i), text_y, text_res)
end do

line_color = (/"red","black"/)
line_dash = (/0,16/)
line_th = (/8.0,5.0/)
line_res = True
line_plot = new(2, graphic)
line_x = new((/2,2/), "float")
line_x(0,:) = (/text_x(8)-1.6,text_x(8)-0.2/)
line_x(1,:) = (/text_x(9)-1.6,text_x(9)-0.2/)
line_y = (/-2.7,-2.7/)
do i = 0, 1
  line_res@gsLineColor = line_color(i)
  line_res@gsLineDashPattern = line_dash(i)
  line_res@gsLineThicknessF = line_th(i)
  line_plot(i) = gsn_add_polyline(wks, plot_without_xgrid, line_x(i,:), line_y, line_res)
end do

draw(plot_without_xgrid)
draw(plot_without_xgrid_xy)
frame(wks)

end