
begin

yrStrt = 1979
yrLast = 2019
year = ispan(yrStrt, yrLast, 1)
nyear = dimsizes(year)
season = "JJA"

LEVs = 1000
LEVe = 300

; latstr_east = 35
; latend_east = 55
; lonstr_east = 80
; lonend_east = 120
latstr_east = 35
latend_east = 55
lonstr_east = 84 ;80
lonend_east = 115 ;120

;=============================================================================================================
diri = "/home/yangsong3/data-observation/ERA5-monthly/pressure/"
fili = "t_mon_1x1_1979_2019.nc"

f = addfile(diri+fili, "r")

Time = f->time
YYYY = cd_calendar(Time, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

t = short2flt(f->t(iYYYY,:,:,:))
printVarSummary(t)


filips = "/home/yangsong3/data-observation/ERA5-monthly/surface/sur_pressure_mon_1x1_1979_2019.nc"
fps = addfile(filips, "r")
ps = short2flt(fps->sp(iYYYY,:,:))
; printVarSummary(ps)

plev = tofloat(t&level)
plev = plev*100
plev!0 = "lev"
plev&lev = t&level
plev@units = "Pa"
; print(plev)

dp = dpres_plevel_Wrap(plev, ps, min(plev), 0)
; printVarSummary(dp)

tdp = t*dp
copy_VarCoords(t, tdp)
; printVarSummary(tdp)

t_tropo = dim_sum_n_Wrap(tdp(:,{LEVs:LEVe},:,:), 1)/dim_sum_n_Wrap(dp(:,{LEVs:LEVe},:,:), 1)
copy_VarCoords(t(:,0,:,:), t_tropo)

t_tropo_JJA = month_to_season(t_tropo, season)
t_series = wgt_areaave_Wrap(t_tropo_JJA(:,{latstr_east:latend_east},{lonstr_east:lonend_east}), 1., 1., 0)
; printVarSummary(t_series)
print(t_series)

t_base = avg(t_series(0:9))
t_new = t_series - t_base
print(t_new)
t_runmean = runave(t_new, 10, 0)
t_runmean!0 = "year"
t_runmean&year = year
print(t_runmean)
;=============================================================================================================

nvname = (/"albedo","wv","cloud","cloud_sw","cloud_lw","co2","o3","solar",\
           "dyn","atm_dyn","ocean_dyn","shflx","lhflx"/)  ; ,"error"
nn = 31

lat  = f->latitude ;(::-1)
lon  = f->longitude
level  = f->level
nlev = 38
lev = new(nlev, integer, "No_FillValue")
lev(0:nlev-2) = level
lev(nlev-1) = 1001
print(lev)

delete(f)

nlat = dimsizes(lat)
nlon = dimsizes(lon)
nlev = dimsizes(lev)
nvars = dimsizes(nvname)

setfileoption("bin","ReadByteOrder","Native")

UNDEF = -999.
x = new((/nn,nvars,nlev,nlat,nlon/),float,UNDEF)
x!0 = "series"
x!1 = "var"
x!2 = "lev"
x!3 = "lat"
x!4 = "lon"
x&series = ispan(1983, 1983+nn-1, 1) ;ispan(1983, 2013, 1)
x&var = nvname
x&lat = lat
x&lon = lon
x&lev = lev
x&lat@units = "degrees_north"
x&lon@units = "degrees_east"
x&lev@units = "hPa"
; x@long_name = "forcing"
; x@units     = "W/m2"
x@long_name = "partial T"
x@units     = "K"

do nnn = 0, nn-1
    fili2 = "/home/yangsong3/wuyt/sysu/scripts/TP_surrounding/cfram20210212/partial_T_"+(nnn+1)+".grd"
    print(fili2)
    do nv = 0,nvars-1
       do nl = 0,nlev-1
          x(nnn,nv,nl,:,:) = fbindirread(fili2,nv*nlev+nl,(/nlat,nlon/),"float")
       end do
    end do
end do

if (any(isnan_ieee(x))) then   ; if x got a element is Not A Number, then the Fill_value
  value = -999.
  replace_ieeenan (x,value,0)
  x@_FillValue = value
end if
printVarSummary(x)
; printMinMax(x(0,:,:,:),1)

x_atmos = x(:,:,0:nlev-2,:,:)
printVarSummary(x_atmos)
;===========================================================
; vertical integrated and area averaged
;===========================================================

; base ps for integrated

ps_jja = month_to_season(ps, season)
ps_base = dim_avg_n_Wrap(ps_jja(0:9,:,:), 0)
printVarSummary(ps_base)

dP_base = dpres_plevel_Wrap(plev, ps_base, min(plev), 0)
printVarSummary(dP_base)
dp_base = conform_dims(dimsizes(x_atmos), dP_base, (/2,3,4/))
copy_VarCoords(x_atmos, dp_base)

xdp = x_atmos*dp_base
copy_VarCoords(x_atmos, xdp)
printVarSummary(xdp)

x_tropo = dim_sum_n_Wrap(xdp(:,:,{LEVs:LEVe},:,:), 2)/dim_sum_n_Wrap(dp_base(:,:,{LEVs:LEVe},:,:), 2)
copy_VarCoords(x(:,:,0,:,:), x_tropo)

x_series = wgt_areaave_Wrap(x_tropo(:,:,{latstr_east:latend_east},{lonstr_east:lonend_east}), 1., 1., 0)
printVarSummary(x_series)

nvname2 = (/"albedo","wv","cloud","co2","atm_dyn","ocean_dyn","shflx","lhflx"/)  ; ,"error"

x_final = x_series(:,0:7)
x_final(:,0:2) = x_series(:,0:2)
x_final(:,3) = x_series(:,5)
x_final(:,4:7) = x_series(:,9:12)
x_final&var = nvname2
printVarSummary(x_final)

x_total = dim_sum_n_Wrap(x_final, 1)
printVarSummary(x_total)


;========================================================
zero = new(nyear, float)
zero = 0.

y = new((/2,nyear/), float)
y(0,:) = t_runmean
y(1,4:4+nn-1) = (/x_total/)

d1 = x_final(:,0)
d2 = x_final(:,1)
d3 = x_final(:,2)
d4 = x_final(:,3)
d5 = x_final(:,4)
d6 = x_final(:,5)
d7 = x_final(:,6)
d8 = x_final(:,7)
; d9 = x_final(:,8)
; d0 = x_final(:,9)

;==================================================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_surrounding/cfram20210212/ncl/"
figname = "east"

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
res@trYMinF                  = -2.2
res@trYMaxF                  = 2.2;(max(d2) + max(d1)+max(d3) + max(d4))*1.2  
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
res@gsnLeftString = "(b) East region: ("+lonstr_east+"-"+lonend_east+"E,"+latstr_east+"-"+latend_east+"N)"
res@gsnRightString = " "
res@tiXAxisString = "Year"
res@tiYAxisFontHeightF = 0.015
res@tiYAxisString = "Partial T (K)"

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
  bres@gsFillIndex          = 17
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
  bres@gsFillIndex          = 5
  bres@gsEdgeThicknessF     = 2.5

  dum7(i-1) = gsn_add_polygon(wks,plot_without_xgrid,xbar,ybar,bres)

;*************************************************
  if (d8(i-1).ge.0) then
    ybar_p_max = ybar_p_max + d8(i-1)
    ybar_positive = (/ybar_p_min,ybar_p_max,ybar_p_max,ybar_p_min,ybar_p_min/)
    ybar = ybar_positive
    ybar_p_min = ybar_p_max
  end if

  if (d8(i-1).lt.0) then
    ybar_n_max = ybar_n_max + d8(i-1)
    ybar_negative = (/ybar_n_min,ybar_n_max,ybar_n_max,ybar_n_min,ybar_n_min/)
    ybar = ybar_negative
    ybar_n_min = ybar_n_max
  end if

  bres@gsFillColor          = "magenta1"
  bres@gsEdgeColor          = "magenta1"
  bres@gsFillIndex          = 8
  bres@gsEdgeThicknessF     = 2.5

  dum8(i-1) = gsn_add_polygon(wks,plot_without_xgrid,xbar,ybar,bres)


;*************************************************
  delete([/ybar_p_min, ybar_p_max, ybar_n_min, ybar_n_max,\
           ybar_positive, ybar_negative, ybar/])

end do



; Create some labelbars for a legend
lb_colors = (/"coral2","goldenrod2","darkviolet","lightskyblue3","green3","dodgerblue","wheat3","magenta1"/)
lb_dash = (/0,0,16,1,17,0,5,8/)

lb_x = new((/10,5/), "float")
; lb_x(0,:) = (/0.5,0.5,2.,2.,0.5/)
lb_x(0,:) = (/1984.5,1984.5,1986.,1986.,1984.5/)
do i = 1, 9
  lb_x(i,:) = lb_x(i-1,:) + 3.0
end do

lb_y = (/-1.9, -1.7, -1.7, -1.9, -1.9/)

dum_lb = new(10, graphic)
lb_res = True
lb_res@gsEdgesOn = True
lb_res@gsEdgeThicknessF = 2.5
do i = 0, 7
  lb_res@gsFillColor = lb_colors(i)
  lb_res@gsEdgeColor = lb_colors(i)
  lb_res@gsFillIndex = lb_dash(i)
  dum_lb(i) = gsn_add_polygon(wks, plot_without_xgrid, lb_x(i,:), lb_y, lb_res)
end do

text = (/"AL","WV","CLD","CO2","ATM","HS","SH","LH","Sum","ERA5"/) ;(/"AL","WV","CLD","ATM","OCN","HF","Sum","Glb Mean WPS"/)
text_res = True
text_res@txFontHeightF = 0.012
text_res@txJust             = "CenterLeft"
text_plot = new(10, graphic)
text_x = new(10, "float")
text_x(0) = 1986.1
do i = 1, 7
  text_x(i) = text_x(i-1) + 3.0
end do
do i = 8, 9
  text_x(i) = text_x(i-1) + 2.8
end do
text_y = -1.8
do i = 0, 9
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
line_y = (/-1.8,-1.8/)
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