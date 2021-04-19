
begin

lonstr = -20
lonend = 160
latstr = 0
latend = 70

yrStrt = 1979
yrLast = 2019
year = ispan(yrStrt, yrLast, 1)
nyear = dimsizes(year)
season = "JJA"

; LEVs = 1000
; LEVe = 300

latstr_east = 35
latend_east = 45
lonstr_east = 105
lonend_east = 140

latstr_west = 30
latend_west = 55
lonstr_west = -10
lonend_west = 60

latstr_south = 10
latend_south = 25
lonstr_south = 50
lonend_south = 100

latstr_north = 45
latend_north = 60
lonstr_north = 85
lonend_north = 120

latstr_tp = 25 ;28
latend_tp = 40 ;38
lonstr_tp = 75 ;80
lonend_tp = 105

;=============================================================================================================
diri = "/home/yangsong3/data-observation/ERA5-monthly/surface/"
fili = "t2m_mon_1x1_1979_2019.nc"

f = addfile(diri+fili, "r")

Time = f->time
YYYY = cd_calendar(Time, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

t = short2flt(f->t2m(iYYYY,:,:))
printVarSummary(t)


t_tropo_JJA = month_to_season(t, season)


fa      = addfile ("/home/yangsong3/wuyt/data/lsm_1279l4_0.1x0.1.grb_v4_unpack.nc","r")
; lsdata = doubletofloat(a->lsm(0,:,:))
lsdata = fa->lsm(0,:,:)
printVarSummary(lsdata)
lsm = linint2_Wrap(lsdata&longitude, lsdata&latitude, lsdata, True, t&longitude, t&latitude, 0)
printVarSummary(lsm)
; delete(lsdata)
; fa = addfile("/home/yangsong3/wuyt/data/landsea.nc", "r")
; lsdata = tofloat(fa->LSMASK(::-1,:))
; lsdata!0 = "lat"
; lsdata!1 = "lon"
; lsdata&lat = fspan(89.5, -89.5, 180)
; lsdata&lon = fspan(0.5, 359.5, 360)
; lsm = linint2_Wrap(lsdata&lon, lsdata&lat, lsdata, True, t&longitude, t&latitude, 0)
; printVarSummary(lsm)


lsm1 = conform_dims(dimsizes(t_tropo_JJA), lsm, (/1,2/))
t_tropo_JJA = where(lsm1.ge.0.8, t_tropo_JJA, t@_FillValue) ;

t_base = dim_avg_n_Wrap(t_tropo_JJA(0:9,:,:), 0)
t_warm = dim_avg_n_Wrap(t_tropo_JJA(31:40,:,:), 0)
t_new = t_warm - t_base
copy_VarCoords(t_base, t_new)
printVarSummary(t_new)

;=============================================================================================================

nvname = (/"albedo","wv","cloud","cloud_sw","cloud_lw","co2","o3","solar",\
           "dyn","atm_dyn","ocean_dyn","shflx","lhflx"/)  ; ,"error"
nn = 31

f = addfile("/home/yangsong3/data-observation/ERA5-monthly/pressure/t_mon_1x1_1979_2019.nc", "r")
lat  = f->latitude ;(::-1)
lon  = f->longitude
level  = f->level
nlev = 38
lev = new(nlev, integer, "No_FillValue")
lev(0:nlev-2) = level
lev(nlev-1) = 1001
; print(lev)

delete(f)

nlat = dimsizes(lat)
nlon = dimsizes(lon)
nlev = dimsizes(lev)
nvars = dimsizes(nvname)

setfileoption("bin","ReadByteOrder","Native")

UNDEF = -999.
x = new((/1,nvars,nlev,nlat,nlon/),float,UNDEF)
x!0 = "series"
x!1 = "var"
x!2 = "lev"
x!3 = "lat"
x!4 = "lon"
x&series = 2013 ;ispan(1983, 1983+nn-1, 1) ;ispan(1983, 2013, 1)
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

; do nnn = 0, nn-1
    fili2 = "/home/yangsong3/wuyt/sysu/scripts/TP_surrounding/cfram20210212/partial_T_"+nn+".grd"
    print(fili2)
    do nv = 0,nvars-1
       do nl = 0,nlev-1
          x(0,nv,nl,:,:) = fbindirread(fili2,nv*nlev+nl,(/nlat,nlon/),"float")
       end do
    end do
; end do

if (any(isnan_ieee(x))) then   ; if x got a element is Not A Number, then the Fill_value
  value = -999.
  replace_ieeenan (x,value,0)
  x@_FillValue = value
end if
printVarSummary(x)
; printMinMax(x(0,:,:,:),1)

x_atmos = x(:,:,nlev-1,:,:)
printVarSummary(x_atmos)
;===========================================================
; vertical integrated and area averaged
;===========================================================

x_tropo = x_atmos

lsm2 = conform_dims(dimsizes(x_tropo), lsm, (/2,3/))
x_tropo = where(lsm2.ge.0.8, x_tropo, UNDEF) ;

; x_series = wgt_areaave_Wrap(x_tropo(:,:,{latstr_south:latend_south},{lonstr_south:lonend_south}), 1., 1., 0)
; printVarSummary(x_series)

nvname2 = (/"albedo","wv","cloud","co2","atm_dyn","ocean_dyn","shflx","lhflx"/)  ; ,"error"

x_final = x_tropo(:,0:7,:,:)
x_final(:,0:2,:,:) = x_tropo(:,0:2,:,:)
x_final(:,3,:,:) = x_tropo(:,5,:,:)
x_final(:,4:7,:,:) = x_tropo(:,9:12,:,:)
x_final&var = nvname2
printVarSummary(x_final)

x_total = dim_sum_n_Wrap(x_final, 1)
printVarSummary(x_total)


x_inn = x_tropo(0,:,:,:)
;=======================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_surrounding/20210415/"
figname = "spatial_pattern_process_partial_T_31"

wks = gsn_open_wks("pdf", figpath+figname)
gsn_define_colormap(wks, "temp_diff_18lev")
plot = new(12, graphic)
plot_topo = plot

dum1 = plot
dum2 = plot
dum3 = plot
dum4 = plot
;=======================================
; this is for outline of TP
;=======================================

ff = "~/wuyt/data/TOPO_1500.nc"
topo_in = addfile(ff,"r")
topo = topo_in ->TP

mlat = topo_in->LATIXY(:,0)
mlon = topo_in->LONGXY(0,:)
topo!0 = "lsmlat"
topo!1 = "lsmlon"
topo&lsmlat = mlat
topo&lsmlon = mlon

;----------------------
res_topo            = True
res_topo@cnFillOn   = False
res_topo@cnLinesOn  = True
res_topo@gsnFrame   = False
res_topo@gsnDraw    = False
res_topo@cnLineLabelsOn         = False
res_topo@cnLineThicknessF       = 2
res_topo@cnLineColor            = "grey32"
;res_topo@cnLineDashPattern      = "1"
res_topo@cnInfoLabelOn          = False
res_topo@cnLevelSelectionMode   ="ExplicitLevels"
res_topo@cnLevels               = (/1500/)



;==========================================================

res = True
res@gsnDraw   = False
res@gsnFrame  = False
;res@gsnSpreadColors = True
res@gsnAddCyclic = True

res@mpFillOn  = True
;res@mpMinLonF = 180.
;res@mpMaxLonF = 360.
res@mpMinLatF = latstr
res@mpMaxLatF = latend
res@mpMinLonF = lonstr
res@mpMaxLonF = lonend
;res@mpCenterLonF = 85.
res@mpGeophysicalLineThicknessF = 2

res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False

;res@lbLabelAutoStride = True
;res@lbLabelFontHeightF = 0.01
;res@pmLabelBarOrthogonalPosF = 0.16  ; 水平位置
;res@pmLabelBarHeightF = 0.02    ;色标高度
;res@pmLabelBarWidthF = 0.7
;res@pmLabelBarHeightF = 0.15
;res@lbLabelFontHeightF = 0.016
res@lbLabelBarOn         = False 
res@lbLabelsOn           = False

res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025
res@tmXTOn = False
res@tmYROn = False

;res@gsnLeftStringFontHeightF  = 0.02
;res@gsnRightStringFontHeightF  = 0.02
res@gsnStringFontHeightF = 0.035
res@gsnRightString = "K" ;W/m~S~2~N~

;-------------------------------------------------------

res@cnLevelSelectionMode = "ExplicitLevels"
res@cnLevels = (/-3,-2.5,-2,-1.5,-1,-0.5,-0.1,0.1,0.5,1,1.5,2,2.5,3/)
res@cnFillColors = (/4,5,6,7,8,9,10,11,12,13,14,15,16,17,18/)

; res@cnLevelSelectionMode = "ManualLevels"
; res@cnMaxLevelValF = 5
; res@cnMinLevelValF = -5
; res@cnLevelSpacingF = 0.5


res@gsnLeftString = "(a) Heat Storage"
plot(0) = gsn_csm_contour_map(wks,x_inn(10,:,:),res)
plot_topo(0) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(0),plot_topo(0))

res@gsnLeftString = "(b) Surf.LH"
plot(1) = gsn_csm_contour_map(wks,x_inn(12,:,:),res)
plot_topo(1) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(1),plot_topo(1))

res@gsnLeftString = "(c) Surf.SH"
plot(2) = gsn_csm_contour_map(wks,x_inn(11,:,:),res)
plot_topo(2) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(2),plot_topo(2))


;**************************************************
;delete(res@cnMaxLevelValF)
;delete(res@cnMinLevelValF)
;res@cnMinLevelValF = -1.5
;res@cnMaxLevelValF = 1.5
;res@cnLevelSpacingF = 0.1

res@gsnLeftString = "(d) Albedo"
plot(3) = gsn_csm_contour_map(wks,x_inn(0,:,:),res)
plot_topo(3) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(3),plot_topo(3))


res@gsnLeftString = "(e) Cloud"
plot(4) = gsn_csm_contour_map(wks,x_inn(2,:,:),res)
plot_topo(4) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(4),plot_topo(4))

res@gsnLeftString = "(f) Atm_Dyn"
plot(5) = gsn_csm_contour_map(wks,x_inn(9,:,:),res)
plot_topo(5) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(5),plot_topo(5))

;***************************************************

;res@cnMinLevelValF = -1
;res@cnMaxLevelValF = 1
;res@cnLevelSpacingF = 0.05

res@gsnLeftString = "(g) WV"
plot(6) = gsn_csm_contour_map(wks,x_inn(1,:,:),res)
plot_topo(6) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(6),plot_topo(6))

res@gsnLeftString = "(h) CO2"
plot(7) = gsn_csm_contour_map(wks,x_inn(5,:,:),res)
plot_topo(7) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(7),plot_topo(7))

res@gsnLeftString = "(i) sum of (a)-(h)"
plot(8) = gsn_csm_contour_map(wks,x_total(0,:,:),res)
plot_topo(8) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(8),plot_topo(8))



res@gsnLeftString = "(j) Obs total"
plot(9) = gsn_csm_contour_map(wks,t_new,res)
plot_topo(9) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(9),plot_topo(9))

error = t_new - x_total(0,:,:)
copy_VarCoords(t_new, error)
res@gsnLeftString = "(k) error: (j) - (i)"
plot(10) = gsn_csm_contour_map(wks,error,res)
plot_topo(10) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(10),plot_topo(10))


;------------------------------------------------------
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

do i = 0,10
    dum1(i) = gsn_add_polyline(wks, plot(i), glon1, glat1, gres)
    dum2(i) = gsn_add_polyline(wks, plot(i), glon2, glat2, gres)
    dum3(i) = gsn_add_polyline(wks, plot(i), glon3, glat3, gres)
    dum4(i) = gsn_add_polyline(wks, plot(i), glon4, glat4, gres)
end do

;-------------------------------------------------------

pres = True
pres@gsnMaximize = True
pres@gsnFrame    = True
pres@gsnDraw     = True
pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = "JJA T2m"

pres@lbLabelFontHeightF = 0.015
pres@pmLabelBarWidthF = 0.8    ;plot manager 
pres@pmLabelBarHeightF = 0.06
gsn_panel(wks,plot,(/3,4/),pres)


end