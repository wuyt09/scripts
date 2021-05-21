
begin

lonstr = -20
lonend = 160
latstr = 0
latend = 70

yrStrt = 1979
yrLast = 2019
year = ispan(yrStrt, yrLast, 1)
season = "JJA"

diri = "/home/yangsong3/data-observation/ERA5-monthly/surface/"
fili0 = "surface_thermal_downwards_mon_1x1_1979_2019.nc"
fili1 = "surface_thermal_downwards_clear_sky_mon_1x1_1979_2019.nc"
fili2 = "sur_pressure_mon_1x1_1979_2019.nc"

f1 = addfile(diri+fili0, "r")
f2 = addfile(diri+fili1, "r")
f3 = addfile(diri+fili2, "r")

Time = f1->time
YYYY = cd_calendar(Time, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

strd = short2flt(f1->strd(iYYYY,:,:))
strd = strd/86400.
strdc = short2flt(f2->strdc(iYYYY,:,:))
strdc = strdc/86400.

strd_JJA = month_to_season(strd, season)
strd_JJA&time = year
strdc_JJA = month_to_season(strdc, season)
copy_VarCoords(strd_JJA, strdc_JJA)
cre = strd_JJA - strdc_JJA
copy_VarCoords(strd_JJA, cre)

yrs1 = 1979
yre1 = 1988

yrs = ispan(1980, 2010, 1)
yre = ispan(1989, 2019, 1)
nn = 23
yrs2 = yrs(nn-1)
yre2 = yre(nn-1)
print("========================================")
era_sr = dim_avg_n_Wrap(strd_JJA({yrs2:yre2},:,:), 0) - dim_avg_n_Wrap(strd_JJA({yrs1:yre1},:,:), 0)
cre_sr = dim_avg_n_Wrap(cre({yrs2:yre2},:,:), 0) - dim_avg_n_Wrap(cre({yrs1:yre1},:,:), 0)
copy_VarCoords(strd(0,:,:), era_sr)
copy_VarCoords(strd(0,:,:), cre_sr)

lat = strd&latitude
lon = strd&longitude
nlat = dimsizes(lat)
nlon = dimsizes(lon)

;===================================================================================================
nvname = (/"rad_conv","lwu","lwd","swu","swd","lwu_toa","lwd_toa","swu_toa","swd_toa"/)  ; ,"error"
nvars = dimsizes(nvname)

diri2 = "/home/yangsong3/wuyt/sysu/scripts/TP_surrounding/cfram20210513_clearsky/"
; diri2 = "/home/yangsong3/wuyt/sysu/scripts/TP_surrounding/cfram20210212/"
fili2 = diri2+"baseline_radsfc_ranc_1.grd"
; fili3 = diri2+"wv_radsfc_ranc_31.grd"
; fili4 = diri2+"ta_radsfc_ranc_31.grd"
; fili5 = diri2+"albedo_radsfc_ranc_31.grd"
; fili6 = diri2+"co2_radsfc_ranc_31.grd"
; fili7 = diri2+"solar_radsfc_ranc_31.grd"
; fili8 = diri2+"o3_radsfc_ranc_31.grd"
; fili9 = diri2+"warm_radsfc_ranc_31.grd"
; fili10 = diri2+"cloud_radranc_10.grd"
fili3 = diri2+"wv_radsfc_ranc_"+nn+".grd"
fili4 = diri2+"ta_radsfc_ranc_"+nn+".grd"
fili5 = diri2+"albedo_radsfc_ranc_"+nn+".grd"
fili6 = diri2+"co2_radsfc_ranc_"+nn+".grd"
fili7 = diri2+"solar_radsfc_ranc_"+nn+".grd"
fili8 = diri2+"o3_radsfc_ranc_"+nn+".grd"
fili9 = diri2+"ts_radsfc_ranc_"+nn+".grd"
fili10 = diri2+"warm_radsfc_ranc_"+nn+".grd"

UNDEF = -999.
x = new((/nvars,nlat,nlon/),float,UNDEF)
x!0 = "var"
x!1 = "lat"
x!2 = "lon"
x&var = nvname
x&lat = lat
x&lon = lon
x&lat@units = "degrees_north"
x&lon@units = "degrees_east"
x@long_name = "forcing"
x@units     = "W/m2"

wv = x
t = x
al = x
co2 = x
solar = x
o3 = x
ts = x
y = x
cld = x

do nv = 0,nvars-1
    x(nv,:,:) = fbindirread(fili2,nv,(/nlat,nlon/),"float")
    wv(nv,:,:) = fbindirread(fili3,nv,(/nlat,nlon/),"float")
    t(nv,:,:) = fbindirread(fili4,nv,(/nlat,nlon/),"float")
    al(nv,:,:) = fbindirread(fili5,nv,(/nlat,nlon/),"float")
    co2(nv,:,:) = fbindirread(fili6,nv,(/nlat,nlon/),"float")
    solar(nv,:,:) = fbindirread(fili7,nv,(/nlat,nlon/),"float")
    o3(nv,:,:) = fbindirread(fili8,nv,(/nlat,nlon/),"float")
    ts(nv,:,:) = fbindirread(fili9,nv,(/nlat,nlon/),"float")
    y(nv,:,:) = fbindirread(fili10,nv,(/nlat,nlon/),"float")
    ; cld(nv,:,:) = fbindirread(fili10, nv*nlev+ilev, (/nlat,nlon/), "float")
end do
printVarSummary(x)
; printMinMax(x(0,:,:,:),1)


process = (/"AL","WV","Ta","Ts","CO2","Solar","O3","CRE","Warm - Base","CFRAM total","Obs total"/)
nx = 2
x_all = new((/dimsizes(process),nlat,nlon/), float)
x_all!0 = "var"
x_all!1 = "lat"
x_all!2 = "lon"
x_all&var = process
x_all&lat = lat
x_all&lon = lon
x_all&lat@units = "degrees_north"
x_all&lon@units = "degrees_east"
x_all@long_name = "forcing"
x_all@units     = "W/m2"

x_all(0,:,:) = al(nx,:,:) - x(nx,:,:)
x_all(1,:,:) = wv(nx,:,:) - x(nx,:,:)
x_all(2,:,:) = t(nx,:,:) - x(nx,:,:)
x_all(3,:,:) = ts(nx,:,:) - x(nx,:,:)
x_all(4,:,:) = co2(nx,:,:) - x(nx,:,:)
x_all(5,:,:) = solar(nx,:,:) - x(nx,:,:)
x_all(6,:,:) = o3(nx,:,:) - x(nx,:,:)
x_all(7,:,:) = cre_sr
x_all(8,:,:) = y(nx,:,:) - x(nx,:,:)
x_all(9,:,:) = dim_sum_n(x_all(0:6,:,:), 0)
x_all(10,:,:) = era_sr
copy_VarCoords(x(0,:,:), x_all(0,:,:))
printVarSummary(x_all)

;+ x_all(0,:,:) + x_all(3,:,:) + x_all(5,:,:)
; diff_cfram = x_all(1,:,:) + x_all(2,:,:) + x_all(4,:,:) + x_all(6,:,:) - x_all(8,:,:)
; copy_VarCoords(x_all(0,:,:), diff_cfram)


diff_cfram = x_all(9,:,:) - x_all(8,:,:)
diff_obs = x_all(10,:,:) - x_all(9,:,:) - x_all(7,:,:)
copy_VarCoords(x_all(0,:,:), diff_cfram)
copy_VarCoords(x_all(0,:,:), diff_obs)

diff = wgt_areaave_Wrap(diff_cfram({25:40},{75:105}), 1, 1, 0)
; diff = wgt_areaave_Wrap(x_all(5,{25:40},{75:105}), 1, 1, 0)
; print(diff)
; exit
;=======================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_surrounding/20210423/"
figname = "check_surf_downR_clearsky_"+nn

Left = (/"(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)","(j)","(k)","(l)"/)

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

res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025
res@tmXTOn = False
res@tmYROn = False

;res@gsnLeftStringFontHeightF  = 0.02
;res@gsnRightStringFontHeightF  = 0.02
res@gsnStringFontHeightF = 0.035

; res@cnLevelSelectionMode = "ManualLevels"
; res@cnMaxLevelValF = 3
; res@cnMinLevelValF = -3
; res@cnLevelSpacingF = 0.5
res@lbLabelBarOn         = False
res@cnLevelSelectionMode = "ExplicitLevels"
res@cnLevels = (/-3,-2.5,-2,-1.5,-1,-0.5,-0.1,0.1,0.5,1,1.5,2,2.5,3/)
res@cnFillColors = (/4,5,6,7,8,9,10,11,12,13,14,15,16,17,18/)


do i = 0, 10
    res@gsnLeftString = Left(i)+" "+process(i)
    plot(i) = gsn_csm_contour_map(wks,x_all(i,:,:),res)
    plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(i),plot_topo(i))
end do

res@gsnLeftString = Left(11)+" k - (j+h)"
plot(11) = gsn_csm_contour_map(wks,y(nx,:,:),res)
plot_topo(11) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(11),plot_topo(11))

;-------------------------------------------------------

pres = True
pres@gsnMaximize = True
pres@gsnFrame    = True
pres@gsnDraw     = True
pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = "surface downwelling LW:"+yrs2+"-"+yre2+" minus "+yrs1+"-"+yre1

pres@lbLabelFontHeightF = 0.015
pres@pmLabelBarWidthF = 0.8    ;plot manager 
pres@pmLabelBarHeightF = 0.06
gsn_panel(wks,plot,(/3,4/),pres)


end