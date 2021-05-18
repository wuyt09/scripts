
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
fili0 = "surface_solar_downwards_mon_1x1_1979_2019.nc"
fili1 = "surface_solar_downwards_clear_sky_mon_1x1_1979_2019.nc"

f1 = addfile(diri+fili0, "r")
f2 = addfile(diri+fili1, "r")

Time = f1->time
YYYY = cd_calendar(Time, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

ssrd = short2flt(f1->ssrd(iYYYY,:,:))
ssrd = ssrd/86400.
ssrdc = short2flt(f2->ssrdc(iYYYY,:,:))
ssrdc = ssrdc/86400.

strd_JJA = month_to_season(ssrd, season)
strd_JJA&time = year
strdc_JJA = month_to_season(ssrdc, season)
copy_VarCoords(strd_JJA, strdc_JJA)
cre = strd_JJA - strdc_JJA
copy_VarCoords(strd_JJA, cre)

yrs1 = 1979
yre1 = 1988
yrs2 = 1982 ;1983
yre2 = 1991 ;1992
print("========================================")
era_sr = dim_avg_n_Wrap(strd_JJA({yrs2:yre2},:,:), 0) - dim_avg_n_Wrap(strd_JJA({yrs1:yre1},:,:), 0)
cre_sr = dim_avg_n_Wrap(cre({yrs2:yre2},:,:), 0) - dim_avg_n_Wrap(cre({yrs1:yre1},:,:), 0)
copy_VarCoords(ssrd(0,:,:), era_sr)
copy_VarCoords(ssrd(0,:,:), cre_sr)


lat = ssrd&latitude
lon = ssrd&longitude
nlat = dimsizes(lat)
nlon = dimsizes(lon)

;===================================================================================================
nvname = (/"rad_conv","lwu","lwd","swu","swd","lwu_toa","lwd_toa","swu_toa","swd_toa"/)  ; ,"error"
nvars = dimsizes(nvname)

diri2 = "/home/yangsong3/wuyt/sysu/scripts/TP_surrounding/cfram20210423_clearsky/"
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
fili3 = diri2+"wv_radsfc_ranc_3.grd"
fili4 = diri2+"ta_radsfc_ranc_3.grd"
fili5 = diri2+"albedo_radsfc_ranc_3.grd"
fili6 = diri2+"co2_radsfc_ranc_3.grd"
fili7 = diri2+"solar_radsfc_ranc_3.grd"
fili8 = diri2+"o3_radsfc_ranc_3.grd"
fili9 = diri2+"warm_radsfc_ranc_3.grd"

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
    y(nv,:,:) = fbindirread(fili9,nv,(/nlat,nlon/),"float")
    ; cld(nv,:,:) = fbindirread(fili10, nv*nlev+ilev, (/nlat,nlon/), "float")
end do
printVarSummary(x)
; printMinMax(x(0,:,:,:),1)


process = (/"wv","t","al","co2","solar","o3","cre","warm","total"/)
nx = 4
x_all = new((/dimsizes(process),nlat,nlon/), float)
x_all(0,:,:) = wv(nx,:,:) - x(nx,:,:)
x_all(1,:,:) = t(nx,:,:) - x(nx,:,:)
x_all(2,:,:) = al(nx,:,:) - x(nx,:,:)
x_all(3,:,:) = co2(nx,:,:) - x(nx,:,:)
x_all(4,:,:) = solar(nx,:,:) - x(nx,:,:)
x_all(5,:,:) = o3(nx,:,:) - x(nx,:,:)
x_all(6,:,:) = cre_sr
x_all(7,:,:) = y(nx,:,:) - x(nx,:,:)
x_all(8,:,:) = dim_sum_n(x_all(0:6,:,:), 0)
copy_VarCoords(x(0,:,:), x_all(0,:,:))
printVarSummary(x_all)

;=======================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_surrounding/20210423/"
figname = "check_surf_downS_clearsky"


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

res@cnLevelSelectionMode = "ManualLevels"
res@cnMaxLevelValF = 3
res@cnMinLevelValF = -3
res@cnLevelSpacingF = 0.5
res@lbLabelBarOn         = False

res@gsnLeftString = "(a) ERA5"
plot(0) = gsn_csm_contour_map(wks,era_sr,res)
plot_topo(0) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(0),plot_topo(0))

; res@gsnLeftString = "(b) CFRAM S-R (warm - base)"
; plot(1) = gsn_csm_contour_map(wks,x_all(7,:,:),res)
; plot_topo(1) = gsn_csm_contour(wks, topo, res_topo)
; overlay(plot(1),plot_topo(1))

res@gsnLeftString = "(b) CFRAM total + CRE"
plot(1) = gsn_csm_contour_map(wks,x_all(8,:,:),res)
plot_topo(1) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(1),plot_topo(1))

diff = x_all(8,:,:) - era_sr
copy_VarCoords(x_all(8,:,:),diff)
res@gsnLeftString = "(c) b - a"
plot(2) = gsn_csm_contour_map(wks,diff,res)
plot_topo(2) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(2),plot_topo(2))

; res@cnMaxLevelValF = 5
; res@cnMinLevelValF = -5
; res@cnLevelSpacingF = 1
; res@gsnLeftString = "(d) Cloud Radiative Effect"
; plot(3) = gsn_csm_contour_map(wks,cre_sr,res)
; plot_topo(3) = gsn_csm_contour(wks, topo, res_topo)
; overlay(plot(3),plot_topo(3))


res@gsnLeftString = "(e) WV"
plot(4) = gsn_csm_contour_map(wks,x_all(0,:,:),res)
plot_topo(4) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(4),plot_topo(4))

res@gsnLeftString = "(f) T"
plot(5) = gsn_csm_contour_map(wks,x_all(1,:,:),res)
plot_topo(5) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(5),plot_topo(5))

res@gsnLeftString = "(g) AL"
plot(6) = gsn_csm_contour_map(wks,x_all(2,:,:),res)
plot_topo(6) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(6),plot_topo(6))

res@gsnLeftString = "(h) CO2"
plot(7) = gsn_csm_contour_map(wks,x_all(3,:,:),res)
plot_topo(7) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(7),plot_topo(7))

; res@gsnLeftString = "(i) Solar"
; plot(8) = gsn_csm_contour_map(wks,x_all(4,:,:),res)
; plot_topo(8) = gsn_csm_contour(wks, topo, res_topo)
; overlay(plot(8),plot_topo(8))
res@gsnLeftString = "(i) Cloud Radiative Effect"
plot(8) = gsn_csm_contour_map(wks,cre_sr,res)
plot_topo(8) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(8),plot_topo(8))

res@gsnLeftString = "(d) o3"
plot(3) = gsn_csm_contour_map(wks,x_all(5,:,:),res)
plot_topo(3) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(3),plot_topo(3))

res@gsnLeftString = "(j) Solar"
plot(9) = gsn_csm_contour_map(wks,x_all(4,:,:),res)
plot_topo(9) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(9),plot_topo(9))

; res@gsnLeftString = "(k) Cloud"
; plot(10) = gsn_csm_contour_map(wks,cld(nx,:,:),res)
; plot_topo(10) = gsn_csm_contour(wks, topo, res_topo)
; overlay(plot(10),plot_topo(10))

;-------------------------------------------------------

pres = True
pres@gsnMaximize = True
pres@gsnFrame    = True
pres@gsnDraw     = True
pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = "surface downwelling SW:"+yrs2+"-"+yre2+" minus "+yrs1+"-"+yre1

pres@lbLabelFontHeightF = 0.015
pres@pmLabelBarWidthF = 0.8    ;plot manager 
pres@pmLabelBarHeightF = 0.06
gsn_panel(wks,plot,(/4,3/),pres)


end