
begin

nn = 3

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

f1 = addfile(diri+fili0, "r")
Time = f1->time
YYYY = cd_calendar(Time, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
strd = short2flt(f1->strd(iYYYY,:,:))
strd = strd/86400.

lat = strd&latitude
lon = strd&longitude
nlat = dimsizes(lat)
nlon = dimsizes(lon)
;===================================================================================================
nvname = (/"rad_conv","lwu","lwd","swu","swd","lwu_toa","lwd_toa","swu_toa","swd_toa"/)  ; ,"error"
nvars = dimsizes(nvname)

diri1 = "/home/yangsong3/wuyt/sysu/scripts/TP_surrounding/cfram20210212/"
diri2 = "/home/yangsong3/wuyt/sysu/scripts/TP_surrounding/cfram20210513_clearsky/"
filib1 = diri1+"baseline_radsfc_ranc_1.grd"
filib2 = diri2+"baseline_radsfc_ranc_1.grd"
filit1 = diri1+"ts_radsfc_ranc_"+nn+".grd"
filit2 = diri2+"ts_radsfc_ranc_"+nn+".grd"

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

y = x
ts1 = x
ts2 = x

do nv = 0,nvars-1
    x(nv,:,:) = fbindirread(filib1,nv,(/nlat,nlon/),"float")
    y(nv,:,:) = fbindirread(filib2,nv,(/nlat,nlon/),"float")
    ts1(nv,:,:) = fbindirread(filit1,nv,(/nlat,nlon/),"float")
    ts2(nv,:,:) = fbindirread(filit2,nv,(/nlat,nlon/),"float")
end do
printVarSummary(x)

nx = 1
dR_full = ts1(nx,:,:) - x(nx,:,:)
dR_clear = ts2(nx,:,:) - x(nx,:,:)
copy_VarCoords(x(nx,:,:), dR_full)
copy_VarCoords(x(nx,:,:), dR_clear)

diff = dR_full - dR_clear
copy_VarCoords(dR_full, diff)
;=======================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_surrounding/20210423/"
figname = "check_surf_upR_clearsky_fullsky_"+nn

wks = gsn_open_wks("pdf", figpath+figname)
gsn_define_colormap(wks, "temp_diff_18lev")
plot = new(3, graphic)
plot_topo = plot


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


res@gsnLeftString = "(a) Full sky upR"
plot(0) = gsn_csm_contour_map(wks,dR_full,res)
plot_topo(0) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(0),plot_topo(0))

res@gsnLeftString = "(b) Clear sky upR"
plot(1) = gsn_csm_contour_map(wks,dR_clear,res)
plot_topo(1) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(1),plot_topo(1))

res@gsnLeftString = "(c) a-b"
plot(2) = gsn_csm_contour_map(wks,diff,res)
plot_topo(2) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(2),plot_topo(2))

;-------------------------------------------------------

pres = True
pres@gsnMaximize = True
pres@gsnFrame    = True
pres@gsnDraw     = True
pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = "surface upwelling LW"

pres@lbLabelFontHeightF = 0.015
pres@pmLabelBarWidthF = 0.8    ;plot manager 
pres@pmLabelBarHeightF = 0.06
gsn_panel(wks,plot,(/1,3/),pres)

end