
begin

lonstr = 60
lonend = 150
latstr = 15
latend = 50


yrStrt = 198102
yrLast = 200601


casename = (/"Ctrl_FAMIP","Hist_FAMIP"/) ;,"Hist_IO","Hist_PO","Hist_AO","Hist_PIO","Hist_PAO","Hist_IAO"
ncase = dimsizes(casename)

season = "DJF"
nlat = 96
nlon = 144

LEV = 300

VT_bar = new((/ncase,nlat,nlon/), float)
VT_bar!0 = "case"
VT_bar&case = casename
VT_ano = VT_bar
VT_tot = VT_bar

UT_bar = VT_bar
UT_ano = VT_bar
UT_tot = VT_bar

do i = 0, ncase-1
    diri = "/home/yangsong3/data-model/wuyt/TPbias/TPbias_"+casename(i)+"/a20191206/"

    fili1 = casename(i)+".cam.h0.T.1979-2005.nc"
    fili2 = casename(i)+".cam.h0.V.1979-2005.nc"    
    fili3 = casename(i)+".cam.h0.U.1979-2005.nc"
    f1 = addfile(diri+fili1, "r")
    f2 = addfile(diri+fili2, "r")
    f3 = addfile(diri+fili3, "r")

    TIME = f1->time
    YYYY = cd_calendar(TIME, -1)
    iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

    T = f1->T(iYYYY,{LEV},:,:)
    V = f2->V(iYYYY,{LEV},:,:)
    U = f3->U(iYYYY,{LEV},:,:)

    T_season = month_to_season(T, season)
    V_season = month_to_season(V, season)
    U_season = month_to_season(U, season)

    T_bar = dim_avg_n_Wrap(T_season, 0)
    V_bar = dim_avg_n_Wrap(V_season, 0)
    U_bar = dim_avg_n_Wrap(U_season, 0)
    VT_bar(i,:,:) = V_bar*T_bar
    UT_bar(i,:,:) = U_bar*T_bar
    copy_VarCoords(V_bar, VT_bar(i,:,:))
    copy_VarCoords(V_bar, UT_bar(i,:,:))

    T_ano = T_season - conform_dims(dimsizes(T_season), T_bar, (/1,2/))
    V_ano = V_season - conform_dims(dimsizes(V_season), V_bar, (/1,2/))
    U_ano = U_season - conform_dims(dimsizes(U_season), V_bar, (/1,2/))
    VT_ano_med = T_ano*V_ano
    VT_ano(i,:,:) = (/dim_avg_n_Wrap(VT_ano_med, 0)/)
    copy_VarCoords(V_bar, VT_ano(i,:,:))
    UT_ano_med = T_ano*U_ano
    UT_ano(i,:,:) = (/dim_avg_n_Wrap(UT_ano_med, 0)/)
    copy_VarCoords(U_bar, UT_ano(i,:,:))

    VT_tot(i,:,:) = VT_bar(i,:,:)+VT_ano(i,:,:)
    copy_VarCoords(V_bar, VT_tot(i,:,:))
    UT_tot(i,:,:) = UT_bar(i,:,:)+UT_ano(i,:,:)
    copy_VarCoords(U_bar, UT_tot(i,:,:))

end do

printVarSummary(VT_bar)
printVarSummary(VT_ano)
printVarSummary(VT_tot)

; diff_bar = VT_bar(0,:,:)
; ; diff_ano = VT_ano(0,:,:)
; diff_bar = VT_bar(1,:,:) - VT_bar(0,:,:)
; diff_ano = VT_ano(1,:,:) - VT_ano(0,:,:)

diff_vt_tot = VT_tot(0,:,:)
diff_vt_tot = VT_tot(1,:,:) - VT_tot(0,:,:)
diff_ut_tot = UT_tot(0,:,:)
diff_ut_tot = UT_tot(1,:,:) - UT_tot(0,:,:)

;==========================================================================
; -d（VT）/dy; -d（UT）/dx
;==========================================================================
lat = tofloat(VT_bar&lat)
lon = tofloat(VT_bar&lon)
nlat = dimsizes(lat)
nlon = dimsizes(lon)

dut_dx = VT_bar(0,:,:)
dvt_dy = dut_dx

dlon = (lon(2)-lon(1))*0.0174533
do i=0,nlat-1                      ; loop over each latitude
    dx = 6378388.*cos(0.0174533*lat(i))*dlon  ; constant at this latitude
    dut_dx(i,:) = center_finite_diff_n (diff_ut_tot(i,:), dx , True,0,0) ; be careful with the dimison downscale
end do
dut_dx = dut_dx*(-1.)

dlat = (lat(2)-lat(1))*0.0174533
dy   = 6378388.*dlat 
dvt_dy = center_finite_diff_n (diff_vt_tot, dy , False,0,0)
dvt_dy = dvt_dy*(-1.)

printVarSummary(dut_dx)
printVarSummary(dvt_dy)
;===========================================================================

nvname = (/"rad_conv","lwu","lwd","swu","swd","lwu_toa","lwd_toa","swu_toa","swd_toa"/)  ; ,"error"
nvars = dimsizes(nvname)

diri2 = "/home/yangsong3/wuyt/sysu/scripts/TP_coldbias/model/a20200910/cfram_decomposeT/"
fili2 = diri2+"baseline_radsfc_ranc_1.grd"
fili5 = diri2+"t_dyn_radsfc_ranc_1.grd"

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

dyn = x

do nv = 0,nvars-1
    x(nv,:,:) = fbindirread(fili2,nv,(/nlat,nlon/),"float")
    dyn(nv,:,:) = fbindirread(fili5,nv,(/nlat,nlon/),"float")
end do

nx = 2

diff_dyn = x(0,:,:)
diff_dyn = dyn(nx,:,:)-x(nx,:,:)
printVarSummary(diff_dyn)

;=========================================================
; outline of tp
;=========================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)

;=========================================================
; plot
;=========================================================

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/a20200910/"
figname = season+"_Ctrl_diff_heat_flux_VT_UT"


wks = gsn_open_wks("pdf", figpath+figname)
plot = new(8, graphic)
plot_topo = plot

 res_topo            = True
 res_topo@cnFillOn   = False
 res_topo@cnLinesOn  = True
 res_topo@gsnFrame   = False
 res_topo@gsnDraw    = False
 res_topo@cnLineLabelsOn         = False
 res_topo@cnLineThicknessF       = 3
 res_topo@cnLineColor            = "purple"
 res_topo@cnInfoLabelOn          = False
 res_topo@cnLevelSelectionMode   ="ExplicitLevels"
 res_topo@cnLevels               = (/1500/)
 res_topo@gsnLeftString = ""
 res_topo@gsnRightString = ""

res                       = True     ; plot mods desired
res@gsnDraw               = False
res@gsnFrame              = False

res@gsnMaximize           = True     ; Make plot fit the frame
res@cnFillOn              = True     ; turn on color fill
res@cnLinesOn             = False    ; turn of contour lines
res@cnLineLabelsOn        = False
res@cnInfoLabelOn         = False
res@cnFillPalette         = "temp_diff_18lev"
;res@lbLabelBarOn          = False
res@pmLabelBarOrthogonalPosF = 0.2

res@mpMinLatF            = latstr      ; range to zoom in on
res@mpMaxLatF            = latend
res@mpMinLonF            = lonstr
res@mpMaxLonF            = lonend
; res@mpCenterLonF         = 180.
res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025

res@gsnStringFontHeightF = 0.035


; res@gsnRightString = "Ctrl"
res@gsnLeftString = "(a) Ctrl <VT>"
res@cnFillPalette = "temp_diff_18lev"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -2000
res@cnMaxLevelValF = 2000
res@cnLevelSpacingF = 200.
plot(0) = gsn_csm_contour_map(wks, VT_tot(0,:,:), res)
plot_topo(0) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(0), plot_topo(0))

res@gsnLeftString = "(b) Ctrl <UT>"
res@cnFillPalette = "WhiteBlueGreenYellowRed"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = 0
res@cnMaxLevelValF = 20000
res@cnLevelSpacingF = 1000.
plot(1) = gsn_csm_contour_map(wks, UT_tot(0,:,:), res)
plot_topo(1) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(1), plot_topo(1))


; res@gsnRightString = "Diff"
res@gsnLeftString = "(c) Diff. <VT>"
res@cnFillPalette = "temp_diff_18lev"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -200
res@cnMaxLevelValF = 200
res@cnLevelSpacingF = 20.
plot(2) = gsn_csm_contour_map(wks, diff_vt_tot, res)
plot_topo(2) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(2), plot_topo(2))

res@gsnLeftString = "(d) Diff. <UT>"
res@cnFillPalette = "temp_diff_18lev"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -1000
res@cnMaxLevelValF = 1000
res@cnLevelSpacingF = 100.
plot(3) = gsn_csm_contour_map(wks, diff_ut_tot, res)
plot_topo(3) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(3), plot_topo(3))

res@gsnLeftString = "(e) diff. -d<VT>/dy"
res@cnFillPalette = "temp_diff_18lev"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -1
res@cnMaxLevelValF = 1
res@cnLevelSpacingF = 0.1
dvt_dy = dvt_dy*1.e+3
res@gsnRightString = "*10~S~-3~N~"
plot(4) = gsn_csm_contour_map(wks, dvt_dy, res)
plot_topo(4) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(4), plot_topo(4))

res@gsnLeftString = "(f) diff. -d<UT>/dx"
res@cnFillPalette = "temp_diff_18lev"
; res@cnLevelSelectionMode = "ExplicitLevels"
; res@cnMinLevelValF = -1000
; res@cnMaxLevelValF = 1000
; res@cnLevelSpacingF = 100
dut_dx = dut_dx*1.e+3
plot(5) = gsn_csm_contour_map(wks, dut_dx, res)
plot_topo(5) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(5), plot_topo(5))

delete(res@gsnRightString)
res@gsnLeftString = "(g) fc_atm_dyn"
res@cnFillPalette = "temp_diff_18lev"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -10
res@cnMaxLevelValF = 10
res@cnLevelSpacingF = 1.
plot(6) = gsn_csm_contour_map(wks, diff_dyn, res)
plot_topo(6) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(6), plot_topo(6))

pres = True
pres@gsnMaximize = True
;pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = season+" mean heat flux: "+LEV+" hPa "

gsn_panel(wks, plot, (/4,2/), pres)

end