
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

diff_vt_bar = VT_bar(0,:,:)
diff_vt_ano = VT_ano(0,:,:)
diff_vt_bar = VT_bar(1,:,:) - VT_bar(0,:,:)
diff_vt_ano = VT_ano(1,:,:) - VT_ano(0,:,:)

diff_ut_bar = UT_bar(0,:,:)
diff_ut_ano = UT_ano(0,:,:)
diff_ut_bar = UT_bar(1,:,:) - UT_bar(0,:,:)
diff_ut_ano = UT_ano(1,:,:) - UT_ano(0,:,:)

; diff_vt_tot = VT_tot(0,:,:)
; diff_vt_tot = VT_tot(1,:,:) - VT_tot(0,:,:)
; diff_ut_tot = UT_tot(0,:,:)
; diff_ut_tot = UT_tot(1,:,:) - UT_tot(0,:,:)

;==========================================================================
; -d（VT）/dy; -d（UT）/dx
;==========================================================================
lat = tofloat(VT_bar&lat)
lon = tofloat(VT_bar&lon)
nlat = dimsizes(lat)
nlon = dimsizes(lon)

dutbar_dx = VT_bar(0,:,:)
dutano_dx = dutbar_dx
dvtbar_dy = dutbar_dx
dvtano_dy = dutbar_dx

dlon = (lon(2)-lon(1))*0.0174533
do i=0,nlat-1                      ; loop over each latitude
    dx = 6378388.*cos(0.0174533*lat(i))*dlon  ; constant at this latitude
    dutbar_dx(i,:) = center_finite_diff_n (diff_ut_bar(i,:), dx , True,0,0) ; be careful with the dimison downscale
    dutano_dx(i,:) = center_finite_diff_n (diff_ut_ano(i,:), dx , True,0,0)
end do
dutbar_dx = dutbar_dx*(-1000.) 
dutano_dx = dutano_dx*(-1000.)

dlat = (lat(2)-lat(1))*0.0174533
dy   = 6378388.*dlat 
dvtbar_dy = center_finite_diff_n (diff_vt_bar, dy , False,0,0)
dvtbar_dy = dvtbar_dy*(-1000.)
dvtano_dy = center_finite_diff_n (diff_vt_ano, dy , False,0,0)
dvtano_dy = dvtano_dy*(-1000.)
; print(dvtano_dy({30:40},{80:90}))
; print(dutano_dx({30:40},{80:90}))

; exit

printVarSummary(dutano_dx)
printVarSummary(dvtano_dy)

duttot_dx = dutano_dx+dutbar_dx
dvttot_dy = dvtano_dy+dvtbar_dy
copy_VarCoords(VT_bar(0,:,:), duttot_dx)
copy_VarCoords(VT_bar(0,:,:), dvttot_dy)

dbar_tot = dutbar_dx+dvtbar_dy
dano_tot = dutano_dx+dvtano_dy
dtot_tot = duttot_dx+dvttot_dy
copy_VarCoords(VT_bar(0,:,:), dbar_tot)
copy_VarCoords(VT_bar(0,:,:), dano_tot)
copy_VarCoords(VT_bar(0,:,:), dtot_tot)

; ;===========================================================================

; nvname = (/"rad_conv","lwu","lwd","swu","swd","lwu_toa","lwd_toa","swu_toa","swd_toa"/)  ; ,"error"
; nvars = dimsizes(nvname)

; diri2 = "/home/yangsong3/wuyt/sysu/scripts/TP_coldbias/model/a20200910/cfram_decomposeT/"
; fili2 = diri2+"baseline_radsfc_ranc_1.grd"
; fili5 = diri2+"t_dyn_radsfc_ranc_1.grd"

; UNDEF = -999.
; x = new((/nvars,nlat,nlon/),float,UNDEF)
; x!0 = "var"
; x!1 = "lat"
; x!2 = "lon"
; x&var = nvname
; x&lat = lat
; x&lon = lon
; x&lat@units = "degrees_north"
; x&lon@units = "degrees_east"
; x@long_name = "forcing"
; x@units     = "W/m2"

; dyn = x

; do nv = 0,nvars-1
;     x(nv,:,:) = fbindirread(fili2,nv,(/nlat,nlon/),"float")
;     dyn(nv,:,:) = fbindirread(fili5,nv,(/nlat,nlon/),"float")
; end do

; nx = 2

; diff_dyn = x(0,:,:)
; diff_dyn = dyn(nx,:,:)-x(nx,:,:)
; printVarSummary(diff_dyn)

;=========================================================
; outline of tp
;=========================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)

;=========================================================
; plot
;=========================================================

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/a20200910/"
figname = season+"_Partial_xy_heat_flux_VT_UT"


wks = gsn_open_wks("pdf", figpath+figname)
plot = new(6, graphic)
plotc = plot
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
res@lbLabelBarOn          = False
; res@pmLabelBarOrthogonalPosF = 0.2

res@mpMinLatF            = latstr      ; range to zoom in on
res@mpMaxLatF            = latend
res@mpMinLonF            = lonstr
res@mpMaxLonF            = lonend
; res@mpCenterLonF         = 180.
res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025

res@gsnStringFontHeightF = 0.035


res@cnFillPalette = "temp_diff_18lev"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -0.5
res@cnMaxLevelValF = 0.5
res@cnLevelSpacingF = 0.05
res@gsnRightString = "*10~S~-3~N~"

res@gsnLeftString = "(a) diff. -d[<V>*<T>]/dy"
plot(0) = gsn_csm_contour_map(wks, dvtbar_dy, res)
plot_topo(0) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(0), plot_topo(0))

res@gsnLeftString = "(b) diff. -d[<V'*T'>]/dy"
plot(1) = gsn_csm_contour_map(wks, dvtano_dy, res)
plot_topo(1) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(1), plot_topo(1))

res@gsnLeftString = "(c) a+b"
plot(2) = gsn_csm_contour_map(wks, dvttot_dy, res)
plot_topo(2) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(2), plot_topo(2))

res@gsnLeftString = "(d) diff. -d[<U>*<T>]/dx"
plot(3) = gsn_csm_contour_map(wks, dutbar_dx, res)
plot_topo(3) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(3), plot_topo(3))

res@gsnLeftString = "(e) diff. -d[<U'*T'>]/dx"
plot(4) = gsn_csm_contour_map(wks, dutano_dx, res)
plot_topo(4) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(4), plot_topo(4))

res@gsnLeftString = "(f) d+e"
plot(5) = gsn_csm_contour_map(wks, duttot_dx, res)
plot_topo(5) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(5), plot_topo(5))

res@gsnLeftStringFontHeightF = 0.025
res@gsnLeftString = "(a) diff. -{d[<V>*<T>]/dy + d[<U>*<T>]/dx}"
; res@cnFillPalette = "temp_diff_18lev"
; res@cnLevelSelectionMode = "ManualLevels"
; res@cnMinLevelValF = -10
; res@cnMaxLevelValF = 10
; res@cnLevelSpacingF = 1.
plot(6) = gsn_csm_contour_map(wks, dbar_tot, res)
plot_topo(6) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(6), plot_topo(6))


res@gsnLeftString = "(b) -{d[<V'*T'>]/dy + d[<U'*T'>]/dx}"
plot(7) = gsn_csm_contour_map(wks, dano_tot, res)
plot_topo(7) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(7), plot_topo(7))

res@gsnLeftString = "(c) a+b"
;"(c) -{d[<V>*<T>]/dy + d[<U>*<T>]/dx + d[<V'*T'>]/dy + d[<U>*<T>]/dx}"
plot(8) = gsn_csm_contour_map(wks, dtot_tot, res)
plot_topo(8) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(8), plot_topo(8))

pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
pres@pmLabelBarOrthogonalPosF = -0.02
pres@lbLabelFontHeightF = 0.015
pres@gsnPanelMainString = season+" mean: "+LEV+" hPa "

gsn_panel(wks, plot(0:5), (/2,3/), pres)
gsn_panel(wks, plot(6:8), (/1,3/), pres)

end