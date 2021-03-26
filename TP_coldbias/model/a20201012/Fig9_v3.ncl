
begin

lonstr = 60
lonend = 150
latstr = 15
latend = 50


yrStrt = 198102
yrLast = 200601


casename = (/"Ctrl_FAMIP","Hist_FAMIP"/) ;,"Hist_PO","Hist_FAMIP_CMIP6","Hist_IO","Hist_AO","Hist_PIO","Hist_PAO","Hist_IAO"
ncase = dimsizes(casename)

season = "DJF"
nlat = 96
nlon = 144

LEV = 300

VT_season = new((/ncase,nlat,nlon/), float)
VT_season!0 = "case"
VT_season&case = casename
UT_season = VT_season

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

    ; T_season = month_to_season(T, season)
    ; V_season = month_to_season(V, season)
    ; U_season = month_to_season(U, season)

    ; T_bar = dim_avg_n_Wrap(T_season, 0)
    ; V_bar = dim_avg_n_Wrap(V_season, 0)
    ; U_bar = dim_avg_n_Wrap(U_season, 0)
    ; VT_season(i,:,:) = V_bar*T_bar
    ; UT_season(i,:,:) = U_bar*T_bar

    UT = U*T
    VT = V*T
    copy_VarCoords(U, UT)
    copy_VarCoords(V, VT)

    UT_season_med = month_to_season(UT, season)
    VT_season_med = month_to_season(VT, season)

    UT_season(i,:,:) = dim_avg_n_Wrap(UT_season_med, 0)
    VT_season(i,:,:) = dim_avg_n_Wrap(VT_season_med, 0)

end do

copy_VarCoords(U(0,:,:), UT_season(0,:,:))
copy_VarCoords(V(0,:,:), VT_season(0,:,:))
printVarSummary(VT_season)

diff_vt = VT_season(0:ncase-2,:,:)
diff_ut = UT_season(0:ncase-2,:,:)
do i = 0, ncase-2
    diff_vt(i,:,:) = VT_season(i+1,:,:) - VT_season(0,:,:)
    diff_ut(i,:,:) = UT_season(i+1,:,:) - UT_season(0,:,:)
end do
tot_season = UT_season + VT_season
copy_VarCoords(UT_season, tot_season)

diff_tot = diff_ut+diff_vt
copy_VarCoords(diff_ut, diff_tot)
;==========================================================================
; -d（VT）/dy; -d（UT）/dx
;==========================================================================
lat = tofloat(V&lat)
lon = tofloat(V&lon)
nlat = dimsizes(lat)
nlon = dimsizes(lon)

dut_dx = VT_season(0:ncase-1,:,:)
dvt_dy = dut_dx

dlon = (lon(2)-lon(1))*0.0174533
do ll=0,ncase-1
    do i=0,nlat-1                      ; loop over each latitude
        dx = 6378388.*cos(0.0174533*lat(i))*dlon  ; constant at this latitude
        dut_dx(ll,i,:) = center_finite_diff_n (UT_season(ll,i,:), dx , True,0,0) ; be careful with the dimison downscale
end do
end do
dut_dx = dut_dx*(-1000.)

dlat = (lat(2)-lat(1))*0.0174533
dy   = 6378388.*dlat
do ll=0,ncase-1
    dvt_dy(ll,:,:) = center_finite_diff_n (VT_season(ll,:,:), dy , False,0,0)
end do
dvt_dy = dvt_dy*(-1000.)
printVarSummary(dut_dx)
printVarSummary(dvt_dy)

dtot = dut_dx+dvt_dy
copy_VarCoords(VT_season(0:ncase-1,:,:), dtot)

diff_dvtdy = dvt_dy(0:ncase-2,:,:)
diff_dutdx = dut_dx(0:ncase-2,:,:)
do i = 0, ncase-2
    diff_dutdx(i,:,:) = dut_dx(i+1,:,:) - dut_dx(0,:,:)
    diff_dvtdy(i,:,:) = dvt_dy(i+1,:,:) - dvt_dy(0,:,:)
end do
copy_VarCoords(dut_dx(0:ncase-2,:,:), diff_dutdx)
copy_VarCoords(dut_dx(0:ncase-2,:,:), diff_dvtdy)
diff_dtot = diff_dutdx+diff_dvtdy
copy_VarCoords(diff_dutdx, diff_dtot)
;=========================================================
; outline of tp
;=========================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)
topo = smth9_Wrap(topo, 0.5, 0.25, True)
;=========================================================
; plot
;=========================================================

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/a20201012/"
figname = "Fig9_v3"


wks = gsn_open_wks("pdf", figpath+figname)
plot = new(6, graphic)
plotc = plot
plot_topo = plot
plotc = plot

 res_topo            = True
 res_topo@cnFillOn   = False
 res_topo@cnLinesOn  = True
 res_topo@gsnFrame   = False
 res_topo@gsnDraw    = False
 res_topo@cnLineLabelsOn         = False
 res_topo@cnLineThicknessF       = 3
 res_topo@cnLineColor            = "grey40"
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

;====================================================================================
cnres = True
cnres@gsnDraw  = False
cnres@gsnFrame = False

cnres@cnFillOn = False
cnres@cnLinesOn = True
cnres@cnLineLabelsOn = False
cnres@cnLineLabelFontColor = "red"
cnres@cnLineColor = "purple"
cnres@cnLineThicknessF = 4
cnres@cnInfoLabelOn = True
cnres@cnInfoLabelOrthogonalPosF = -0.15
cnres@cnInfoLabelFontHeightF = 0.02
cnres@cnLabelDrawOrder = "PostDraw"
cnres@gsnContourNegLineDashPattern = 16
cnres@gsnContourPosLineDashPattern = 0
;====================================================================================


res@cnFillPalette = "temp_diff_18lev"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -3.
res@cnMaxLevelValF = 3.
res@cnLevelSpacingF = 0.5
res@gsnRightString = "*10~S~-3~N~"

cnres@cnLevelSelectionMode = "ManualLevels"
cnres@cnMinLevelValF = -2000
cnres@cnMaxLevelValF = 2000
cnres@cnLevelSpacingF = 1000

res@gsnLeftString = "(a) CTRL. -d(V*T)/dy" ;C/~F34~6~Cy
plot(0) = gsn_csm_contour_map(wks, dvt_dy(0,:,:), res)
plot_topo(0) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(0), plot_topo(0))
plotc(0) = gsn_csm_contour(wks, VT_season(0,:,:), cnres)
overlay(plot(0), plotc(0))

cnres@cnLevelSelectionMode = "ManualLevels"
cnres@cnMinLevelValF = -20000
cnres@cnMaxLevelValF = 20000
cnres@cnLevelSpacingF = 5000
res@gsnLeftString = "(b) CTRL. -d(U*T)/dx"
plot(1) = gsn_csm_contour_map(wks, dut_dx(0,:,:), res)
plot_topo(1) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(1), plot_topo(1))
plotc(1) = gsn_csm_contour(wks, UT_season(0,:,:), cnres)
overlay(plot(1), plotc(1))

res@gsnLeftString = "(c) a+b"
plot(2) = gsn_csm_contour_map(wks, dtot(0,:,:), res)
plot_topo(2) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(2), plot_topo(2))
; plotc(2) = gsn_csm_contour(wks, tot_season(0,:,:), cnres)
; overlay(plot(2), plotc(2))



res@gsnRightString = "*10~S~-3~N~"
cnres@cnLevelSelectionMode = "ManualLevels"
cnres@cnMinLevelValF = -200
cnres@cnMaxLevelValF = 200
cnres@cnLevelSpacingF = 100
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -0.5
res@cnMaxLevelValF = 0.5
res@cnLevelSpacingF = 0.05
res@gsnLeftString = "(d) diff. -d(V*T)/dy"
plot(3) = gsn_csm_contour_map(wks, diff_dvtdy(0,:,:), res)
plot_topo(3) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(3), plot_topo(3))
plotc(3) = gsn_csm_contour(wks, diff_vt(0,:,:), cnres)
overlay(plot(3), plotc(3))

cnres@cnMinLevelValF = -1000
cnres@cnMaxLevelValF = 1000
cnres@cnLevelSpacingF = 500
res@gsnLeftString = "(e) diff. -d(U*T)/dx"
plot(4) = gsn_csm_contour_map(wks, diff_dutdx(0,:,:), res)
plot_topo(4) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(4), plot_topo(4))
plotc(4) = gsn_csm_contour(wks, diff_ut(0,:,:), cnres)
overlay(plot(4), plotc(4))

res@gsnLeftString = "(f) d+e"
plot(5) = gsn_csm_contour_map(wks, diff_dtot(0,:,:), res)
plot_topo(5) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(5), plot_topo(5))
; plotc(5) = gsn_csm_contour(wks, diff_tot(0,:,:), cnres)
; overlay(plot(5), plotc(5))


; pres = True
; pres@gsnMaximize = True
; pres@gsnPanelLabelBar = True
; ; pres@pmLabelBarOrthogonalPosF = -0.02
; pres@lbLabelFontHeightF = 0.015
; ; pres@gsnPanelMainString = season+" mean: "+LEV+" hPa "

; gsn_panel(wks, plot, (/2,3/), pres)

pres = True
; pres@gsnMaximize = True
pres@gsnDraw  = False
pres@gsnFrame = False

height = 0.2
pres@gsnDraw = True
pres@gsnPanelLabelBar = True
pres@lbLabelFontHeightF = 0.012
pres@gsnPanelTop = 0.7
pres@gsnPanelBottom = pres@gsnPanelTop - height
gsn_panel(wks, plot(0:2), (/1,3/), pres)


pres@gsnPanelLabelBar = True
pres@gsnPanelTop = pres@gsnPanelBottom
pres@gsnPanelBottom = pres@gsnPanelTop - height
pres@gsnPanelLabelBar = True
pres@lbLabelFontHeightF = 0.012
gsn_panel(wks, plot(3:5), (/1,3/), pres)

end