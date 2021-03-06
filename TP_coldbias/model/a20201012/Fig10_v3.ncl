
begin

lonstr = 60
lonend = 110
latstr = 20
latend = 50


yrStrt = 198102
yrLast = 200601


casename = (/"Ctrl_FAMIP","Hist_FAMIP","Hist_FAMIP_CMIP6","Hist_PO"/) ;"Hist_PO","Hist_IO","Hist_PO","Hist_AO","Hist_PIO","Hist_PAO","Hist_IAO"
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
; exit

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


; ;===========================================================================

; yrStrt = 1979
; yrLast = 2005

; TS = new((/4,324,nlat,nlon/), float)

; dirisst = "/home/yangsong3/wuyt/sysu/scripts/TP_coldbias/data/"
; fsst1 = "sst_HadOIBl_bc_1.9x2.5_1850_2012_c130411.nc"
; fsst2 = "Hist_bias_sst_cmip5_bc_1.9x2.5_c20191204.nc"
; fsst3 = "Hist_bias_sst_cmip6_bc_1.9x2.5_c20201004.nc"
; fsst4 = "PO_sst_cmip5_bc_1.9x2.5_c20191204.nc"
; fs1 = addfile(dirisst+fsst1, "r")
; fs2 = addfile(dirisst+fsst2, "r")
; fs3 = addfile(dirisst+fsst3, "r")
; fs4 = addfile(dirisst+fsst4, "r")

; TIMEs = fs1->time
; YYYYs = cd_calendar(TIMEs, -1)/100
; iYYYYs = ind(YYYYs.ge.yrStrt .and. YYYYs.le.yrLast)

; TS(0,:,:,:) = fs1->SST_cpl(iYYYYs,:,:)
; TS(1,:,:,:) = fs2->SST_cpl(iYYYYs,:,:)
; TS(2,:,:,:) = fs3->SST_cpl(iYYYYs,:,:)
; TS(3,:,:,:) = fs4->SST_cpl(iYYYYs,:,:)

; sst_diff = TS(0,:,:,:)
; sst_DJF = TS(0:2,27,:,:)
; do i = 0, 2
;     sst_diff = TS(i+1,:,:,:) - TS(0,:,:,:)
;     sst_DJFx = month_to_season(sst_diff, "DJF")
;     sst_DJF(i,:,:) = dim_avg_n_Wrap(sst_DJFx, 0)
; end do

; printVarSummary(sst_DJF)

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
figname = "Fig10_v3"


wks = gsn_open_wks("pdf", figpath+figname)
plot = new(4, graphic)
plotc = plot
plot_topo = plot

 res_topo            = True
 res_topo@cnFillOn   = False
 res_topo@cnLinesOn  = True
 res_topo@gsnFrame   = False
 res_topo@gsnDraw    = False
 res_topo@cnLineLabelsOn         = False
 res_topo@cnLineThicknessF       = 3
 res_topo@cnLineColor            = "gray40"
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
; res@lbLabelBarOn          = False
res@pmLabelBarOrthogonalPosF = 0.12

res@mpMinLatF            = latstr      ; range to zoom in on
res@mpMaxLatF            = latend
res@mpMinLonF            = lonstr
res@mpMaxLonF            = lonend
; res@mpCenterLonF         = 180.
res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025
res@gsnStringFontHeightF = 0.025

;====================================================================================
cnres = True
cnres@gsnDraw  = False
cnres@gsnFrame = False

cnres@cnFillOn = False
cnres@cnLinesOn = True
cnres@cnLineLabelsOn = False
cnres@cnLineLabelFontColor = "red"
cnres@cnLineColor = "dodgerblue"
cnres@cnLineThicknessF = 2
cnres@cnInfoLabelOn = True
cnres@cnInfoLabelOrthogonalPosF = -0.15
cnres@cnInfoLabelFontHeightF = 0.02
cnres@cnLabelDrawOrder = "PostDraw"
cnres@gsnContourNegLineDashPattern = 16
cnres@gsnContourPosLineDashPattern = 0

; cnres@cnLevelSelectionMode = "ManualLevels"
; cnres@cnMinLevelValF = -5.
; cnres@cnMaxLevelValF = 5.
; cnres@cnLevelSpacingF = 0.5
cnres@cnLevelSelectionMode = "ExplicitLevels"
cnres@cnLevels = (/-5,-4,-3,-2,-1,-0.5,0.5,1,2,3,4,5/)
; cnres@cnLineColors = (/5,6,7,8,9,10,11,12,13,14,15,16,17/)
cnres@gsnLeftString = ""
cnres@gsnRightString = ""
;====================================================================================

res@gsnRightString = "*10~S~-3~N~"
; res@gsnLeftString = "(a) Diff: -{d[<V>*<T>]/dy + d[<U>*<T>]/dy}"
res@gsnLeftString = "(b) SEN-ALL - CTRL"
res@cnFillPalette = "temp_diff_18lev"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -0.5
res@cnMaxLevelValF = 0.5
res@cnLevelSpacingF = 0.05
plot(1) = gsn_csm_contour_map(wks, diff_dtot(0,:,:), res)
; plotc(1) = gsn_csm_contour(wks, sst_DJF(0,:,:), cnres)
; overlay(plot(1), plotc(1))
plot_topo(1) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(1), plot_topo(1))

res@gsnLeftString = "(d) SEN-ALL6 - CTRL"
plot(3) = gsn_csm_contour_map(wks, diff_dtot(1,:,:), res)
; plotc(3) = gsn_csm_contour(wks, sst_DJF(1,:,:), cnres)
; overlay(plot(3), plotc(3))
plot_topo(3) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(3), plot_topo(3))

res@gsnLeftString = "(c) SEN-PAC - CTRL"
plot(2) = gsn_csm_contour_map(wks, diff_dtot(2,:,:), res)
; plotc(2) = gsn_csm_contour(wks, sst_DJF(2,:,:), cnres)
; overlay(plot(2), plotc(2))
plot_topo(2) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(2), plot_topo(2))

res@gsnLeftString = "(a) Atm. dyn."
res@gsnRightString = "W/m~S~2~N~"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -10
res@cnMaxLevelValF = 10
res@cnLevelSpacingF = 1
plot(0) = gsn_csm_contour_map(wks, diff_dyn, res)
plot_topo(0) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(0), plot_topo(0))


pres = True
pres@gsnMaximize = True
; pres@gsnPanelLabelBar = True
; pres@pmLabelBarOrthogonalPosF = -0.02
; pres@lbLabelFontHeightF = 0.012
; pres@gsnPanelMainString = season+" mean "

gsn_panel(wks, plot, (/2,2/), pres)

end