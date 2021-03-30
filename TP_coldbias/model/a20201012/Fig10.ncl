
begin

lonstr = 80
lonend = 105
latstr = 28
latend = 36


yrStrt = 198102
yrLast = 200601


casename = (/"Ctrl_FAMIP","Hist_FAMIP"/) ;,"Hist_PO","Hist_FAMIP_CMIP6","Hist_IO","Hist_AO","Hist_PIO","Hist_PAO","Hist_IAO"
ncase = dimsizes(casename)

season = "DJF"
nlat = 96
nlon = 144
nlev = 17
; LEV = 200

VT_season = new((/ncase,nlev,nlat,nlon/), float)
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

    T = f1->T(iYYYY,:,:,:)
    V = f2->V(iYYYY,:,:,:)
    U = f3->U(iYYYY,:,:,:)

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

    UT_season(i,:,:,:) = dim_avg_n_Wrap(UT_season_med, 0)
    VT_season(i,:,:,:) = dim_avg_n_Wrap(VT_season_med, 0)

end do

copy_VarCoords(UT(0,:,:,:), UT_season(0,:,:,:))
copy_VarCoords(VT(0,:,:,:), VT_season(0,:,:,:))
printVarSummary(VT_season)

diff_vt = VT_season(0:ncase-2,:,:,:)
diff_ut = UT_season(0:ncase-2,:,:,:)
do i = 0, ncase-2
    diff_vt(i,:,:,:) = VT_season(i+1,:,:,:) - VT_season(0,:,:,:)
    diff_ut(i,:,:,:) = UT_season(i+1,:,:,:) - UT_season(0,:,:,:)
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

dut_dx = VT_season(0:ncase-1,:,:,:)
dvt_dy = dut_dx

dlon = (lon(2)-lon(1))*0.0174533
do ll=0,ncase-1
    do i=0,nlat-1                      ; loop over each latitude
        dx = 6378388.*cos(0.0174533*lat(i))*dlon  ; constant at this latitude
        dut_dx(ll,:,i,:) = center_finite_diff_n (UT_season(ll,:,i,:), dx , True,0,1) ; be careful with the dimison downscale
end do
end do
dut_dx = dut_dx*(-1000.)


dlat = (lat(2)-lat(1))*0.0174533
dy   = 6378388.*dlat
do ll=0,ncase-1
    dvt_dy(ll,:,:,:) = center_finite_diff_n (VT_season(ll,:,:,:), dy , False,0,1)
end do
dvt_dy = dvt_dy*(-1000.)
printVarSummary(dut_dx)
printVarSummary(dvt_dy)

dutdx_profile = dim_avg_n_Wrap(dut_dx(:,:,{latstr:latend},:), 2)
dvtdy_profile = dim_avg_n_Wrap(dvt_dy(:,:,{latstr:latend},:), 2)
printVarSummary(dutdx_profile)

dtot = dutdx_profile+dvtdy_profile
copy_VarCoords(VT_season(0:ncase-1,:,{latstr},:), dtot)

diff_dvtdy = dvtdy_profile(0:ncase-2,:,:)
diff_dutdx = dutdx_profile(0:ncase-2,:,:)
do i = 0, ncase-2
    diff_dutdx(i,:,:) = dutdx_profile(i+1,:,:) - dutdx_profile(0,:,:)
    diff_dvtdy(i,:,:) = dvtdy_profile(i+1,:,:) - dvtdy_profile(0,:,:)
end do
copy_VarCoords(dutdx_profile(0:ncase-2,:,:), diff_dutdx)
copy_VarCoords(dvtdy_profile(0:ncase-2,:,:), diff_dvtdy)
printVarSummary(diff_dutdx)

diff_dtot = diff_dutdx+diff_dvtdy
copy_VarCoords(diff_dutdx, diff_dtot)
printVarSummary(diff_dtot)
;=========================================================
; outline of tp
;=========================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)
topo_z = dim_avg_n_Wrap(topo({latstr:latend},:), 0)
topo_zonal = 1013.25*(1-topo_z*0.0065/288.15)^5.25145
copy_VarCoords(topo({latstr},:), topo_zonal)
;=========================================================
;=========================================================
; plot
;=========================================================

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/a20201012/"
figname = "Fig10"


wks = gsn_open_wks("pdf", figpath+figname)
plot = new(1, graphic)
plotc = plot
plot_topo = plot
plotc = plot

;=====================================================
    res_topo                   = True

    res_topo@gsnDraw  = False           ;先别直接画出来
    res_topo@gsnFrame = False           ;也不要翻页

    res_topo@gsnLeftString  = " "              ;不要标示出变量名和单位等信息
    res_topo@gsnRightString = " "

    ;res_topo@gsnMaximize        = True
    ;res_topo@trYMaxF            = 1000.0
    res_topo@gsnYRefLine        = 1000 ;(/100,1000/)          ; create a reference line 
    res_topo@gsnBelowYRefLineColor = "black" ;(/"black","black"/) ; above ref line fill grey
    ;res_topo@gsnXYBelowFillColors = "goldenrod2"
    res_topo@gsnAboveYRefLineColor = "blue"
    res_topo@xyLineColors = "black"
    res_topo@xyLineThicknesses = 4.
;===================================================================
;===================================================================

    res = True
    res@gsnDraw  = False
    res@gsnFrame = False
    res@vpHeightF = 0.4
    res@vpWidthF  = 0.7

    res@cnFillOn = True
    res@cnLinesOn = False
    res@cnLineLabelsOn = False
    res@cnInfoLabelOn = False

    res@gsnStringFontHeightF = 0.03
    res@tmXBLabelFontHeightF = 0.02
    res@tmYLLabelFontHeightF = 0.02

    res@pmLabelBarOrthogonalPosF = -0.06
    ; res@lbLabelBarOn = False
    ;res@lbLabelFontHeightF = 0.02
    res@tiYAxisString = ""
    res@cnFillDrawOrder = "PreDraw"

    res@tmYLMode      = "Explicit"
    res@tmYLValues    = (/"1000","850","750","500","300","250"/) ;,"150"
    res@tmYLLabels    = (/"1000","850","750","500","300","250"/) ;,"150"

    res@trYMinF   = 250
    res@trYMaxF   = 1000
    res@gsnPresHgtHeightLabelOn = False
    res@tmYROn = False

;=====================================================

;----------------------------------------------
; draw prob
;----------------------------------------------
 res3  = True
 res3@gsnDraw   = False
 res3@gsnFrame  = False

 ;res3@cnFillOn  = False 
 res3@cnLinesOn = False
 res3@cnLineLabelsOn = False
 res3@cnInfoLabelOn  = False
 ;res3@lbLabelBarOn   = False

 res3@gsnRightString = ""

 res3@cnLevelSelectionMode = "ManualLevels"
 res3@cnMinLevelValF = 0.00
 res3@cnMaxLevelValF = 0.10
 res3@cnLevelSpacingF = 0.01

 ; res3@cnFillPattern = 17 ;(/17/)
 res3@cnFillColor  =  "black"
 res3@gsnLeftString = ""
 res3@cnFillDotSizeF = 0.003
 res3@cnFillScaleF  = 1.2
 ;----------------------------------------------
 ;----------------------------------------------

res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 0.02
res@cnMaxLevelValF = 0.2
res@cnMinLevelValF = -0.2
res@cnFillPalette = "temp_diff_18lev" ; "MPL_RdBu"

res@gsnLeftString = "heat flux diff"
res@gsnRightString = "lat: "+latstr+"~S~o~N~N - "+latend+"~S~o~N~N"
plot(0) = gsn_csm_pres_hgt(wks, diff_dtot(0,:,{lonstr:lonend}), res)

plot_topo(0) = gsn_csm_xy(wks, topo&lon, topo_zonal, res_topo)
overlay(plot(0), plot_topo(0))

draw(plot(0))


end