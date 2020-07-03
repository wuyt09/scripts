
begin

lonstr = 75
lonend = 105
latstr = 25
latend = 40
yrStrt = 198102
yrLast = 200601

casename = (/"Ctrl_FAMIP","Hist_FAMIP"/)
ncase = dimsizes(casename)

season = "DJF"
nyear = 25
U_X = new((/ncase,nyear,17,144/), float)
U_Y = new((/ncase,nyear,17,96/), float)
V_X = U_X
V_Y = U_Y
W_X = U_X
W_Y = U_Y
Z_X = U_X
Z_Y = U_Y
DV_X = U_X
DV_Y = U_Y
dudt_X = U_X
dudt_Y = U_Y

do i = 0, ncase-1
    diri = "/home/yangsong3/data-model/wuyt/TPbias/TPbias_"+casename(i)+"/a20191206/"
    filiw = casename(i)+".cam.h0.OMEGA.1979-2005.nc"
    filiu = casename(i)+".cam.h0.U.1979-2005.nc"
    filiv = casename(i)+".cam.h0.V.1979-2005.nc"
    filiz = casename(i)+".cam.h0.Z3.1979-2005.nc"
    
    fu = addfile(diri+filiu, "r")
    fv = addfile(diri+filiv, "r")
    fw = addfile(diri+filiw, "r")
    fz = addfile(diri+filiz, "r")

    TIME = fu->time
    YYYY = cd_calendar(TIME, -1)
    iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

    OMEGAX = fw->OMEGA(iYYYY,:,:,:)
    OMEGA = month_to_season(OMEGAX, season)
    OMEGA = OMEGA*100.
    OMEGA@units = "10^-2 Pa/s"    
    W_X(i,:,:,:) = dim_avg_n_Wrap(OMEGA(:,:,{latstr:latend},:), 2)
    W_Y(i,:,:,:) = dim_avg_n_Wrap(OMEGA(:,:,:,{lonstr:lonend}), 3)

    UX = fu->U(iYYYY,:,:,:)
    U = month_to_season(UX, season)
    U_X(i,:,:,:) = dim_avg_n_Wrap(U(:,:,{latstr:latend},:), 2)
    U_Y(i,:,:,:) = dim_avg_n_Wrap(U(:,:,:,{lonstr:lonend}), 3)

    VX = fv->V(iYYYY,:,:,:)
    V = month_to_season(VX, season)
    V_X(i,:,:,:) = dim_avg_n_Wrap(V(:,:,{latstr:latend},:), 2)
    V_Y(i,:,:,:) = dim_avg_n_Wrap(V(:,:,:,{lonstr:lonend}), 3)

    ZX = fz->Z3(iYYYY,:,:,:)
    Z = month_to_season(ZX, season)
    Z@units = "gpm"    
    Z_X(i,:,:,:) = dim_avg_n_Wrap(Z(:,:,{latstr:latend},:), 2)
    Z_Y(i,:,:,:) = dim_avg_n_Wrap(Z(:,:,:,{lonstr:lonend}), 3)

    DVX = uv2dv_cfd(UX, VX, UX&lat, UX&lon, 1)
    copy_VarCoords(UX, DVX)
    DVX = DVX * 1.e+6
    printVarSummary(DVX)
    DV = month_to_season(DVX, season)
    DV_X(i,:,:,:) = dim_avg_n_Wrap(DV(:,:,{latstr:latend},:), 2)
    DV_Y(i,:,:,:) = dim_avg_n_Wrap(DV(:,:,:,{lonstr:lonend}), 3)



    ;***********************************************************
    ; du/dt = du/dt + u*du/dx + v*du/dy + w*du/dp
    ;***********************************************************
    g_lvl   =   tofloat(UX&level*100) 
    g_lat   =   tofloat(UX&lat)
    g_lon   =   tofloat(UX&lon)
    g_nlat  =   dimsizes(g_lat)
    g_nlon  =   dimsizes(g_lon)
    g_nlvl  =   dimsizes(g_lvl)
    ;===================================
    ; du/dt
    ;===================================
    dudt = UX
    dudt = center_finite_diff_n(UX,86400,False,0,0)

    ;***********TERM [2], this is a little tricky*******
    adv_u_x =  UX    ;u*dudx
    adv_u_y =  UX    ;v*dudy

    ;u*ds0dx
    dlon = (g_lon(2)-g_lon(1))*0.0174533 ; convert to radians
    
    do k=0,g_nlat-1                      ; loop over each latitude
        dx = 6378388.*cos(0.0174533*g_lat(k))*dlon  ; constant at this latitude
        adv_u_x(:,:,k,:) = center_finite_diff_n (UX(:,:,k,:), dx , True,0,2) ; be careful with the dimison downscale
    end do
    adv_u_x =   UX*adv_u_x

    ;v*ds0dy
    dlat = (g_lat(2)-g_lat(1))*0.0174533 ; convert to radians
    dy   = 6378388.*dlat   
    adv_u_y = center_finite_diff_n (UX, dy , False,0,2)
    adv_u_y = VX*adv_u_y
    
    adv_u   =   adv_u_x ; copy metadata
    adv_u   =   adv_u_x + adv_u_y
    ;***************************************************


    ;*****************TERM [3], easy********************
    dudp    =   UX
    dudp    =   center_finite_diff_n(UX,g_lvl,False,0,1)
    dudp    =   OMEGAX*dudp
    ;***************************************************

    dUdt      =   dudp ;copy metadata
    dUdt      =   dudt + adv_u + dudp
    dUdt      =   dUdt*1.e+5

    dUdt@long_name    =   "monthly dUdt on Pressure Levels"
    dUdt@units    =   "m/s2"
    dUdt@var_desc    =   "dUdt"
    dUdt!1        =   "lev"
    ; print(nyr+" Calculate done!")

    dudt_Sea = month_to_season(dUdt, season)    
    dudt_X(i,:,:,:) = dim_avg_n_Wrap(dudt_Sea(:,:,{latstr:latend},:), 2)
    dudt_Y(i,:,:,:) = dim_avg_n_Wrap(dudt_Sea(:,:,:,{lonstr:lonend}), 3)
end do
printVarSummary(U_X)

;=======================================================================
; ttest
;=======================================================================

ux_ave = dim_avg_n_Wrap(U_X, 1)
ux_var = dim_variance_n_Wrap(U_X, 1)
probux = ttest(ux_ave(0,:,:), ux_var(0,:,:), nyear, ux_ave(1,:,:), ux_var(1,:,:), nyear, False, False)
copy_VarCoords(ux_ave(0,:,:), probux)
printVarSummary(probux)
uy_ave = dim_avg_n_Wrap(U_Y, 1)
uy_var = dim_variance_n_Wrap(U_Y, 1)
probuy = ttest(uy_ave(0,:,:), uy_var(0,:,:), nyear, uy_ave(1,:,:), uy_var(1,:,:), nyear, False, False)
copy_VarCoords(uy_ave(0,:,:), probuy)
printVarSummary(probuy)


vx_ave = dim_avg_n_Wrap(V_X, 1)
vx_var = dim_variance_n_Wrap(V_X, 1)
probvx = ttest(vx_ave(0,:,:), vx_var(0,:,:), nyear, vx_ave(1,:,:), vx_var(1,:,:), nyear, False, False)
copy_VarCoords(vx_ave(0,:,:), probvx)
printVarSummary(probvx)
vy_ave = dim_avg_n_Wrap(V_Y, 1)
vy_var = dim_variance_n_Wrap(V_Y, 1)
probvy = ttest(vy_ave(0,:,:), vy_var(0,:,:), nyear, vy_ave(1,:,:), vy_var(1,:,:), nyear, False, False)
copy_VarCoords(vy_ave(0,:,:), probvy)
printVarSummary(probvy)

wx_ave = dim_avg_n_Wrap(W_X, 1)
wx_var = dim_variance_n_Wrap(W_X, 1)
probwx = ttest(wx_ave(0,:,:), wx_var(0,:,:), nyear, wx_ave(1,:,:), wx_var(1,:,:), nyear, False, False)
copy_VarCoords(wx_ave(0,:,:), probwx)
printVarSummary(probwx)
wy_ave = dim_avg_n_Wrap(W_Y, 1)
wy_var = dim_variance_n_Wrap(W_Y, 1)
probwy = ttest(wy_ave(0,:,:), wy_var(0,:,:), nyear, wy_ave(1,:,:), wy_var(1,:,:), nyear, False, False)
copy_VarCoords(wy_ave(0,:,:), probwy)
printVarSummary(probwy)

zx_ave = dim_avg_n_Wrap(Z_X, 1)
zx_var = dim_variance_n_Wrap(Z_X, 1)
probzx = ttest(zx_ave(0,:,:), zx_var(0,:,:), nyear, zx_ave(1,:,:), zx_var(1,:,:), nyear, False, False)
copy_VarCoords(zx_ave(0,:,:), probzx)
printVarSummary(probzx)
zy_ave = dim_avg_n_Wrap(Z_Y, 1)
zy_var = dim_variance_n_Wrap(Z_Y, 1)
probzy = ttest(zy_ave(0,:,:), zy_var(0,:,:), nyear, zy_ave(1,:,:), zy_var(1,:,:), nyear, False, False)
copy_VarCoords(zy_ave(0,:,:), probzy)
printVarSummary(probzy)

dudtx_ave = dim_avg_n_Wrap(dudt_X, 1)
dudtx_var = dim_variance_n_Wrap(dudt_X, 1)
probdudtx = ttest(dudtx_ave(0,:,:), dudtx_var(0,:,:), nyear, dudtx_ave(1,:,:), dudtx_var(1,:,:), nyear, False, False)
copy_VarCoords(dudtx_ave(0,:,:), probdudtx)
printVarSummary(probdudtx)
dudty_ave = dim_avg_n_Wrap(dudt_Y, 1)
dudty_var = dim_variance_n_Wrap(dudt_Y, 1)
probdudty = ttest(dudty_ave(0,:,:), dudty_var(0,:,:), nyear, dudty_ave(1,:,:), dudty_var(1,:,:), nyear, False, False)
copy_VarCoords(dudty_ave(0,:,:), probdudty)
printVarSummary(probdudty)

dvx_ave = dim_avg_n_Wrap(DV_X, 1)
dvx_var = dim_variance_n_Wrap(DV_X, 1)
probdvx = ttest(dvx_ave(0,:,:), dvx_var(0,:,:), nyear, dvx_ave(1,:,:), dvx_var(1,:,:), nyear, False, False)
copy_VarCoords(dvx_ave(0,:,:), probdvx)
printVarSummary(probdvx)
dvy_ave = dim_avg_n_Wrap(DV_Y, 1)
dvy_var = dim_variance_n_Wrap(DV_Y, 1)
probdvy = ttest(dvy_ave(0,:,:), dvy_var(0,:,:), nyear, dvy_ave(1,:,:), dvy_var(1,:,:), nyear, False, False)
copy_VarCoords(dvy_ave(0,:,:), probdvy)
printVarSummary(probdvy)


diff_wx = wx_ave(0,:,:)
diff_wx = wx_ave(1,:,:) - wx_ave(0,:,:)
diff_wy = wy_ave(0,:,:)
diff_wy = wy_ave(1,:,:) - wy_ave(0,:,:)

diff_ux = ux_ave(0,:,:)
diff_ux = ux_ave(1,:,:) - ux_ave(0,:,:)
diff_uy = uy_ave(0,:,:)
diff_uy = uy_ave(1,:,:) - uy_ave(0,:,:)

diff_vx = vx_ave(0,:,:)
diff_vx = vx_ave(1,:,:) - vx_ave(0,:,:)
diff_vy = vy_ave(0,:,:)
diff_vy = vy_ave(1,:,:) - vy_ave(0,:,:)

diff_zy = zy_ave(0,:,:)
diff_zy = zy_ave(1,:,:) - zy_ave(0,:,:)
diff_zx = zx_ave(0,:,:)
diff_zx = zx_ave(1,:,:) - zx_ave(0,:,:)

diff_dudty = dudty_ave(0,:,:)
diff_dudty = dudty_ave(1,:,:) - dudty_ave(0,:,:)
diff_dudtx = dudtx_ave(0,:,:)
diff_dudtx = dudtx_ave(1,:,:) - dudtx_ave(0,:,:)

diff_dvy = dvy_ave(0,:,:)
diff_dvy = dvy_ave(1,:,:) - dvy_ave(0,:,:)
diff_dvx = dvx_ave(0,:,:)
diff_dvx = dvx_ave(1,:,:) - dvx_ave(0,:,:)

;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)

;=========================================================
; plot 
;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/a20200408/"
figname = season+"_vertical_cell_U_dUdt_wind_divergence"

; lonstr = 60
; lonend = 120
; latstr = 10
; latend = 50

wks = gsn_open_wks("pdf", figpath+figname)
plot = new(3, graphic)
plotc = plot
plotd = plot
plot_topo = plot
plot_prob = plot

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
    ;res@lbLabelBarOn = False
    ;res@lbLabelFontHeightF = 0.02
    res@tiYAxisString = ""

    res@tmYLMode      = "Explicit"
    res@tmYLValues    = (/"1000","850","750","500","300","200","100"/)
    res@tmYLLabels    = (/"1000","850","750","500","300","200","100"/)

 ;------------------------------------------
     res@trYMinF   = 100
     res@trYMaxF   = 1000
         
     res@vcRefAnnoOrthogonalPosF =  -0.35 ;-1.0             ; move ref vector up
     ;res@vcGlyphStyle            = "CurlyVector"  ;"LineArrow"    ; turn on curly vectors
     ;res@vcMinDistanceF          = 0.017
     res@vcGlyphStyle            = "LineArrow"
     res@vcMinDistanceF          = 0.03
     res@vcRefLengthF            = 0.04            ; define length of vec ref
     
     res@vcLineArrowThicknessF   = 3.0
     ;res@vcLineArrowHeadMinSizeF = 0.008
     ;res@vcLineArrowHeadMaxSizeF = 0.3
     res@vcVectorDrawOrder = "PostDraw" 
     
     res@vcRefAnnoString2On = False
     res@vcRefAnnoFontHeightF = 0.02
     res@gsnPresHgtHeightLabelOn = False
     res@tmYROn = False

;=====================================================

    cnres = True
    cnres@gsnDraw = False
    cnres@gsnFrame = False

    cnres@cnFillOn = False
    cnres@cnLinesOn = True
    cnres@cnLineLabelsOn = True
    cnres@cnLineColor = "green4"
    cnres@cnLineThicknessF = 5
    cnres@cnInfoLabelOn = False
    ;cnres@cnLineDrawOrder = "PostDraw"
    cnres@cnLabelDrawOrder = "PostDraw"
    ;cnres@cnLevelFlag = "LineAndLabel"
    cnres@cnLineLabelInterval = 1
    cnres@gsnLeftString = ""
    cnres@gsnRightString = ""
;-----------------------------------------------------
cnres@gsnContourPosLineDashPattern = 0
cnres@gsnContourNegLineDashPattern = 2

res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 0.5
res@cnMaxLevelValF = 5.
res@cnMinLevelValF = -5.
res@cnFillPalette = "BlueWhiteOrangeRed" ; "MPL_RdBu"
res@vcRefMagnitudeF         = 10.             ; define vector ref mag

; cnres@cnLevelSelectionMode = "ExplicitLevels"
; cnres@cnLevels = (/-10,-7,-4,-2,-0.5,0.5,2,4,7,10/)

res@gsnLeftString = "Ctrl"; (90~S~o~N~E)
plot(0) = gsn_csm_pres_hgt_vector(wks, dvy_ave(0,:,{-20:60}), uy_ave(0,:,{-20:60}), -wy_ave(0,:,{-20:60}), res)

cnres@cnLineColor = "red"
cnres@cnLineLabelFontColor = "red"
cnres@cnLevelSelectionMode = "ManualLevels"
cnres@cnLevelSpacingF = 10.
cnres@cnMaxLevelValF = 50.
cnres@cnMinLevelValF = -50.
plotc(0) = gsn_csm_contour(wks, uy_ave(0,:,{-20:60}), cnres)
overlay(plot(0), plotc(0))

cnres@cnLineColor = "green4"
cnres@cnLineLabelFontColor = "green4"
cnres@cnLevelSelectionMode = "ManualLevels"
cnres@cnLevelSpacingF = 10.
cnres@cnMaxLevelValF = 100.
cnres@cnMinLevelValF = -100.
plotd(0) = gsn_csm_contour(wks, dudty_ave(0,:,{-20:60}), cnres)
overlay(plot(0), plotd(0))

res@gsnLeftString = "Sen"; (90~S~o~N~E)
plot(1) = gsn_csm_pres_hgt_vector(wks, dvy_ave(1,:,{-20:60}), uy_ave(1,:,{-20:60}), -wy_ave(1,:,{-20:60}), res)
cnres@cnLineColor = "red"
cnres@cnLineLabelFontColor = "red"
cnres@cnLevelSelectionMode = "ManualLevels"
cnres@cnLevelSpacingF = 10.
cnres@cnMaxLevelValF = 50.
cnres@cnMinLevelValF = -50.
plotc(1) = gsn_csm_contour(wks, uy_ave(1,:,{-20:60}), cnres)
overlay(plot(1), plotc(1))

cnres@cnLineColor = "green4"
cnres@cnLineLabelFontColor = "green4"
cnres@cnLevelSelectionMode = "ManualLevels"
cnres@cnLevelSpacingF = 10.
cnres@cnMaxLevelValF = 100.
cnres@cnMinLevelValF = -100.
plotd(1) = gsn_csm_contour(wks, dudty_ave(1,:,{-20:60}), cnres)
overlay(plot(1),plotd(1))
;==========================================================================
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
 res3@cnFillColor  =  "purple"
 res3@gsnLeftString = ""
 res3@cnFillDotSizeF = 0.003
 res3@cnFillScaleF  = 1.2
 ;----------------------------------------------
 ;----------------------------------------------

res@cnFillPalette = "BlueWhiteOrangeRed" ; "MPL_RdBu"
res@vcRefMagnitudeF         = 2.             ; define vector ref mag

cnres@cnLevelSelectionMode = "ManualLevels"
cnres@cnLevelSpacingF = 1.
cnres@cnMaxLevelValF = 5.
cnres@cnMinLevelValF = -5.

res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 0.1
res@cnMaxLevelValF = 1.0
res@cnMinLevelValF = -1.0

res@gsnLeftString = "Diff"; (90~S~o~N~E)
; diff_dvy = mask(diff_dvy, probdvy.le.0.055, True)
; diff_wy = mask(diff_wy, probwy.le.0.055, True)
; diff_uy = mask(diff_uy, probuy.le.0.055, True)
plot(2) = gsn_csm_pres_hgt_vector(wks, diff_dvy(:,{-20:60}), diff_uy(:,{-20:60}), -diff_wy(:,{-20:60}), res)
cnres@cnLineColor = "red"
cnres@cnLineLabelFontColor = "red"
cnres@cnLevelSelectionMode = "ManualLevels"
cnres@cnLevelSpacingF = 2.
cnres@cnMaxLevelValF = 10.
cnres@cnMinLevelValF = -10.
plotc(2) = gsn_csm_contour(wks, diff_uy(:,{-20:60}), cnres)
overlay(plot(2), plotc(2))
cnres@cnLineColor = "green4"
cnres@cnLineLabelFontColor = "green4"
cnres@cnLevelSelectionMode = "ManualLevels"
cnres@cnLevelSpacingF = 2.
cnres@cnMaxLevelValF = 10.
cnres@cnMinLevelValF = -10.
plotd(2) = gsn_csm_contour(wks, diff_dudty(:,{-20:60}), cnres)
overlay(plot(2), plotd(2))

; plot_prob(1) = gsn_csm_contour(wks,probzy,res3)
; plot_prob(1) = ShadeLtContour(plot_prob(1),0.055,17)
; overlay(plot(2), plot_prob(1))


pres = True
pres@gsnMaximize = True
;pres@gsnPanelLabelBar  = True
pres@gsnPanelMainString = season+" mean: wind Div (shading) & U (red contour) & dUdt (green contour), 75~S~o~N~E - 105~S~o~N~E"
gsn_panel(wks, plot, (/1,3/), pres)

end