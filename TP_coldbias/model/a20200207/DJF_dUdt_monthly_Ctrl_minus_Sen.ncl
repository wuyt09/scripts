
begin
lonstr = 80
lonend = 105
latstr = 30
latend = 40
yrStrt = 198102
yrLast = 200601

casename = (/"Ctrl_FAMIP","Hist_FAMIP"/)
ncase = dimsizes(casename)

season = "DJF"
LEV = 200
dUdT = new((/ncase,96,144/), float)

do j = 0, ncase-1
    diri = "/home/yangsong3/data-model/wuyt/TPbias/TPbias_"+casename(j)+"/a20191206/"
    filiu = casename(j)+".cam.h0.U.1979-2005.nc"
    filiv = casename(j)+".cam.h0.V.1979-2005.nc"
    filiw = casename(j)+".cam.h0.OMEGA.1979-2005.nc"

    fu = addfile(diri+filiu, "r")
    fv = addfile(diri+filiv, "r")
    fw = addfile(diri+filiw, "r")

    TIME = fu->time
    YYYY = cd_calendar(TIME, -1)
    iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

    UX = fu->U(iYYYY,:,:,:)
    VX = fv->V(iYYYY,:,:,:)
    WX = fw->OMEGA(iYYYY,:,:,:)

    g_lvl   =   tofloat(UX&level*100) 
    g_lat   =   tofloat(UX&lat)
    g_lon   =   tofloat(UX&lon)
    g_nlat  =   dimsizes(g_lat)
    g_nlon  =   dimsizes(g_lon)
    g_nlvl  =   dimsizes(g_lvl)

    ; U0 = month_to_season(UX, season)
    ; U = dim_avg_n_Wrap(U0, 0)

    ; V0 = month_to_season(VX, season)
    ; V = dim_avg_n_Wrap(V0, 0)

    ; W0 = month_to_season(WX, season)
    ; W = dim_avg_n_Wrap(W0, 0)
    U = UX
    V = VX
    W = WX
    printVarSummary(U)

    ;***********************************************************
    ; du/dt = du/dt + u*du/dx + v*du/dy + w*du/dp
    ;***********************************************************

    ;===================================
    ; du/dt
    ;===================================
    dudt = U
    dudt = center_finite_diff_n(U,86400,False,0,0)

    ;***********TERM [2], this is a little tricky*******
    adv_u_x =  U    ;u*dudx
    adv_u_y =  U    ;v*dudy

    ;u*ds0dx
    dlon = (g_lon(2)-g_lon(1))*0.0174533 ; convert to radians
    
    do i=0,g_nlat-1                      ; loop over each latitude
        dx = 6378388.*cos(0.0174533*g_lat(i))*dlon  ; constant at this latitude
        adv_u_x(:,:,i,:) = center_finite_diff_n (U(:,:,i,:), dx , True,0,2) ; be careful with the dimison downscale
    end do
    adv_u_x =   U*adv_u_x

    ;v*ds0dy
    dlat = (g_lat(2)-g_lat(1))*0.0174533 ; convert to radians
    dy   = 6378388.*dlat   
    adv_u_y = center_finite_diff_n (U, dy , False,0,2)
    adv_u_y = V*adv_u_y
    
    adv_u   =   adv_u_x ; copy metadata
    adv_u   =   adv_u_x + adv_u_y
    ;***************************************************


    ;*****************TERM [3], easy********************
    dudp    =   U
    dudp    =   center_finite_diff_n(U,g_lvl,False,0,1)
    dudp    =   W*dudp
    ;***************************************************

    dUdt      =   dudp ;copy metadata
    dUdt      =   dudt + adv_u + dudp

    dUdt@long_name    =   "monthly dUdt on Pressure Levels"
    dUdt@units    =   "m/s2"
    dUdt@var_desc    =   "dUdt"
    dUdt!1        =   "lev"
    ; print(nyr+" Calculate done!")

    dudt_Sea = month_to_season(dUdt(:,{LEV},:,:), season)
    dudt_sea = dim_avg_n_Wrap(dudt_Sea, 0)

    dUdT(j,:,:) = dudt_sea
; print("     "+i)
; delete([/dlon,dx/])
end do
printVarSummary(dUdT)
dUdT = dUdT*1.e+4

diff_dudt = dUdT(0,:,:)
diff_dudt = dUdT(1,:,:) - dUdT(0,:,:)
printMinMax(dUdT(:,{latstr:latend},{lonstr:lonend}), True)
printMinMax(diff_dudt({latstr:latend},{lonstr:lonend}), True)

;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)

;=========================================================
; plot 
;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/a20200207/"
figname = season+"_dUdt_monthly_Ctrl_minus_Sen"

lonstr = 60
lonend = 120
latstr = 10
latend = 50

wks = gsn_open_wks("pdf", figpath+figname)
plot = new(3, graphic)
plot_topo = plot

 res_topo            = True
 res_topo@cnFillOn   = False
 res_topo@cnLinesOn  = True
 res_topo@gsnFrame   = False
 res_topo@gsnDraw    = False
 res_topo@cnLineLabelsOn         = False
 res_topo@cnLineThicknessF       = 3
 res_topo@cnLineColor            = "grey"
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
res@cnFillPalette         = "BlueWhiteOrangeRed" ;"MPL_RdBu"
;res@lbLabelBarOn          = False

res@mpMinLatF            = latstr      ; range to zoom in on
res@mpMaxLatF            = latend
res@mpMinLonF            = lonstr
res@mpMaxLonF            = lonend
; res@mpCenterLonF         = 180.
res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025

res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -5
res@cnMaxLevelValF = 5
res@cnLevelSpacingF = 0.5
; res@cnMinLevelValF = -2
; res@cnMaxLevelValF = 2
; res@cnLevelSpacingF = 0.2
res@gsnStringFontHeightF = 0.035
;res@gsnRightString       = "K"


res@gsnLeftString = "Ctrl run"
plot(0) = gsn_csm_contour_map(wks, dUdT(0,:,:), res)
plot_topo(0) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(0), plot_topo(0))

res@gsnLeftString = "Sen run"
plot(1) = gsn_csm_contour_map(wks, dUdT(1,:,:), res)
plot_topo(1) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(1), plot_topo(1))


res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -1
res@cnMaxLevelValF = 1
res@cnLevelSpacingF = 0.1
res@gsnLeftString = "diff"
plot(2) = gsn_csm_contour_map(wks, diff_dudt, res)
plot_topo(2) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(2), plot_topo(2))


pres = True
pres@gsnMaximize = True
;pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = season+" mean: dU/dt"
; pres@pmLabelBarOrthogonalPosF = -0.02
; pres@lbLabelFontHeightF = 0.015

gsn_panel(wks, plot, (/1,3/), pres)

end