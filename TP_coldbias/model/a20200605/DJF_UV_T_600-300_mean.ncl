
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
u = new((/ncase,nyear,96,144/), float)
v = u
t = u
levs = 600
levt = 300

do i = 0, ncase-1
    diri = "/home/yangsong3/data-model/wuyt/TPbias/TPbias_"+casename(i)+"/a20191206/"
    filiu = casename(i)+".cam.h0.U.1979-2005.nc"
    filiv = casename(i)+".cam.h0.V.1979-2005.nc"
    filiz = casename(i)+".cam.h0.T.1979-2005.nc"
    filip = casename(i)+".cam.h0.PS.1979-2005.nc"
    
    fu = addfile(diri+filiu, "r")
    fv = addfile(diri+filiv, "r")
    fz = addfile(diri+filiz, "r")
    fp = addfile(diri+filip, "r")

    TIME = fu->time
    YYYY = cd_calendar(TIME, -1)
    iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

    UX = fu->U(iYYYY,:,:,:)
    VX = fv->V(iYYYY,:,:,:)
    TX = fz->T(iYYYY,:,:,:)
    PS = fp->PS(iYYYY,:,:)

    lev = UX&level
    plev = lev
    plev = plev *100
    plev!0 = "lev"
    plev&lev = lev
    plev@units = "Pa"

    ptop = min(plev)
    psfc = PS
    dP = dpres_plevel_Wrap(plev, psfc, ptop, 0)
    ; dP(:,{600},{30},{90}) = dP@_FillValue
    dP = where(ismissing(UX), dP@_FillValue, dP)
    printVarSummary(dP)

    UdP = UX*dP
    copy_VarCoords(UX, UdP)
    VdP = VX*dP
    copy_VarCoords(VX, VdP)
    TdP = TX*dP
    copy_VarCoords(TX, TdP)

    U = dim_sum_n_Wrap(UdP(:,{levs:levt},:,:), 1)/dim_sum_n_Wrap(dP(:,{levs:levt},:,:), 1)
    copy_VarCoords(UX(:,0,:,:), U)
    V = dim_sum_n_Wrap(VdP(:,{levs:levt},:,:), 1)/dim_sum_n_Wrap(dP(:,{levs:levt},:,:), 1)
    copy_VarCoords(VX(:,0,:,:), V)
    T = dim_sum_n_Wrap(TdP(:,{levs:levt},:,:), 1)/dim_sum_n_Wrap(dP(:,{levs:levt},:,:), 1)
    copy_VarCoords(TX(:,0,:,:), T)

    u(i,:,:,:) = month_to_season(U, season)
    v(i,:,:,:) = month_to_season(V, season)
    t(i,:,:,:) = month_to_season(T, season)
end do

u_ave = dim_avg_n_Wrap(u, 1)
u_var = dim_variance_n_Wrap(u, 1)
probu = ttest(u_ave(0,:,:), u_var(0,:,:), nyear, u_ave(1,:,:), u_var(1,:,:), nyear, False, False)
copy_VarCoords(u_ave(0,:,:), probu)

v_ave = dim_avg_n_Wrap(v, 1)
v_var = dim_variance_n_Wrap(v, 1)
probv = ttest(v_ave(0,:,:), v_var(0,:,:), nyear, v_ave(1,:,:), v_var(1,:,:), nyear, False, False)
copy_VarCoords(v_ave(0,:,:), probv)

t_ave = dim_avg_n_Wrap(t, 1)
t_var = dim_variance_n_Wrap(t, 1)
probt = ttest(t_ave(0,:,:), t_var(0,:,:), nyear, t_ave(1,:,:), t_var(1,:,:), nyear, False, False)
copy_VarCoords(t_ave(0,:,:), probt)
printVarSummary(probt)


u_diff = u_ave(1,:,:) - u_ave(0,:,:)
copy_VarCoords(u_ave(0,:,:), u_diff)
v_diff = v_ave(1,:,:) - v_ave(0,:,:)
copy_VarCoords(v_ave(0,:,:), v_diff)
t_diff = t_ave(1,:,:) - t_ave(0,:,:)
copy_VarCoords(t_ave(0,:,:), t_diff)
printVarSummary(t_diff) 


;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)

;=========================================================
; plot 
;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/a20200605/"
figname = season+"_UV_T_600-300_mean"

lonstr = 60
lonend = 150
latstr = 20
latend = 60

wks = gsn_open_wks("pdf", figpath+figname)
plot = new(3, graphic)
plotv = plot
plot_topo = plot
plot_prob = plot

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
res@cnFillPalette         = "CBR_coldhot"
;res@lbLabelBarOn          = False

res@mpMinLatF            = latstr      ; range to zoom in on
res@mpMaxLatF            = latend
res@mpMinLonF            = lonstr
res@mpMaxLonF            = lonend
; res@mpCenterLonF         = 180.
res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025

res@gsnStringFontHeightF = 0.035


vcres = True
vcres@gsnDraw  = False
vcres@gsnFrame = False
vcres@vcMinDistanceF            = 0.025
vcres@vcRefAnnoOrthogonalPosF   = -1.0            ; move ref vector
vcres@vcRefAnnoArrowLineColor   = "black"         ; change ref vector color
;vcres@vcRefAnnoArrowUseVecColor = False
vcres@vcRefAnnoString2On = False
vcres@vcRefAnnoFontHeightF = 0.02

vcres@gsnLeftString = ""
vcres@gsnRightString = ""

vcres@vcGlyphStyle            = "LineArrow"     ; turn on curly vectors CurlyVector
;vcres@vcLineArrowColor        = "white"           ; change vector color
vcres@vcLineArrowThicknessF   = 3.0               ; change vector thickness
;vcres@vcVectorDrawOrder       = "PostDraw"        ; draw vectors last


res@gsnRightString = "K"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = 220
res@cnMaxLevelValF = 260
res@cnLevelSpacingF = 2.
res@gsnLeftString = "Ctrl"
vcres@vcRefMagnitudeF           = 20.0             ; define vector ref mag
vcres@vcRefLengthF              = 0.045           ; define length of vec ref
plot(0) = gsn_csm_contour_map(wks, t_ave(0,:,:), res)
plotv(0) = gsn_csm_vector(wks, u_ave(0,:,:), v_ave(0,:,:), vcres)
overlay(plot(0), plotv(0))
plot_topo(0) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(0), plot_topo(0))

res@gsnLeftString = "Sen"
plot(1) = gsn_csm_contour_map(wks, t_ave(1,:,:), res)
plotv(1) = gsn_csm_vector(wks, u_ave(1,:,:), v_ave(1,:,:), vcres)
overlay(plot(1), plotv(1))
plot_topo(1) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(1), plot_topo(1))

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
 res3@cnFillColor  =  "red"
 res3@gsnLeftString = ""
 res3@cnFillDotSizeF = 0.003
 res3@cnFillScaleF  = 1.2
 ;----------------------------------------------
 ;----------------------------------------------

res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -2
res@cnMaxLevelValF = 2
res@cnLevelSpacingF = 0.2
res@gsnLeftString = "Diff"
vcres@vcRefMagnitudeF           = 2.0             ; define vector ref mag
vcres@vcRefLengthF              = 0.045           ; define length of vec ref
plot(2) = gsn_csm_contour_map(wks, t_diff, res)
plotv(2) = gsn_csm_vector(wks, u_diff, v_diff, vcres)
overlay(plot(2), plotv(2))
plot_topo(2) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(2), plot_topo(2))
plot_prob(2) = gsn_csm_contour(wks, probt, res3)
plot_prob(2) = ShadeLtContour(plot_prob(2),0.055,17)
overlay(plot(2), plot_prob(2))

pres = True
pres@gsnMaximize = True
;pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = season+" mean UV & T (600hPa-300hPa) " ;(500hPa to 200hPa)
; pres@pmLabelBarOrthogonalPosF = -0.02
; pres@lbLabelFontHeightF = 0.015

gsn_panel(wks, plot, (/3,1/), pres)

end