
begin

runs = (/"CTRL", "NOQ1","NOTP_NOQ1"/)
nm = dimsizes(runs)
u = new((/nm,96,144/), float)
u!0 = "run"
u&run = runs
v = u
t = u
z = new((/nm,37,96,144/), float)
z!0 = "run"
z&run = runs

do i = 0, nm-1

    diri = "/home/yangsong3/data-model/wuyt/TP_"+runs(i)+"_F2000_F19G16_CAM4/pro/"

    filiu = "F2000_"+runs(i)+".plevel.U.mon.0601-3012.nc"
    filiv = "F2000_"+runs(i)+".plevel.V.mon.0601-3012.nc"
    filiz = "F2000_"+runs(i)+".plevel.Z3.mon.0601-3012.nc"
    filip = "F2000_"+runs(i)+".PS.mon.0601-3012.nc"
    filit = "F2000_"+runs(i)+".plevel.T.mon.0601-3012.nc"


    fu = addfile(diri+filiu, "r")
    fv = addfile(diri+filiv, "r")
    fz = addfile(diri+filiz, "r")
    fp = addfile(diri+filip, "r")
    ft = addfile(diri+filit, "r")

    lat = fu->lat
    lon = fu->lon
    lev = fu->lev

    U = fu->U(:,{200},:,:)
    V = fv->V(:,{200},:,:)
    Z = fz->Z3
    PS = fp->PS
    T = ft->T
    printVarSummary(PS)

    plev = T&lev
    plev&lev = T&lev
    psfc = max(plev) ;PS ;
    ptop = min(plev)

    dP = dpres_plevel_Wrap(plev, psfc, ptop, 0)
    printVarSummary(dP)
    dp = conform_dims(dimsizes(T), dP, 1)
    Tdp = T*dp
    copy_VarCoords(T, Tdp)
    TT = dim_sum_n_Wrap(Tdp(:,{500:200},:,:), 1)/sum(dP({500:200}))
    copy_VarCoords(T(:,0,:,:), TT)
  
    season = "SON"

    ux = month_to_season(U, season)
    vx = month_to_season(V, season)
    zx = month_to_season(Z, season)
    tx = month_to_season(TT, season)
    u(i,:,:) = dim_avg_n_Wrap(ux, 0)
    v(i,:,:) = dim_avg_n_Wrap(vx, 0)
    z(i,:,:,:) = dim_avg_n_Wrap(zx, 0)
    t(i,:,:) = dim_avg_n_Wrap(tx, 0)

delete([/ux,vx,zx,tx/])

end do
printVarSummary(u)


tpmaskpath = "/home/yangsong3/wuyt/sysu/scripts/TP_CESM_F/model_need_input_data/"
tpmaskname = "TPmask_1000m_C20190410.nc"
ftp = addfile(tpmaskpath+tpmaskname, "r")

tpmask = ftp->tp
printVarSummary(tpmask)

;=========================================================
; plot
;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/JJA/"
figname = "draw_clim_JJA_UV200_SAH_TropoT_CTRL_S1_S2"
wks = gsn_open_wks("pdf", figpath+figname)
plot = new(nm, graphic)
plotv = plot
plot_topo = plot
plotc = plot

lonstr = 20 ;20
lonend = 160 ;220
latstr = -20
latend = 70

 res_topo            = True
 res_topo@cnFillOn   = False
 res_topo@cnLinesOn  = True
 res_topo@gsnFrame   = False
 res_topo@gsnDraw    = False
 res_topo@cnLineLabelsOn         = False
 res_topo@cnLineThicknessF       = 3
 res_topo@cnLineColor            = "blue"
 res_topo@cnInfoLabelOn          = False
 res_topo@cnLevelSelectionMode   ="ExplicitLevels"
 res_topo@cnLevels               = (/0.9/)
 res_topo@gsnLeftString = ""
 res_topo@gsnRightString = ""


res   = True
res@gsnDraw  = False
res@gsnFrame = False
;res@gsnAddCyclic = False

res@cnFillOn      = True
res@cnLinesOn     = False
res@cnInfoLabelOn = False
res@cnLineLabelsOn = False

res@mpOutlineOn   = True
res@mpFillOn      = False
res@mpMinLatF     = latstr
res@mpMaxLatF     = latend
res@mpMinLonF     = lonstr
res@mpMaxLonF     = lonend

res@gsnStringFontHeightF = 0.04

res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025

;res@lbLabelFontHeightF   = 0.03
;res@pmLabelBarHeightF    = 0.1
;res@pmLabelBarWidthF     = 0.65
res@pmLabelBarOrthogonalPosF = 0.12
;res@o=pmLabelBarParallelPosF = 


;===========================================================
cnres = True
cnres@gsnDraw  = False
cnres@gsnFrame = False

cnres@cnFillOn = False
cnres@cnLinesOn = True
cnres@cnLineLabelsOn = True
cnres@cnLineLabelFontColor = "red"
cnres@cnLineColor = "green4"
cnres@cnLineThicknessF = 3
cnres@cnInfoLabelOn = False
;cnres@cnLineDrawOrder = "PostDraw"
cnres@cnLabelDrawOrder = "PostDraw"
;cnres@cnLevelFlag = "LineAndLabel"
cnres@cnLineLabelInterval = 1

cnres@cnLevelSelectionMode = "ManualLevels"
cnres@cnLevelSpacingF = 20.
cnres@cnMaxLevelValF = 12550.
cnres@cnMinLevelValF = 12450.

;===========================================================

;===========================================================
 vecres  = True
 vecres@gsnDraw  = False
 vecres@gsnFrame = False

 vecres@vcGlyphStyle     = "LineArrow"   ; curly vectors"Fillarrow" 
 vecres@vcRefLengthF     = 0.03            ; define length of vec ref
 vecres@vcMinDistanceF   = 0.03
 vecres@vcLineArrowThicknessF = 3.0
 vecres@tiXAxisString    = " "             ; turn off axis label
 vecres@vcRefAnnoOn      = True            ; draw the reference vector annotation.
 vecres@vcRefAnnoString2On =  False          ;display a string below or to the right
 vecres@vcRefAnnoOrientation = "horizontal"
 vecres@vcRefAnnoString2  = " "
 vecres@vcRefAnnoOrthogonalPosF = -1.0  ; move ref vector into plot
 vecres@vcFillArrowsOn       = True
 vecres@vcLineArrowHeadMinSizeF = 0.008
 vecres@vcLineArrowHeadMaxSizeF = 0.3
 vecres@gsnLeftString  = ""
 ;vecres@gsnRightString  = "g*s/Kg"
 vecres@gsnRightString  = ""
 vecres@gsnRightStringFontHeightF = 0.04
 
 vecres@vcLineArrowColor = "black"
 vecres@vcRefAnnoFontHeightF = 0.02
 vecres@vcRefMagnitudeF  = 10.             ; define vector ref mag
 vecres@vcRefAnnoString1  = "10 m/s"
;=========================================================== 

res@gsnRightString = "gpm"
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 0.5
res@cnMaxLevelValF = 247. ;255
res@cnMinLevelValF = 240. ;245
res@cnFillPalette = "sunshine_9lev" ;BlueWhiteOrangeRed
res@lbLabelBarOn = False

do i = 0, nm-1
    res@gsnLeftString = runs(i)
    plot(i) = gsn_csm_contour_map(wks, t(i,:,:), res)
    plot_topo(i) = gsn_csm_contour(wks, tpmask, res_topo)
    plotv(i) = gsn_csm_vector(wks, u(i,:,:), v(i,:,:), vecres)
    plotc(i) = gsn_csm_contour(wks, z(i,{200},:,:), cnres)
    overlay(plot(i), plotv(i))
    overlay(plot(i), plotc(i))
    overlay(plot(i), plot_topo(i))
end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
gsn_panel(wks, plot, (/1,3/), pres)


t_th = t(0,:,:)
t_me = t_th
t_th = t(0,:,:)-t(1,:,:)
t_me = t(1,:,:)-t(2,:,:)

z_th = z(0,:,:,:)
z_me = z_th
z_th = z(0,:,:,:)-z(1,:,:,:)
z_me = z(1,:,:,:)-z(2,:,:,:)

u_th = u(0,:,:)
u_me = u_th
u_th = u(0,:,:)-u(1,:,:)
u_me = u(1,:,:)-u(2,:,:)

v_th = v(0,:,:)
v_me = v_th
v_th = v(0,:,:)-v(1,:,:)
v_me = v(1,:,:)-v(2,:,:)


res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 0.5
res@cnMaxLevelValF = 5.
res@cnMinLevelValF = -5.
res@cnFillPalette = "BlueWhiteOrangeRed"
vecres@vcRefMagnitudeF         = 3.0             ; define vector ref mag
vecres@vcRefLengthF            = 0.025            ; define length of vec ref
cnres@cnLevelSelectionMode = "ManualLevels"
cnres@cnLevelSpacingF = 10.
cnres@cnMaxLevelValF = 80.
cnres@cnMinLevelValF = -80.

res@gsnLeftString = "thermal"
plot(0) = gsn_csm_contour_map(wks, t_th, res)
plot_topo(0) = gsn_csm_contour(wks, tpmask, res_topo)
plotv(0) = gsn_csm_vector(wks, u_th, v_th, vecres)
plotc(0) = gsn_csm_contour(wks, z_th({200},:,:), cnres)
overlay(plot(0), plot_topo(0))
overlay(plot(0), plotv(0))
overlay(plot(0), plotc(0))

res@gsnLeftString = "mechanical"
plot(1) = gsn_csm_contour_map(wks, t_me, res)
plot_topo(1) = gsn_csm_contour(wks, tpmask, res_topo)
plotv(1) = gsn_csm_vector(wks, u_me, v_me, vecres)
plotc(1) = gsn_csm_contour(wks, z_me({200},:,:), cnres)
overlay(plot(1), plot_topo(1))
overlay(plot(1), plotv(1))
overlay(plot(1), plotc(1))

gsn_panel(wks, plot(0:1), (/1,2/), pres)

end