
begin

runs = (/"CTRL", "NOQ1","NOTP_NOQ1"/)
nm = dimsizes(runs)
qu = new((/nm,96,144/), float)
qu!0 = "run"
qu&run = runs
qv = new((/nm,96,144/), float)
qdiv = new((/nm,96,144/), float)

do i = 0, nm-1

    diri = "/home/yangsong3/data-model/wuyt/TP_"+runs(i)+"_F2000_F19G16_CAM4/pro/"

    filiu = "F2000_"+runs(i)+".plevel.U.mon.0601-3012.nc"
    filiv = "F2000_"+runs(i)+".plevel.V.mon.0601-3012.nc"
    filiq = "F2000_"+runs(i)+".plevel.Q.mon.0601-3012.nc"
    filip = "F2000_"+runs(i)+".PS.mon.0601-3012.nc"

    fu = addfile(diri+filiu, "r")
    fv = addfile(diri+filiv, "r")
    fq = addfile(diri+filiq, "r")
    fp = addfile(diri+filip, "r")

    lat = fu->lat
    lon = fu->lon
    lev = fu->lev

    U = fu->U
    V = fv->V
    Q = fq->Q
    PS = fp->PS
    printVarSummary(PS)
    ; print(U(6,{850},{20:30},{110:120}))
    ; print(V(6,{850},{20:30},{110:120}))
    ; print(Q(6,{850},{20:30},{110:120}))


    QU = Q
    QU = Q*U
    QV = Q
    QV = Q*V

    ;==================================================
    ; column integrated
    ;==================================================
    plev = lev
    plev = lev*100
    plev@units = "Pa"
    printVarSummary(plev)
    psfc = PS ; units Pa
    ptop = min(plev)
    dp = dpres_plevel_Wrap(plev, psfc, ptop, 0)

    qudp = QU
    qudp = QU*dp
    qvdp = QV
    qvdp = QV*dp
    qu_lay = dim_sum_n_Wrap(qudp, 1)
    qv_lay = dim_sum_n_Wrap(qvdp, 1)
    qu_lay = qu_lay/9.8
    qv_lay = qv_lay/9.8  ; units Kg/(m*s)
    copy_VarCoords(Q(:,0,:,:), qu_lay)
    copy_VarCoords(Q(:,0,:,:), qv_lay)
    printVarSummary(qu_lay)

    season = "DJF"

    Qu = month_to_season(qu_lay, season)
    Qv = month_to_season(qv_lay, season)
    qu(i,:,:) = dim_avg_n_Wrap(Qu, 0)
    qv(i,:,:) = dim_avg_n_Wrap(Qv, 0)
    
    ; QQ_div = QQ_div*1000.
    ; copy_VarCoords(qu_lay, QQ_div)
    ; Qq_div = month_to_season(QQ_div, season)
    ; qdiv(i,:,:) = dim_avg_n_Wrap(Qq_div, 0)
end do
qdiv = uv2dv_cfd(qu, qv, lat, lon, 1)
qdiv = qdiv*1000
copy_VarCoords(qu, qdiv)

; print(qu(0,{20:30},{110:120}))
; print(qv(0,{20:30},{110:120}))
; print(qdiv(0,{20:30},{110:120}))

tpmaskpath = "/home/yangsong3/wuyt/sysu/scripts/TP_CESM_F/model_need_input_data/"
tpmaskname = "TPmask_1000m_C20190410.nc"
ftp = addfile(tpmaskpath+tpmaskname, "r")

tpmask = ftp->tp
printVarSummary(tpmask)

;=========================================================
; plot
;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/JJA/"
figname = "draw_clim_JJA_WVflux_divergence_CTRL_S1_S2"
wks = gsn_open_wks("pdf", figpath+figname)
plot = new(nm, graphic)
plotv = plot
plot_topo = plot
plot2 = plot

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
res@cnFillPalette = "CBR_drywet"

;===========================================================
 vecres  = True
 vecres@gsnDraw  = False
 vecres@gsnFrame = False

 vecres@vcGlyphStyle     = "LineArrow"   ; curly vectors"Fillarrow" 
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
 vecres@vcRefMagnitudeF  = 200.             ; define vector ref mag
 vecres@vcRefAnnoString1  = "200 Kg/m/s"
 vecres@vcRefLengthF     = 0.04            ; define length of vec ref
 vecres@vcMinDistanceF   = 0.03
;=========================================================== 

res@gsnSpreadColors = True
res@gsnSpreadColorEnd   = 0
res@gsnSpreadColorStart = 10
res@gsnRightString = "g/m~S~2~N~/s"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMaxLevelValF  = 0.16
res@cnMinLevelValF  = -0.16
res@cnLevelSpacingF = 0.02
res@lbLabelBarOn = False

do i = 0, nm-1
    res@gsnLeftString = runs(i)
    plot(i) = gsn_csm_contour_map(wks, qdiv(i,:,:), res)
    plot_topo(i) = gsn_csm_contour(wks, tpmask, res_topo)
    plotv(i) = gsn_csm_vector(wks, qu(i,:,:), qv(i,:,:), vecres)
    overlay(plot(i), plotv(i))
    overlay(plot(i), plot_topo(i))
end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
gsn_panel(wks, plot, (/1,3/), pres)

qdiv_th = qdiv(0,:,:)
qdiv_me = qdiv_th
qdiv_th = qdiv(0,:,:)-qdiv(1,:,:)
qdiv_me = qdiv(1,:,:)-qdiv(2,:,:)

qu_th = qu(0,:,:)
qu_me = qu_th
qu_th = qu(0,:,:)-qu(1,:,:)
qu_me = qu(1,:,:)-qu(2,:,:)

qv_th = qv(0,:,:)
qv_me = qv_th
qv_th = qv(0,:,:)-qv(1,:,:)
qv_me = qv(1,:,:)-qv(2,:,:)

res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 0.01
res@cnMaxLevelValF = 0.1
res@cnMinLevelValF = -0.1
res@cnFillPalette = "CBR_drywet"
vecres@vcRefMagnitudeF  = 100.             ; define vector ref mag
vecres@vcRefAnnoString1  = "100 Kg/m/s"


res@gsnLeftString = "thermal"
plot(0) = gsn_csm_contour_map(wks, qdiv_th, res)
plot_topo(0) = gsn_csm_contour(wks, tpmask, res_topo)
plotv(0) = gsn_csm_vector(wks, qu_th, qv_th, vecres)
overlay(plot(0), plot_topo(0))
overlay(plot(0), plotv(0))

res@gsnLeftString = "mechanical"
plot(1) = gsn_csm_contour_map(wks, qdiv_me, res)
plot_topo(1) = gsn_csm_contour(wks, tpmask, res_topo)
plotv(1) = gsn_csm_vector(wks, qu_me, qv_me, vecres)
overlay(plot(1), plot_topo(1))
overlay(plot(1), plotv(1))

gsn_panel(wks, plot(0:1), (/1,2/), pres)

end