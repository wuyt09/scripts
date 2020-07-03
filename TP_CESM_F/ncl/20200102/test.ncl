
begin

lonW = 75
lonE = 105
latS = 25
latN = 40

casename = (/"CTRL","NOQ1","NOTP_NOQ1"/) ;,"NOTP"
ncase = dimsizes(casename)
; month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
; nm = dimsizes(month)
season = (/"JJA"/)

LEV = (/925/)
u = new((/ncase,37,96,144/), float)
v = u
w = u
t = u

do i = 0, dimsizes(casename)-1
    diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(i)+"_F2000_F19G16_CAM4/pro/"

    filiu = "F2000_"+casename(i)+".plevel.U.mon.0601-3012.nc"
    fu = addfile(diri+filiu, "r")
    U1 = fu->U
    u1 = month_to_season(U1, season)
    u(i,:,:,:) = dim_avg_n_Wrap(u1, 0)

    filiv = "F2000_"+casename(i)+".plevel.V.mon.0601-3012.nc"
    fv = addfile(diri+filiv, "r")
    V1 = fv->V
    v1 = month_to_season(V1, season)
    v(i,:,:,:) = dim_avg_n_Wrap(v1, 0)

    filiw = "F2000_"+casename(i)+".plevel.OMEGA.mon.0601-3012.nc"
    fw = addfile(diri+filiw, "r")
    W1 = fw->OMEGA
    w1 = month_to_season(W1, season)
    w(i,:,:,:) = dim_avg_n_Wrap(w1, 0)

    filit = "F2000_"+casename(i)+".plevel.T.mon.0601-3012.nc"
    ft = addfile(diri+filit, "r")
    T1 = ft->T
    t1 = month_to_season(T1, season)
    t(i,:,:,:) = dim_avg_n_Wrap(t1, 0)
end do
w = w*(-100)
w@units = "10^-2*s^-1"

lat = u&lat
lon = u&lon

ux = u(0:1,:,:,:)
ux(0,:,:,:) = u(0,:,:,:) - u(1,:,:,:)
ux(1,:,:,:) = u(1,:,:,:) - u(2,:,:,:)

vx = v(0:1,:,:,:)
vx(0,:,:,:) = v(0,:,:,:) - v(1,:,:,:)
vx(1,:,:,:) = v(1,:,:,:) - v(2,:,:,:)

wx = w(0:1,:,:,:)
wx(0,:,:,:) = w(0,:,:,:) - w(1,:,:,:)
wx(1,:,:,:) = w(1,:,:,:) - w(2,:,:,:)

tx = t(0:1,:,:,:)
tx(0,:,:,:) = t(0,:,:,:) - t(1,:,:,:)
tx(1,:,:,:) = t(1,:,:,:) - t(2,:,:,:)


ux_zonal = dim_avg_n_Wrap(ux(:,:,{latS:latN},:), 2)
vx_zonal = dim_avg_n_Wrap(vx(:,:,{latS:latN},:), 2)
wx_zonal = dim_avg_n_Wrap(wx(:,:,{latS:latN},:), 2)
tx_zonal = dim_avg_n_Wrap(tx(:,:,{latS:latN},:), 2)

ux_merid = dim_avg_n_Wrap(ux(:,:,:,{lonW:lonE}), 3)
vx_merid = dim_avg_n_Wrap(vx(:,:,:,{lonW:lonE}), 3)
wx_merid = dim_avg_n_Wrap(wx(:,:,:,{lonW:lonE}), 3)
tx_merid = dim_avg_n_Wrap(tx(:,:,:,{lonW:lonE}), 3)
printVarSummary(tx_merid)

;----------------------------------------------------------------
ftopo = addfile("./USGS-gtopo30_1.9x2.5_remap_c050602.nc", "r")
topo = ftopo->PHIS
topo = topo/9.8
topo@units = "gpm"
printVarSummary(topo)
topo_z = dim_avg_n_Wrap(topo({latS:latN},:), 0)
topo_m = dim_avg_n_Wrap(topo(:,{lonW:lonE}), 1)
topo_zonal = topo_z
topo_merid = topo_m
topo_zonal = 1013.25*(1-topo_z*0.0065/288.15)^5.25145 
topo_merid = 1013.25*(1-topo_m*0.0065/288.15)^5.25145
copy_VarCoords(topo(0,:), topo_zonal)
printVarSummary(topo_zonal)
print(topo_zonal({90}))

; ux_zonal = mask(ux_zonal, topo_zonal.le.1000, True)
; vx_zonal = mask(vx_zonal, topo_zonal.le.1000, True)
; wx_zonal = mask(wx_zonal, topo_zonal.le.1000, True)
; tx_zonal = mask(tx_zonal, topo_zonal.le.1000, True)

; ux_zonal = where(ismissing(ux_zonal), 0, ux_zonal)
; vx_zonal = where(ismissing(vx_zonal), 0, vx_zonal)
; wx_zonal = where(ismissing(wx_zonal), 0, wx_zonal)
; tx_zonal = where(ismissing(tx_zonal), 0, tx_zonal)


;----------------------------------------------------------------

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20200102/"

figname = "JJA_Vertical_cell_T"

wks = gsn_open_wks("pdf", figpath+figname)
nplot = 4
plot = new(nplot, graphic)
plotc = plot
plot_topo = plot

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
    res_topo@trYReverse = True
;=====================================================

    res = True
    res@gsnDraw  = False
    res@gsnFrame = False
    res@vpHeightF = 0.4
    res@vpWidthF  = 0.7

    res@cnFillOn = True
    res@cnLinesOn = False
    res@cnLineLabelsOn = False
    res@cnInfoLabelOn = True

    res@gsnStringFontHeightF = 0.03
    res@tmXBLabelFontHeightF = 0.02
    res@tmYLLabelFontHeightF = 0.02

    ;res@pmLabelBarOrthogonalPosF = 0.12
    res@lbLabelBarOn = False
    ;res@lbLabelFontHeightF = 0.02
    res@cnFillPalette = "BlueWhiteOrangeRed"
    res@tiYAxisString = ""


 ;------------------------------------------
     res@trYMinF   = 100
     res@trYMaxF   = 1000
      
     res@tmYLMode      = "Explicit"
     res@tmYLValues    = (/"1000","850","750","500","300","200","100"/)
     res@tmYLLabels    = (/"1000","850","750","500","300","200","100"/)

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

;------------------------------------------

;=====================================================

    cnres = True
    cnres@gsnDraw = False
    cnres@gsnFrame = False

    cnres@cnFillOn = False
    cnres@cnLinesOn = True
    cnres@cnLineLabelsOn = True
    cnres@cnLineLabelFontColor = "red"
    cnres@cnLineColor = "green4"
    cnres@cnLineThicknessF = 5
    cnres@cnInfoLabelOn = False
    ;cnres@cnLineDrawOrder = "PostDraw"
    cnres@cnLabelDrawOrder = "PostDraw"
    ;cnres@cnLevelFlag = "LineAndLabel"
    cnres@cnLineLabelInterval = 1
    cnres@gsnPresHgtHeightLabelOn = False

;=====================================================

res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 2
res@cnMaxLevelValF = 10
res@cnMinLevelValF = -10
res@cnFillPalette = "BlueWhiteOrangeRed"
res@vcRefMagnitudeF         = 10.             ; define vector ref mag

cnres@cnLevelSelectionMode = "ExplicitLevels"
cnres@cnLevels = (/-10,-7,-4,-2,-0.5,0.5,2,4,7,10/)

cnres@gsnContourPosLineDashPattern = 0
cnres@gsnContourNegLineDashPattern = 2

res@gsnLeftString = "Thermal"
res@gsnRightString = latS+"~S~o~N~N-"+latN+"~S~o~N~N"
plot(0) = gsn_csm_pres_hgt_vector(wks, wx_zonal(0,{1000:100},{40:180}),\
 ux_zonal(0,{1000:100},{40:180}), wx_zonal(0,{1000:100},{40:180}), res)

plotc(0) = gsn_csm_pres_hgt(wks, tx_zonal(0,{1000:100},{40:180}), cnres)
overlay(plot(0), plotc(0))

plot_topo(0) = gsn_csm_xy(wks, lon, topo_zonal, res_topo)




pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar  = True
gsn_panel(wks, plot, (/2,2/), pres)


end

