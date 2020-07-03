
begin

casename = (/"CTRL","NOQ1","NOTP_NOQ1"/)
ncase = dimsizes(casename)
LEV = 850
U = new((/ncase,6,96,144/), float)
V = U
phis = U
precip = U

do i = 0, ncase-1
    diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(i)+"_F2000_F19G16_CAM4/pro/"
    filiu = "F2000_"+casename(i)+".plevel.U.mon.0601-3012.nc"
    filiv = "F2000_"+casename(i)+".plevel.V.mon.0601-3012.nc"
    filis = "F2000_"+casename(i)+".PHIS.mon.0601-3012.nc"
    filip1 = "F2000_"+casename(i)+".PRECL.mon.0601-3012.nc"
    filip2 = "F2000_"+casename(i)+".PRECC.mon.0601-3012.nc"

    fu = addfile(diri+filiu, "r")
    fv = addfile(diri+filiv, "r")
    fs = addfile(diri+filis, "r")
    fp1 = addfile(diri+filip1, "r")
    fp2 = addfile(diri+filip2, "r")

    UX = fu->U(:,{LEV},:,:)
    Ux = clmMonTLL(UX)
    VX = fv->V(:,{LEV},:,:)
    Vx = clmMonTLL(VX)
    U(i,:,:,:) = Ux(3:8,:,:)
    V(i,:,:,:) = Vx(3:8,:,:)
    PHISX = fs->PHIS
    PHIS = clmMonTLL(PHISX)
    phis(i,:,:,:) = PHIS(3:8,:,:)
    P1 = fp1->PRECL
    P2 = fp2->PRECC
    PRECPX = P2
    PRECPX = P1 + P2
    PRECPX = PRECPX*1000*86400
    PRECP = clmMonTLL(PRECPX)
    precip(i,:,:,:) = PRECP(3:8,:,:)
delete([/diri,filiu,filiv,filis,filip1,filip2,fu,fv,fs,fp1,fp2,UX,Ux,VX,Vx,PHISX,PHIS,P1,P2,PRECPX,PRECP/])
end do
phis = phis/9.81
printVarSummary(U)

U_zonal = dim_avg_n_Wrap(U, 2)
u = U - conform_dims(dimsizes(U), U_zonal, (/0,1,3/))
copy_VarCoords(U, u)
printVarSummary(u)

V_zonal = dim_avg_n_Wrap(V, 2)
v = V - conform_dims(dimsizes(V), V_zonal, (/0,1,3/))
copy_VarCoords(V, v)
printVarSummary(v)

;==============================================================
; plot
;==============================================================

lonstr = 20 ;20
lonend = 160 ;220
latstr = -20
latend = 70

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/AMJJAS/"
figname = "AMJJAS_850UV_deviations"
wks = gsn_open_wks("pdf", figpath+figname)
plot = new(6, graphic)
plot_topo = plot
plotc = plot

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
 res_topo@cnLevels               = (/1500/)
 res_topo@gsnLeftString = ""
 res_topo@gsnRightString = ""

res = True
res@gsnDraw  = False
res@gsnFrame = False

res@mpOutlineOn = True
res@mpFillOn = True
res@mpGeophysicalLineThicknessF = 2
res@mpGeophysicalLineColor = "grey45" ;"grey25"
res@mpFillDrawOrder = "PreDraw"
res@mpMinLonF = lonstr
res@mpMaxLonF = lonend
res@mpMinLatF = latstr
res@mpMaxLatF = latend

res@gsnStringFontHeightF = 0.035
res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025

res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False

res@lbLabelBarOn = False
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 2.
res@cnMaxLevelValF = 20.
res@cnMinLevelValF = 0.
res@cnFillPalette = "CBR_wet"

; res@stArrowLengthF = 0.01
; ; res@stArrowStride = 1
; res@stStepSizeF = 0.01
; ; res@stMinArrowSpacingF = 0.035
; res@stLineThicknessF = 3

;------------------------------------------

  vcres                         = True             ; plot mods desired
  vcres@gsnFrame                = False            ; so we can draw time stamp
  vcres@gsnDraw                 = False
  vcres@vcRefAnnoOrthogonalPosF = -1.0 ;   -0.35           ; move ref vector up
  ;vcres@vcGlyphStyle            = "CurlyVector"  ;"LineArrow"    ; turn on curly vectors
  ;vcres@vcMinDistanceF          = 0.017
  vcres@vcGlyphStyle            = "LineArrow"
  vcres@vcMinDistanceF          = 0.04
 
  vcres@vcLineArrowThicknessF   = 3.0
  ;vecres@vcLineArrowHeadMinSizeF = 0.008
  ;vecres@vcLineArrowHeadMaxSizeF = 0.3

  vcres@vcRefAnnoString2On = False
  vcres@vcRefAnnoFontHeightF = 0.02
  vcres@vcRefMagnitudeF         = 3.0             ; define vector ref mag
  vcres@vcRefLengthF            = 0.025            ; define length of vec ref
;------------------------------------------
pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True

do i = 0, ncase-1
    do j = 0, 5
        res@gsnLeftString = "Month = "+(j+4)
        ;plot(j) = gsn_csm_streamline_map_ce(wks, u(i,j,:,:), v(i,j,:,:), res)
        plotc(j) = gsn_csm_contour_map(wks, precip(i,j,:,:), res)
        plot(j) = gsn_csm_vector(wks, u(i,j,:,:), v(i,j,:,:), vcres)
        plot_topo(j) = gsn_csm_contour(wks, phis(i,j,:,:), res_topo)
        overlay(plotc(j), plot(j))
        overlay(plotc(j), plot_topo(j))
    end do
    pres@gsnPanelMainString = "Precip & 850hPa zonal deviation wind: "+casename(i)
    gsn_panel(wks, plotc, (/3,2/), pres)
end do



end