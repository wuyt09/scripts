
begin

nlat = 96
nlon = 144

lonstr = 70
lonend = 120
latstr = 5
latend = 28

casename = (/"CTRL","NOQ1","NOTP"/)
; casename = "NOTP"
ncase = dimsizes(casename)
clim_pentad_z = new((/ncase,73,96,144/), float)
clim_pentad_t = clim_pentad_z
clim_pentad_u = clim_pentad_z
clim_pentad_v = clim_pentad_z
LEV = 850

do i = 0, ncase-1
    diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(i)+"_F2000_F19G16_CAM4/pro/daily/"
 
    uname = systemfunc("ls "+diri+"F2000_"+casename(i)+".plevel.U.*.daily.nc")
    vname = systemfunc("ls "+diri+"F2000_"+casename(i)+".plevel.V.*.daily.nc")
    zname = systemfunc("ls "+diri+"F2000_"+casename(i)+".plevel.Z3.*.daily.nc")
    tname = systemfunc("ls "+diri+"F2000_"+casename(i)+".plevel.T.*.daily.nc")
    print(zname)

    fz = addfiles(zname, "r")
    ft = addfiles(tname, "r")
    fu = addfiles(uname, "r")
    fv = addfiles(vname, "r")

    U = fu[:]->U(:,{LEV},:,:)
    u = reshape(U, (/25,365,nlat,nlon/))
    copy_VarCoords(U(0:364,:,:), u(0,:,:,:))

    V = fv[:]->V(:,{LEV},:,:)
    v = reshape(V, (/25,365,nlat,nlon/))
    copy_VarCoords(V(0:364,:,:), v(0,:,:,:))

    Z = fz[:]->Z3(:,{LEV},:,:)
    z = reshape(Z, (/25,365,nlat,nlon/))
    copy_VarCoords(Z(0:364,:,:), z(0,:,:,:))

    T = ft[:]->T(:,{LEV},:,:)
    t = reshape(T, (/25,365,nlat,nlon/))
    copy_VarCoords(T(0:364,:,:), t(0,:,:,:))

    pentad_u = u(:,0:72,:,:)
    pentad_v = v(:,0:72,:,:)
    pentad_z = z(:,0:72,:,:)
    pentad_t = t(:,0:72,:,:)

    do p = 0,72  ; a year has 73 penta
      pentad_u(:,p,:,:) = dim_sum_n_Wrap(u(:,5*p:5*p+4,:,:), 1)/5.0
      pentad_v(:,p,:,:) = dim_sum_n_Wrap(v(:,5*p:5*p+4,:,:), 1)/5.0
      pentad_z(:,p,:,:) = dim_sum_n_Wrap(z(:,5*p:5*p+4,:,:), 1)/5.0
      pentad_t(:,p,:,:) = dim_sum_n_Wrap(t(:,5*p:5*p+4,:,:), 1)/5.0
    end do

    clim_pentad_U = dim_avg_n_Wrap(pentad_u, 0)
    clim_pentad_u(i,:,:,:) = clim_pentad_U ;(18:35,:,:,:)
    delete([/U,u,pentad_u,clim_pentad_U/])

    clim_pentad_V = dim_avg_n_Wrap(pentad_v, 0)
    clim_pentad_v(i,:,:,:) = clim_pentad_V ;(18:35,:,:,:)
    delete([/V,v,pentad_v,clim_pentad_V/])

    clim_pentad_Z = dim_avg_n_Wrap(pentad_z, 0)
    clim_pentad_z(i,:,:,:) = clim_pentad_Z ;(18:35,:,:,:)
    delete([/Z,z,pentad_z,clim_pentad_Z/])

    clim_pentad_T = dim_avg_n_Wrap(pentad_t, 0)
    clim_pentad_t(i,:,:,:) = clim_pentad_T ;(18:35,:,:,:)
    delete([/T,t,pentad_t,clim_pentad_T/])

end do
printVarSummary(clim_pentad_z)

;-------------------------------------------------------------------

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20200513/"
figname = "H850_monsoon_trough_P27_P28"


wks = gsn_open_wks("pdf", figpath+figname)
nplot = 6
plot = new(nplot, graphic)
plotc = plot
plotv = plot

res = True
res@gsnDraw  = False
res@gsnFrame = False

res@mpOutlineOn = True
res@mpGeophysicalLineThicknessF = 2
res@mpGeophysicalLineColor = "grey25" ;"red" 
res@mpFillDrawOrder = "PreDraw"
res@mpFillOn = True

res@mpMinLatF = latstr
res@mpMaxLatF = latend
res@mpMinLonF = lonstr
res@mpMaxLonF = lonend

res@cnFillOn = True,
res@cnFillPalette = "BlueWhiteOrangeRed"
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False
res@lbLabelBarOn = False

res@gsnStringFontHeightF = 0.03
res@tmXBLabelFontHeightF = 0.02
res@tmYLLabelFontHeightF = 0.02
;=====================================================
    cnres = True
    cnres@gsnDraw = False
    cnres@gsnFrame = False

    cnres@cnFillOn = False
    cnres@cnLinesOn = True
    cnres@cnLineLabelsOn = False
    cnres@cnLineLabelFontColor = "red"
    cnres@cnLineThicknessF = 5
    cnres@cnInfoLabelOn = True
    cnres@cnInfoLabelFontHeightF = 0.02
    cnres@cnInfoLabelOrthogonalPosF = -0.1
    ;cnres@cnLineDrawOrder = "PostDraw"
    cnres@cnLabelDrawOrder = "PostDraw"
    ;cnres@cnLevelFlag = "LineAndLabel"
    cnres@cnLineLabelInterval = 1

;=====================================================

;------------------------------------------

  vcres                         = True             ; plot mods desired
  vcres@gsnFrame                = False            ; so we can draw time stamp
  vcres@gsnDraw                 = False
  vcres@vcRefAnnoOrthogonalPosF = -1.0  ;   -0.35          ; move ref vector up

  vcres@vcGlyphStyle            = "CurlyVector" ; "LineArrow"    ; turn on curly vectors
  ;vcres@vcMinDistanceF          = 0.017
  ;vcres@vcGlyphStyle            = "LineArrow"
  vcres@vcMinDistanceF          = 0.03
  vcres@vcRefLengthF            = 0.045            ; define length of vec ref
  
  vcres@vcLineArrowThicknessF   = 3.0
  ;vecres@vcLineArrowHeadMinSizeF = 0.008
  ;vecres@vcLineArrowHeadMaxSizeF = 0.3

  vcres@vcRefAnnoString2On = False
  vcres@vcRefAnnoFontHeightF = 0.02
;------------------------------------------


vcres@vcRefMagnitudeF         = 5.0  ; define vector ref mag

cnres@cnLevelSelectionMode = "ManualLevels"
cnres@cnMinLevelValF = 1460
cnres@cnMaxLevelValF = 1500
cnres@cnLevelSpacingF = 5

res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = 290
res@cnMaxLevelValF = 300
res@cnLevelSpacingF = 1

pens = 26
do i = 0, ncase-1
    res@gsnLeftString = casename(i)
    res@gsnRightString = "Pentad "+(pens+1)
    cnres@cnLineColor = "red"
    plot(i) = gsn_csm_contour_map(wks, clim_pentad_t(i,pens,:,:), res)
    plotc(i) = gsn_csm_contour(wks, clim_pentad_z(i,pens,:,:), cnres)
    plotv(i) = gsn_csm_vector(wks, clim_pentad_u(i,pens,:,:), clim_pentad_v(i,pens,:,:), vcres)
    overlay(plot(i), plotv(i))
    overlay(plot(i), plotc(i))


    res@gsnRightString = "Pentad "+(pens+2)
    plot(i+3) = gsn_csm_contour_map(wks, clim_pentad_t(i,pens+1,:,:), res)
    plotc(i+3) = gsn_csm_contour(wks, clim_pentad_z(i,pens+1,:,:), cnres)
    plotv(i+3) = gsn_csm_vector(wks, clim_pentad_u(i,pens+1,:,:), clim_pentad_v(i,pens+1,:,:), vcres)
    overlay(plot(i+3), plotv(i+3))
    overlay(plot(i+3), plotc(i+3))

end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = "T850 & Z850 "

gsn_panel(wks, plot, (/2,3/), pres)
end