
begin

nlat = 96
nlon = 144
LEV = 850
season = "JJA"

lonstr = 30
lonend = 150
latstr = 0
latend = 60

casename = (/"CTRL","NOQ1", "NOTPQ1_APR", "NOTPQ1_MAY"/) ;
ncase = dimsizes(casename)
year = ispan(6, 30, 1)
nyear = dimsizes(year)
season = "JJA"

z = new((/ncase,3,nyear,nlat,nlon/), float)
u = z
v = z
p = z

do i = 0, ncase-1
  diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(i)+"_F2000_F19G16_CAM4/pro/"
  pname = "F2000_"+casename(i)+".PREC.mon.0601-3012.nc"
  zname = "F2000_"+casename(i)+".plevel.Z3.mon.0601-3012.nc"
  uname = "F2000_"+casename(i)+".plevel.U.mon.0601-3012.nc"
  vname = "F2000_"+casename(i)+".plevel.V.mon.0601-3012.nc"
  fp = addfile(diri+pname, "r")
  fz = addfile(diri+zname, "r")
  fu = addfile(diri+uname, "r")
  fv = addfile(diri+vname, "r")
  P = fp->PREC
  P = P*1000*86400
  P@units = "mm/day"
  Z = fz->Z3(:,{500},:,:)
  U = fu->U(:,{LEV},:,:)
  V = fv->V(:,{LEV},:,:)

  z(i,0,:,:,:) = Z(3::12,:,:)
  z(i,1,:,:,:) = Z(4::12,:,:)
  z(i,2,:,:,:) = Z(5::12,:,:) ;month_to_season(Z, season)
  copy_VarCoords(Z(0,:,:), z(i,0,0,:,:))

  p(i,0,:,:,:) = P(3::12,:,:)
  p(i,1,:,:,:) = P(4::12,:,:)
  p(i,2,:,:,:) = P(5::12,:,:) ;month_to_season(P, season)
  copy_VarCoords(P(0,:,:), p(i,0,0,:,:))

  u(i,0,:,:,:) = U(3::12,:,:)
  u(i,1,:,:,:) = U(4::12,:,:)
  u(i,2,:,:,:) = U(5::12,:,:) ;month_to_season(U, season)
  ; copy_VarCoords(U(0,:,:), u(i,0,0,:,:))

  v(i,0,:,:,:) = V(3::12,:,:)
  v(i,1,:,:,:) = V(4::12,:,:)
  v(i,2,:,:,:) = V(5::12,:,:) ;month_to_season(V, season)
  ; copy_VarCoords(V(0,:,:), v(i,0,0,:,:))

end do
z!0 = "case"
z&case = casename
z!1 = "time"
title = (/"April","May","June"/)
z&time = title
z!2 = "year"
z&year = year
printVarSummary(z)
copy_VarCoords(z, u)
copy_VarCoords(z, v)
copy_VarCoords(z, p)
printVarSummary(u)

;=========================ttest================================================
u_var = dim_variance_n_Wrap(u, 2)
v_var = dim_variance_n_Wrap(v, 2)
z_var = dim_variance_n_Wrap(z, 2)
p_var = dim_variance_n_Wrap(p, 2)

u_ave = dim_avg_n_Wrap(u, 2)
v_ave = dim_avg_n_Wrap(v, 2)
z_ave = dim_avg_n_Wrap(z, 2)
p_ave = dim_avg_n_Wrap(p, 2)

diff_u = u_ave(0:ncase-2,:,:,:)
diff_v = v_ave(0:ncase-2,:,:,:)
diff_z = z_ave(0:ncase-2,:,:,:)
diff_p = p_ave(0:ncase-2,:,:,:)
prob_u = u_ave(0:ncase-2,:,:,:)
prob_v = v_ave(0:ncase-2,:,:,:)
prob_z = z_ave(0:ncase-2,:,:,:)
prob_p = p_ave(0:ncase-2,:,:,:)

do i = 0, ncase-2
    diff_u(i,:,:,:) = u_ave(0,:,:,:) - u_ave(i+1,:,:,:)
    diff_v(i,:,:,:) = v_ave(0,:,:,:) - v_ave(i+1,:,:,:)
    diff_z(i,:,:,:) = z_ave(0,:,:,:) - z_ave(i+1,:,:,:)
    diff_p(i,:,:,:) = p_ave(0,:,:,:) - p_ave(i+1,:,:,:)

    prob_u(i,:,:,:) = ttest(u_ave(0,:,:,:), u_var(0,:,:,:), nyear, u_ave(i+1,:,:,:), u_var(i+1,:,:,:), nyear, False, False)
    prob_v(i,:,:,:) = ttest(v_ave(0,:,:,:), v_var(0,:,:,:), nyear, v_ave(i+1,:,:,:), v_var(i+1,:,:,:), nyear, False, False)
    prob_z(i,:,:,:) = ttest(z_ave(0,:,:,:), z_var(0,:,:,:), nyear, z_ave(i+1,:,:,:), z_var(i+1,:,:,:), nyear, False, False)
    prob_p(i,:,:,:) = ttest(p_ave(0,:,:,:), p_var(0,:,:,:), nyear, p_ave(i+1,:,:,:), p_var(i+1,:,:,:), nyear, False, False)
end do

sig = 0.05
; diff_u = mask(diff_u, (prob_u.le.sig) .and. (prob_v.le.sig), True)
; diff_v = mask(diff_v, (prob_u.le.sig) .and. (prob_v.le.sig), True)
; diff_z = mask(diff_z, prob_z.le.sig, True)
printVarSummary(diff_z)

;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo ;(::-1,:)
printVarSummary(topo)
tpmask = topo

 ; ff = "~/wuyt/data/TOPO_1500.nc"
 ; topo_in = addfile(ff,"r")
 ; topo = topo_in->TP
 ; xlat = topo_in->LATIXY(:,0)
 ; xlon = topo_in->LONGXY(0,:)
 ; topo!0 = "lsmlat"
 ; topo!1 = "lsmlon"
 ; topo&lsmlat = xlat
 ; topo&lsmlon = xlon
 ; printVarSummary(topo)
 ; tpmask = topo
;-------------------------------------------------------------------
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_surf_temp_climate_effects/20210702/"
figname = "UV850_Precip_Apr_May_JJA_CTRL_thermal"

wks = gsn_open_wks("pdf", figpath+figname)

nplot = 12
plot = new(nplot, graphic)
plotv = plot
plotc = plot
plot_topo = plot

;----------------------------------------------------
 res_topo            = True
 res_topo@cnFillOn   = False
 res_topo@cnLinesOn  = True
 res_topo@gsnFrame   = False
 res_topo@gsnDraw    = False
 res_topo@cnLineLabelsOn         = False
 res_topo@cnLineThicknessF       = 4
 res_topo@cnLineDashPattern      = 0
 res_topo@cnLineColor            = "green4"
 res_topo@cnInfoLabelOn          = False
 res_topo@cnLevelSelectionMode   ="ExplicitLevels"
 res_topo@cnLevels               = 2000 ;(/0.9/)
 res_topo@gsnLeftString = ""
 res_topo@gsnRightString = ""
;----------------------------------------------------

res = True
res@gsnDraw  = False
res@gsnFrame = False

res@mpOutlineOn = True
res@mpGeophysicalLineThicknessF = 2
res@mpGeophysicalLineColor = "grey25" ;"red" 
res@mpFillDrawOrder = "PreDraw"
; res@mpFillOn = True

res@mpMinLatF = latstr
res@mpMaxLatF = latend
res@mpMinLonF = lonstr
res@mpMaxLonF = lonend

res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
; res@lbLabelBarOn = False
res@cnInfoLabelOn = False
;res@cnLineDrawOrder = "PostDraw"
; res@cnLabelDrawOrder = "PostDraw"
;res@cnLevelFlag = "LineAndLabel"
;res@cnLineLabelInterval = 1

res@gsnStringFontHeightF = 0.03
res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025

;------------------------------------------

  vcres                         = True             ; plot mods desired
  vcres@gsnFrame                = False            ; so we can draw time stamp
  vcres@gsnDraw                 = False
  vcres@vcRefAnnoOrthogonalPosF = -0.8          ; move ref vector up

  vcres@vcGlyphStyle            = "LineArrow" ; "CurlyVector"    ; turn on curly vectors
  ;vcres@vcMinDistanceF          = 0.017
  ;vcres@vcGlyphStyle            = "LineArrow"
  vcres@vcMinDistanceF          = 0.025
  vcres@vcRefLengthF            = 0.05            ; define length of vec ref
  
  vcres@vcLineArrowThicknessF   = 3.0
  ;vecres@vcLineArrowHeadMinSizeF = 0.008
  ;vecres@vcLineArrowHeadMaxSizeF = 0.3

  vcres@vcRefAnnoString2On = False
  vcres@vcRefAnnoFontHeightF = 0.03
;------------------------------------------
;===========================================================
cnres = True
cnres@gsnDraw  = False
cnres@gsnFrame = False

cnres@cnFillOn = False
cnres@cnLinesOn = True
cnres@cnLineLabelsOn = True
cnres@cnLineLabelFontColor = "red"
cnres@cnLineColor = "purple"
cnres@cnLineThicknessF = 5
cnres@cnInfoLabelOn = False
;cnres@cnLineDrawOrder = "PostDraw"
cnres@cnLabelDrawOrder = "PostDraw"
;cnres@cnLevelFlag = "LineAndLabel"
cnres@cnLineLabelInterval = 1

; cnres@cnLevelSelectionMode = "ManualLevels"
; cnres@cnLevelSpacingF = 30.
; cnres@cnMaxLevelValF = 12650.
; cnres@cnMinLevelValF = 12450.
cnres@cnLevelSelectionMode = "ExplicitLevels"
;===========================================================

;res@lbLabelBarOn = False

;res@gsnRightString = "mm/day"

  res@pmLabelBarOrthogonalPosF = 0.12

  
  i=0
  res@gsnRightString = "CTRL"
  res@cnFillPalette = "CBR_wet"
  res@cnLevelSelectionMode = "ManualLevels"
  vcres@vcRefMagnitudeF = 8.0
  res@cnLevelSpacingF   = 2.
  res@cnMaxLevelValF    = 20.
  res@cnMinLevelValF    = 2.
  cnres@cnLevels = (/5885/)
  do j = 0, 2
    res@gsnLeftString = title(j)
    plot(i+4*j) = gsn_csm_contour_map(wks, p_ave(0,j,:,:), res)
    plotv(i+4*j) = gsn_csm_vector(wks, u_ave(0,j,:,:), v_ave(0,j,:,:), vcres)
    ; plotc(i+4*j) = gsn_csm_contour(wks, z_ave(0,j,:,:), cnres)
    plot_topo(i+4*j) = gsn_csm_contour(wks, tpmask, res_topo)
    overlay(plot(i+4*j), plotv(i+4*j))
    ; overlay(plot(0), plotc(0))
    overlay(plot(i+4*j), plot_topo(i+4*j))
  end do

  i=1
  res@gsnRightString = "total thermal"
  vcres@vcRefMagnitudeF = 5.0
  res@cnFillPalette = "CBR_drywet"
  res@cnLevelSpacingF   = 0.5
  res@cnMaxLevelValF    = 5.
  res@cnMinLevelValF    = -5.
  do j = 0, 2
    res@gsnLeftString = title(j)
    plot(i+4*j) = gsn_csm_contour_map(wks, diff_p(i-1,j,:,:), res)
    plotv(i+4*j) = gsn_csm_vector(wks, diff_u(i-1,j,:,:), diff_v(i-1,j,:,:), vcres)
    ; plotc(i+4*j) = gsn_csm_contour(wks, diff_z(i-1,j,:,:), cnres)
    plot_topo(i+4*j) = gsn_csm_contour(wks, tpmask, res_topo)
    overlay(plot(i+4*j), plotv(i+4*j))
    ; overlay(plot(0), plotc(0))
    overlay(plot(i+4*j), plot_topo(i+4*j))
  end do

  i=2
  res@gsnRightString = "Apr thermal"
  vcres@vcRefMagnitudeF = 2.0
  ; res@cnFillPalette = "CBR_drywet"
  ; res@cnLevelSpacingF   = 1.
  ; res@cnMaxLevelValF    = 10.
  ; res@cnMinLevelValF    = -10.
  do j = 0, 2
    res@gsnLeftString = title(j)
    plot(i+4*j) = gsn_csm_contour_map(wks, diff_p(i-1,j,:,:), res)
    plotv(i+4*j) = gsn_csm_vector(wks, diff_u(i-1,j,:,:), diff_v(i-1,j,:,:), vcres)
    ; plotc(i+4*j) = gsn_csm_contour(wks, diff_z(i-1,j,:,:), cnres)
    plot_topo(i+4*j) = gsn_csm_contour(wks, tpmask, res_topo)
    overlay(plot(i+4*j), plotv(i+4*j))
    ; overlay(plot(0), plotc(0))
    overlay(plot(i+4*j), plot_topo(i+4*j))
  end do

  i=3
  res@gsnRightString = "May thermal"
  ; vcres@vcRefMagnitudeF = 5.0
  ; res@cnFillPalette = "CBR_drywet"
  ; res@cnLevelSpacingF   = 1.
  ; res@cnMaxLevelValF    = 10.
  ; res@cnMinLevelValF    = -10.
  do j = 0, 2
    res@gsnLeftString = title(j)
    plot(i+4*j) = gsn_csm_contour_map(wks, diff_p(i-1,j,:,:), res)
    plotv(i+4*j) = gsn_csm_vector(wks, diff_u(i-1,j,:,:), diff_v(i-1,j,:,:), vcres)
    ; plotc(i+4*j) = gsn_csm_contour(wks, diff_z(i-1,j,:,:), cnres)
    plot_topo(i+4*j) = gsn_csm_contour(wks, tpmask, res_topo)
    overlay(plot(i+4*j), plotv(i+4*j))
    ; overlay(plot(0), plotc(0))
    overlay(plot(i+4*j), plot_topo(i+4*j))
  end do

  pres = True
  pres@gsnMaximize = True
  ; pres@gsnPanelMainString = "UV850 & Precip"
  ; pres@gsnPanelMainFontHeightF = 0.02
  gsn_panel(wks, plot, (/3,4/), pres)
  gsn_panel(wks, (/plot(1),plot(5),plot(9)/), (/3,1/), pres)
  gsn_panel(wks, (/plot(2),plot(6),plot(10)/), (/3,1/), pres)
  gsn_panel(wks, (/plot(3),plot(7),plot(11)/), (/3,1/), pres)
end