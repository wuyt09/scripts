
begin

nlat = 96
nlon = 144
LEV = 200
season = "JJA"

lonstr = 0
lonend = 180
latstr = 0 ;-10
latend = 60

casename = (/"NOHEAT_CTRL","NOHEAT_ALL","NOHEAT_APR","NOHEAT_MAY"/)
ncase = dimsizes(casename)
year = ispan(6, 30, 1)
nyear = dimsizes(year)
season = "JJA"

STRT_ESM=1
END_ESM=30
ensemble = ispan(STRT_ESM, END_ESM, 1)
member = dimsizes(ensemble)

z = new((/ncase,6,member,nlat,nlon/), float)
z!0 = "case"
z&case = casename
z!1 = "time"
title = (/"April","May","June","July","August","JJA"/) ;"March",
z&time = title
z!2 = "member"
z&member = ensemble

u = z
v = z
t = z

do i = 0, ncase-1
  do j = 0, member-1
    diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(i)+"_F2000_F19G16_CAM4/pro/Ensemble_"+ensemble(j)+"/"
    tname = casename(i)+".T.plevel.mon.51.nc"
    zname = casename(i)+".Z3.plevel.mon.51.nc"
    uname = casename(i)+".U.plevel.mon.51.nc"
    vname = casename(i)+".V.plevel.mon.51.nc"
    fz = addfile(diri+zname, "r")
    ft = addfile(diri+tname, "r")
    fu = addfile(diri+uname, "r")
    fv = addfile(diri+vname, "r")
    Z = fz->Z3(:,{200},:,:)
    T = ft->T(:,{500},:,:)
    U = fu->U(:,{LEV},:,:)
    V = fv->V(:,{LEV},:,:)

    z(i,0:4,j,:,:) = Z(3:7,:,:)
    z_med = month_to_season(Z, season)
    z(i,5,j,:,:) = (/z_med(0,:,:)/)
    copy_VarCoords(Z(0,:,:), z(i,0,0,:,:))

    u(i,0:4,j,:,:) = U(3:7,:,:)
    u_med = month_to_season(U, season)
    u(i,5,j,:,:) = (/u_med(0,:,:)/)
    ; copy_VarCoords(U(0,:,:), u(i,0,0,:,:))

    v(i,0:4,j,:,:) = V(3:7,:,:)
    v_med = month_to_season(V, season)
    v(i,5,j,:,:) = (/v_med(0,:,:)/)
    ; copy_VarCoords(V(0,:,:), v(i,0,0,:,:))

    t(i,0:4,j,:,:) = T(3:7,:,:)
    t_med = month_to_season(T, season)
    t(i,5,j,:,:) = (/t_med(0,:,:)/)
  end do
end do
printVarSummary(u)

;=========================ttest================================================
u_var = dim_variance_n_Wrap(u, 2)
v_var = dim_variance_n_Wrap(v, 2)
z_var = dim_variance_n_Wrap(z, 2)
t_var = dim_variance_n_Wrap(t, 2)

u_ave = dim_avg_n_Wrap(u, 2)
v_ave = dim_avg_n_Wrap(v, 2)
z_ave = dim_avg_n_Wrap(z, 2)
t_ave = dim_avg_n_Wrap(t, 2)

diff_u = u_ave(0:ncase-2,:,:,:)
diff_v = v_ave(0:ncase-2,:,:,:)
diff_z = z_ave(0:ncase-2,:,:,:)
diff_t = t_ave(0:ncase-2,:,:,:)
prob_u = u_ave(0:ncase-2,:,:,:)
prob_v = v_ave(0:ncase-2,:,:,:)
prob_z = z_ave(0:ncase-2,:,:,:)
prob_t = t_ave(0:ncase-2,:,:,:)

do i = 0, ncase-2
    diff_u(i,:,:,:) = u_ave(0,:,:,:) - u_ave(i+1,:,:,:)
    diff_v(i,:,:,:) = v_ave(0,:,:,:) - v_ave(i+1,:,:,:)
    diff_z(i,:,:,:) = z_ave(0,:,:,:) - z_ave(i+1,:,:,:)
    diff_t(i,:,:,:) = t_ave(0,:,:,:) - t_ave(i+1,:,:,:)

    prob_u(i,:,:,:) = ttest(u_ave(0,:,:,:), u_var(0,:,:,:), nyear, u_ave(i+1,:,:,:), u_var(i+1,:,:,:), nyear, False, False)
    prob_v(i,:,:,:) = ttest(v_ave(0,:,:,:), v_var(0,:,:,:), nyear, v_ave(i+1,:,:,:), v_var(i+1,:,:,:), nyear, False, False)
    prob_z(i,:,:,:) = ttest(z_ave(0,:,:,:), z_var(0,:,:,:), nyear, z_ave(i+1,:,:,:), z_var(i+1,:,:,:), nyear, False, False)
    prob_t(i,:,:,:) = ttest(t_ave(0,:,:,:), t_var(0,:,:,:), nyear, t_ave(i+1,:,:,:), t_var(i+1,:,:,:), nyear, False, False)
end do

sig = 0.05
; diff_u = mask(diff_u, (prob_u.le.sig) .and. (prob_v.le.sig), True)
; diff_v = mask(diff_v, (prob_u.le.sig) .and. (prob_v.le.sig), True)
; diff_z = mask(diff_z, prob_z.le.sig, True)
; diff_t = mask(diff_t, prob_t.le.sig, True)
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
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_surf_temp_climate_effects/20210715_model/"
figname = "UV200_T500_SAH_Apr2Aug_casebycase_ensemble"

wks = gsn_open_wks("pdf", figpath+figname)

nplot = 6
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
 res_topo@cnLevels               = (/2000/) ;(/0.9/)
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

res@gsnStringFontHeightF = 0.035
res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025

;------------------------------------------

  vcres                         = True             ; plot mods desired
  vcres@gsnFrame                = False            ; so we can draw time stamp
  vcres@gsnDraw                 = False
  vcres@vcRefAnnoOrthogonalPosF = -1.0          ; move ref vector up

  vcres@vcGlyphStyle            = "LineArrow" ; "CurlyVector"    ; turn on curly vectors
  ;vcres@vcMinDistanceF          = 0.017
  ;vcres@vcGlyphStyle            = "LineArrow"
  vcres@vcMinDistanceF          = 0.025
  vcres@vcRefLengthF            = 0.05            ; define length of vec ref
  
  vcres@vcLineArrowThicknessF   = 3.0
  ;vecres@vcLineArrowHeadMinSizeF = 0.008
  ;vecres@vcLineArrowHeadMaxSizeF = 0.3

  vcres@vcRefAnnoString2On = False
  vcres@vcRefAnnoFontHeightF = 0.025
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

cnres@cnLevelSelectionMode = "ManualLevels"
cnres@cnLevelSpacingF = 50.
cnres@cnMaxLevelValF = 12650.
cnres@cnMinLevelValF = 12450.
;===========================================================

;res@lbLabelBarOn = False

;res@gsnRightString = "mm/day"

  res@pmLabelBarOrthogonalPosF = 0.2

  do i = 0, ncase-2
    res@gsnRightString = "CTRL - "+casename(i+1)
    vcres@vcRefMagnitudeF = 3.0
    res@cnFillPalette = "CBR_coldhot"
    res@cnLevelSelectionMode = "ManualLevels"
    res@cnLevelSpacingF   = 0.2
    res@cnMaxLevelValF    = 2.
    res@cnMinLevelValF    = -2.
    res@lbLabelBarOn = False
    do j = 0, 5
      res@gsnLeftString = title(j)
      plot(j) = gsn_csm_contour_map(wks, diff_t(i,j,:,:), res)
      plotv(j) = gsn_csm_vector(wks, diff_u(i,j,:,:), diff_v(i,j,:,:), vcres)
      plotc(j) = gsn_csm_contour(wks, z_ave(i+1,j,:,:), cnres)
      plot_topo(j) = gsn_csm_contour(wks, tpmask, res_topo)
      overlay(plot(j), plotv(j))
      ; overlay(plot(j), plotc(j))
      overlay(plot(j), plot_topo(j))
    end do
    pres = True
    pres@gsnMaximize = True
    pres@gsnPanelLabelBar = True
    pres@lbLabelFontHeightF = 0.012
    gsn_panel(wks, plot, (/2,3/), pres)
  end do
end