
begin

nlat = 96
nlon = 144
LEV = 850
season = "JJA"

lonstr = 30
lonend = 150
latstr = 0
latend = 60

casename = (/"NOHEAT_APR"/) ;
ncase = dimsizes(casename)
season = "JJA"

STRT_ESM=1
END_ESM=5
ensemble = ispan(STRT_ESM, END_ESM, 1)
member = dimsizes(ensemble)

pens = 18
pene = 26
npen = pene - pens + 1

P = new((/member,365,nlat,nlon/), float)
U = P
V = P
p = new((/ncase,npen,nlat,nlon/), float)
u = p
v = p

do i = 0, ncase-1
  do j = 0, member-1
    diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(i)+"_F2000_F19G16_CAM4/pro/Ensemble_"+ensemble(j)+"/"
    pname = casename(i)+".PREC.daily.51.nc"
    uname = casename(i)+".U.daily.51.nc"
    vname = casename(i)+".V.daily.51.nc"
    fp = addfile(diri+pname, "r")
    fu = addfile(diri+uname, "r")
    fv = addfile(diri+vname, "r")
    P_med = fp->PREC
    P_med = P_med*1000*86400 ;units mm/day
    P(j,:,:,:) = P_med
    copy_VarCoords(P_med, P(j,:,:,:))
    U_med = fu->U(:,{LEV},:,:)
    V_med = fv->V(:,{LEV},:,:)
    U(j,:,:,:) = U_med
    V(j,:,:,:) = V_med
    copy_VarCoords(U_med, U(j,:,:,:))
    copy_VarCoords(V_med, V(j,:,:,:))
  end do
  do ll = pens-1, pene-1
    l = ll - (pens-1)
    p_med = dim_sum_n_Wrap(P(:,5*ll:5*ll+4,:,:), 1)/5.0
    p(i,l,:,:)  = dim_avg_n_Wrap(p_med, 0)
    u_med = dim_sum_n_Wrap(U(:,5*ll:5*ll+4,:,:), 1)/5.0
    u(i,l,:,:)  = dim_avg_n_Wrap(u_med, 0)
    v_med = dim_sum_n_Wrap(V(:,5*ll:5*ll+4,:,:), 1)/5.0
    v(i,l,:,:)  = dim_avg_n_Wrap(v_med, 0)
  end do
  copy_VarCoords(P(0,0,:,:), p(0,0,:,:))
  copy_VarCoords(U(0,0,:,:), u(0,0,:,:))
  copy_VarCoords(V(0,0,:,:), v(0,0,:,:))
end do

p!0 = "case"
p&case = casename
p!1 = "time"
title = ispan(pens, pene, 1)
p&time = title
; z!2 = "year"
; z&year = year
printVarSummary(p)
copy_VarCoords(p, u)
copy_VarCoords(p, v)
printVarSummary(u)

; ;=========================ttest================================================
; u_var = dim_variance_n_Wrap(u, 2)
; v_var = dim_variance_n_Wrap(v, 2)
; z_var = dim_variance_n_Wrap(z, 2)
; p_var = dim_variance_n_Wrap(p, 2)

; u_ave = dim_avg_n_Wrap(u, 2)
; v_ave = dim_avg_n_Wrap(v, 2)
; z_ave = dim_avg_n_Wrap(z, 2)
; p_ave = dim_avg_n_Wrap(p, 2)

; diff_u = u_ave(0:ncase-2,:,:,:)
; diff_v = v_ave(0:ncase-2,:,:,:)
; diff_z = z_ave(0:ncase-2,:,:,:)
; diff_p = p_ave(0:ncase-2,:,:,:)
; prob_u = u_ave(0:ncase-2,:,:,:)
; prob_v = v_ave(0:ncase-2,:,:,:)
; prob_z = z_ave(0:ncase-2,:,:,:)
; prob_p = p_ave(0:ncase-2,:,:,:)

; do i = 0, ncase-2
;     diff_u(i,:,:,:) = u_ave(0,:,:,:) - u_ave(i+1,:,:,:)
;     diff_v(i,:,:,:) = v_ave(0,:,:,:) - v_ave(i+1,:,:,:)
;     diff_z(i,:,:,:) = z_ave(0,:,:,:) - z_ave(i+1,:,:,:)
;     diff_p(i,:,:,:) = p_ave(0,:,:,:) - p_ave(i+1,:,:,:)

;     prob_u(i,:,:,:) = ttest(u_ave(0,:,:,:), u_var(0,:,:,:), nyear, u_ave(i+1,:,:,:), u_var(i+1,:,:,:), nyear, False, False)
;     prob_v(i,:,:,:) = ttest(v_ave(0,:,:,:), v_var(0,:,:,:), nyear, v_ave(i+1,:,:,:), v_var(i+1,:,:,:), nyear, False, False)
;     prob_z(i,:,:,:) = ttest(z_ave(0,:,:,:), z_var(0,:,:,:), nyear, z_ave(i+1,:,:,:), z_var(i+1,:,:,:), nyear, False, False)
;     prob_p(i,:,:,:) = ttest(p_ave(0,:,:,:), p_var(0,:,:,:), nyear, p_ave(i+1,:,:,:), p_var(i+1,:,:,:), nyear, False, False)
; end do

; diff_u = u(0:ncase-2,:,:,:)
; diff_v = v(0:ncase-2,:,:,:)
; diff_z = z(0:ncase-2,:,:,:)
; diff_p = p(0:ncase-2,:,:,:)
; do i = 0, ncase-2
;   diff_u(i,:,:,:) = u(0,:,:,:) - u(i+1,:,:,:)
;   diff_v(i,:,:,:) = v(0,:,:,:) - v(i+1,:,:,:)
;   diff_z(i,:,:,:) = z(0,:,:,:) - z(i+1,:,:,:)
;   diff_p(i,:,:,:) = p(0,:,:,:) - p(i+1,:,:,:)
; end do

sig = 0.05
; diff_u = mask(diff_u, (prob_u.le.sig) .and. (prob_v.le.sig), True)
; diff_v = mask(diff_v, (prob_u.le.sig) .and. (prob_v.le.sig), True)
; diff_z = mask(diff_z, prob_z.le.sig, True)
; printVarSummary(diff_z)

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
figname = "UV850_Precip_Pentad_"+pens+"-"+pene+"_CTRL_thermal_ensemble"

wks = gsn_open_wks("pdf", figpath+figname)

nplot = npen
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

  ; vcres@vcRefMagnitudeF = 2.0
  ; res@cnFillPalette = "CBR_drywet"
  ; res@cnLevelSelectionMode = "ManualLevels"
  ; res@cnLevelSpacingF   = 0.5
  ; res@cnMaxLevelValF    = 5.
  ; res@cnMinLevelValF    = -5.
  vcres@vcRefMagnitudeF = 8.0
  res@cnFillPalette = "CBR_wet"
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnLevelSpacingF   = 2
  res@cnMaxLevelValF    = 20.
  res@cnMinLevelValF    = 1.
  res@lbLabelBarOn = False

  ; do i = 1, ncase-2
  i = 0
    ; res@gsnRightString = "CTRL - "+casename(i+1)
    do j = 0, npen-1
      res@gsnLeftString = "Pentad "+(pens+j)
      plot(j) = gsn_csm_contour_map(wks, p(i,j,:,:), res)
      plotv(j) = gsn_csm_vector(wks, u(i,j,:,:), v(i,j,:,:), vcres)
      plot_topo(j) = gsn_csm_contour(wks, tpmask, res_topo)
      overlay(plot(j), plotv(j))
      overlay(plot(j), plot_topo(j))
    end do
    pres = True
    pres@gsnMaximize = True
    pres@gsnPanelLabelBar = True

    gsn_panel(wks, plot, (/3,3/), pres)

  ; end do

end