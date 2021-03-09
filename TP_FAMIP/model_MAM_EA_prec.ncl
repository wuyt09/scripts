
begin

LEV = 200
nlat = 96
nlon = 144

; lonstr = 108
; lonend = 120
; latstr = 20
; latend = 31
lonstr = 60
lonend = 122
latstr = 10
latend = 40

casename = (/"Ctrl_FAMIP","Hist_FAMIP","Hist_FAMIP_CMIP6","Hist_PO"/) ;
ncase = dimsizes(casename)
year = ispan(1979, 2005, 1)
nyear = dimsizes(year)
month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
nmonth = dimsizes(month)
season = "MAM"

LEV = 850
z = new((/ncase,nyear,96,144/), float)
u = z
v = z

do i = 0, ncase-1
  diri = "/home/yangsong3/data-model/wuyt/TPbias/TPbias_"+casename(i)+"/a20191206/"
  zname = casename(i)+".cam.h0.PREC.1979-2005.nc"
  uname = casename(i)+".cam.h0.U.1979-2005.nc"
  vname = casename(i)+".cam.h0.V.1979-2005.nc"
  fz = addfile(diri+zname, "r")
  fu = addfile(diri+uname, "r")
  fv = addfile(diri+vname, "r")
  Z = fz->PREC
  Z = Z*1000*86400
  Z@units = "mm/day"
  U = fu->U(:,{LEV},:,:)
  V = fv->V(:,{LEV},:,:)
  z(i,:,:,:) = month_to_season(Z, season)
  ; copy_VarCoords(Z(0,:,:), z(i,0,:,:))
  u(i,:,:,:) = month_to_season(U, season)
  v(i,:,:,:) = month_to_season(V, season)
end do
z!0 = "case"
z&case = casename
; z!1 = "year"
; z&year = year
; z!2 = "month"
; z&month = month
printVarSummary(z)
copy_VarCoords(z, u)
copy_VarCoords(z, v)

;------------ttest--------------------------------------
u_var = dim_variance_n_Wrap(u, 1)
v_var = dim_variance_n_Wrap(v, 1)
z_var = dim_variance_n_Wrap(z, 1)

u_ave = dim_avg_n_Wrap(u, 1)
v_ave = dim_avg_n_Wrap(v, 1)
z_ave = dim_avg_n_Wrap(z, 1)



diff_z = z_ave(1:ncase-1,:,:)
diff_u = u_ave(1:ncase-1,:,:)
diff_v = v_ave(1:ncase-1,:,:)
probz = diff_z
probu = diff_u
probv = diff_v


sig = 0.05
do i = 0, ncase-2
  probu(i,:,:) = ttest(u_ave(0,:,:), u_var(0,:,:), nyear, u_ave(i+1,:,:), u_var(i+1,:,:), nyear, False, False)
  probv(i,:,:) = ttest(v_ave(0,:,:), v_var(0,:,:), nyear, v_ave(i+1,:,:), v_var(i+1,:,:), nyear, False, False)
  probz(i,:,:) = ttest(z_ave(0,:,:), z_var(0,:,:), nyear, z_ave(i+1,:,:), z_var(i+1,:,:), nyear, False, False)

  diff_z(i,:,:) = z_ave(i+1,:,:) - z_ave(0,:,:)
  diff_u(i,:,:) = u_ave(i+1,:,:) - u_ave(0,:,:)
  diff_v(i,:,:) = v_ave(i+1,:,:) - v_ave(0,:,:)
end do

; diff_u = mask(diff_u, (probu.lt.sig).or.(probv.lt.sig), True)
; diff_v = mask(diff_v, (probu.lt.sig).or.(probv.lt.sig), True)
; diff_z = mask(diff_z, probz.lt.sig, True)
; printVarSummary(diff_z)
;----------------------------------------------------------------
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)

printVarSummary(topo)
;-------------------------------------------------------------------
;-------------------------------------------------------------------


figpath = "/home/yangsong3/wuyt/sysu/figures/TP_FAMIP/"
figname = "model_"+season+"_EA_prec"

wks = gsn_open_wks("pdf", figpath+figname)

nplot = ncase-1
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
 res_topo@cnLineDashPattern      = 16
 res_topo@cnLineColor            = "gray"
 res_topo@cnInfoLabelOn          = False
 res_topo@cnLevelSelectionMode   ="ExplicitLevels"
 res_topo@cnLevels               = (/1500/)
 res_topo@gsnLeftString = ""
 res_topo@gsnRightString = ""
;----------------------------------------------------

res = True
res@gsnDraw  = False
res@gsnFrame = False
res@gsnMaximize = True

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
res@lbLabelBarOn = False
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
  vcres@vcRefAnnoOrthogonalPosF = -0.25          ; move ref vector up

  vcres@vcGlyphStyle            = "LineArrow" ;"CurlyVector"     ; turn on curly vectors
  ;vcres@vcMinDistanceF          = 0.017
  ;vcres@vcGlyphStyle            = "LineArrow"
  vcres@vcMinDistanceF          = 0.03
  vcres@vcRefLengthF            = 0.045            ; define length of vec ref
  
  vcres@vcLineArrowThicknessF   = 3.0
  ;vecres@vcLineArrowHeadMinSizeF = 0.008
  ;vecres@vcLineArrowHeadMaxSizeF = 0.3

  vcres@vcRefAnnoString2On = False
  vcres@vcRefAnnoFontHeightF = 0.02
  vcres@gsnLeftString = ""
  vcres@gsnRightString = ""
;------------------------------------------

  vcres@vcRefMagnitudeF         = 3.0  ; define vector ref mag
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnLevelSpacingF = 0.2
  res@cnMaxLevelValF = 2.
  res@cnMinLevelValF = -2.
  res@cnFillPalette = "CBR_drywet"

  do i = 0, ncase-2
    res@gsnLeftString = casename(i+1)+" - Ctrl"
    plot(i) = gsn_csm_contour_map(wks, diff_z(i,:,:), res)
    plotv(i) = gsn_csm_vector(wks, diff_u(i,:,:), diff_v(i,:,:), vcres)
    overlay(plot(i), plotv(i))
    plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(i), plot_topo(i))
  end do


pres = True
pres@gsnMaximize = True
pres@gsnPanelMainString = season+"_UV850_Precip"
pres@gsnPanelLabelBar = True
pres@lbLabelFontHeightF = 0.015
pres@pmLabelBarOrthogonalPosF = -0.01

gsn_panel(wks, plot, (/1,ncase-1/), pres)

end