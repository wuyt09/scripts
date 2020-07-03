begin

nlat = 96
nlon = 144

; lonstr = 50
; lonend = 180
; latstr = -10
; latend = 40
lonstr = 50
lonend = 150
latstr = 10
latend = 50

casename = (/"CTRL","NOQ1"/)
ncase = dimsizes(casename)
year = ispan(6, 30, 1)
nyear = dimsizes(year)
month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
nmonth = dimsizes(month)

LEV = 850
z = new((/ncase,nyear,nmonth,96,144/), float)
u = z
v = z
dv = z

do i = 0, ncase-1
  diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(i)+"_F2000_F19G16_CAM4/pro/"
  zname = "F2000_"+casename(i)+".plevel.Q.mon.0601-3012.nc"
  uname = "F2000_"+casename(i)+".plevel.U.mon.0601-3012.nc"
  vname = "F2000_"+casename(i)+".plevel.V.mon.0601-3012.nc"
  fz = addfile(diri+zname, "r")
  fu = addfile(diri+uname, "r")
  fv = addfile(diri+vname, "r")
  Z = fz->Q(:,{LEV},:,:)
  Z = Z*1000 ;units g/Kg
  U = fu->U(:,{LEV},:,:)
  V = fv->V(:,{LEV},:,:)
  QU = Z*U/9.8
  QV = Z*V/9.8

  lat = Z&lat
  lon = Z&lon
  DV = uv2dv_cfd(QU, QV, lat, lon, 1)
  DV = DV*1.e+6

  u(i,:,:,:,:) = reshape(QU, (/nyear,nmonth,96,144/))
  copy_VarCoords(Z(0,:,:), u(i,0,0,:,:))
  v(i,:,:,:,:) = reshape(QV, (/nyear,nmonth,96,144/))
  dv(i,:,:,:,:) = reshape(DV, (/nyear,nmonth,96,144/))
end do

u!0 = "case"
u&case = casename
u!1 = "year"
u&year = year
u!2 = "month"
u&month = month
printVarSummary(u)
copy_VarCoords(u, v)
copy_VarCoords(u, dv)

;------------ttest--------------------------------------
u_var = dim_variance_n_Wrap(u, 1)
v_var = dim_variance_n_Wrap(v, 1)
dv_var = dim_variance_n_Wrap(dv, 1)

u_ave = dim_avg_n_Wrap(u, 1)
v_ave = dim_avg_n_Wrap(v, 1)
dv_ave = dim_avg_n_Wrap(dv, 1)


probu = ttest(u_ave(0,:,:,:), u_var(0,:,:,:), nyear, u_ave(1,:,:,:), u_var(1,:,:,:), nyear, False, False)
probv = ttest(v_ave(0,:,:,:), v_var(0,:,:,:), nyear, v_ave(1,:,:,:), v_var(1,:,:,:), nyear, False, False)
probdv = ttest(dv_ave(0,:,:,:), dv_var(0,:,:,:), nyear, dv_ave(1,:,:,:), dv_var(1,:,:,:), nyear, False, False)

diff_u = u_ave(0,:,:,:)
diff_u = u_ave(0,:,:,:) - u_ave(1,:,:,:)
diff_v = v_ave(0,:,:,:)
diff_v = v_ave(0,:,:,:) - v_ave(1,:,:,:)
diff_dv = dv_ave(0,:,:,:)
diff_dv = dv_ave(0,:,:,:) - dv_ave(1,:,:,:)

sig = 0.05
diff_u = mask(diff_u, (probu.lt.sig).or.(probv.lt.sig), True)
diff_v = mask(diff_v, (probu.lt.sig).or.(probv.lt.sig), True)
diff_dv = mask(diff_dv, probdv.lt.sig, True)
printVarSummary(diff_u)

;----------------------------------------------------------------
tpmaskpath = "/home/yangsong3/wuyt/sysu/scripts/TP_CESM_F/model_need_input_data/"
tpmaskname = "TPmask_1000m_C20190410.nc"
ftp = addfile(tpmaskpath+tpmaskname, "r")

tpmask = ftp->tp
printVarSummary(tpmask)
;-------------------------------------------------------------------
;-------------------------------------------------------------------
Month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20200513/"
figname = "WV850_Apr_Jul_"+casename(0)+"_"+casename(1)

wks = gsn_open_wks("pdf", figpath+figname)

nplot = 6
plot = new(nplot, graphic)
plotv = plot
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
 res_topo@cnLineColor            = "blue"
 res_topo@cnInfoLabelOn          = False
 res_topo@cnLevelSelectionMode   ="ExplicitLevels"
 res_topo@cnLevels               = (/0.9/)
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
res@mpFillOn = True

res@mpMinLatF = latstr
res@mpMaxLatF = latend
res@mpMinLonF = lonstr
res@mpMaxLonF = lonend

res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@lbLabelBarOn = False
res@cnInfoLabelOn = False

res@gsnStringFontHeightF = 0.035
res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025

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
;res@lbLabelBarOn = False

;res@gsnRightString = "mm/day"

  vcres@vcRefMagnitudeF         = 5.0  ; define vector ref mag
  res@cnFillPalette = "CBR_drywet"
  res@gsnSpreadColorEnd = 0
  res@gsnSpreadColorStart = 10
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnLevelSpacingF = 0.5
  res@cnMaxLevelValF = 3.
  res@cnMinLevelValF = -3.


  do i = 0,5
  res@gsnLeftString = Month(i+3)
  
  plot(i) = gsn_csm_contour_map(wks, diff_dv(i+3,:,:), res)
  plotv(i) = gsn_csm_vector(wks, diff_u(i+3,:,:), diff_v(i+3,:,:), vcres)
  overlay(plot(i), plotv(i))
  plot_topo(i) = gsn_csm_contour(wks, tpmask, res_topo)
  overlay(plot(i), plot_topo(i))
end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = "WV850_divergence"
;pres@pmLabelBarOrthogonalPosF = -0.01

gsn_panel(wks, plot, (/2,3/), pres)

end