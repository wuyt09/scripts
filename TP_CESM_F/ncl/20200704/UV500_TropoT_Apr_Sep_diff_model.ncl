
begin

LEV = 200
nlat = 96
nlon = 144

; lonstr = 50
; lonend = 180
; latstr = -10
; latend = 40
lonstr = 50
lonend = 150
latstr = 0 ;-10
latend = 50

casename = (/"CTRL","NOTPQ1","NOTPQ1_MAR","NOTPQ1_APR","NOTPQ1_MAY"/)
ncase = dimsizes(casename)
year = ispan(6, 30, 1)
nyear = dimsizes(year)
month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
nmonth = dimsizes(month)

LEV = 500
z = new((/ncase,nyear,nmonth,96,144/), float)
u = z
v = z
levs = 500
levt = 200

do i = 0, ncase-1
  diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(i)+"_F2000_F19G16_CAM4/pro/"
  pname = "F2000_"+casename(i)+".PS.mon.0601-3012.nc"
  zname = "F2000_"+casename(i)+".plevel.T.mon.0601-3012.nc"
  uname = "F2000_"+casename(i)+".plevel.U.mon.0601-3012.nc"
  vname = "F2000_"+casename(i)+".plevel.V.mon.0601-3012.nc"
  fp = addfile(diri+pname, "r")
  fz = addfile(diri+zname, "r")
  fu = addfile(diri+uname, "r")
  fv = addfile(diri+vname, "r")
  PS = fp->PS
  T = fz->T ;(:,{LEV},:,:)
  plev = T&lev
  plev = plev*100
  plev@units = "Pa"
  dp = dpres_plevel_Wrap(plev, PS, min(plev), 0)
  Tdp = T*dp
  copy_VarCoords(T, Tdp)
  tdp = dim_sum_n_Wrap(Tdp(:,{levs:levt},:,:), 1)/dim_sum_n_Wrap(dp(:,{levs:levt},:,:), 1)
  copy_VarCoords(T(:,0,:,:), tdp)

  U = fu->U(:,{LEV},:,:)
  V = fv->V(:,{LEV},:,:)
  z(i,:,:,:,:) = reshape(tdp, (/nyear,nmonth,96,144/))
  copy_VarCoords(tdp(0,:,:), z(i,0,0,:,:))
  u(i,:,:,:,:) = reshape(U, (/nyear,nmonth,96,144/))
  v(i,:,:,:,:) = reshape(V, (/nyear,nmonth,96,144/))
end do
z!0 = "case"
z&case = casename
z!1 = "year"
z&year = year
z!2 = "month"
z&month = month
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

probu = u_ave(0:ncase-2,:,:,:)
probu!0 = "case"
case2 = "Ctrl -"+casename(1:ncase-1)
probu&case = case2
probv = probu
probz = probu
diff_u = probu
diff_v = probu
diff_z = probz

do i = 0, ncase-2
  probu(i,:,:,:) = ttest(u_ave(0,:,:,:), u_var(0,:,:,:), nyear, u_ave(i+1,:,:,:), u_var(i+1,:,:,:), nyear, False, False)
  probv(i,:,:,:) = ttest(v_ave(0,:,:,:), v_var(0,:,:,:), nyear, v_ave(i+1,:,:,:), v_var(i+1,:,:,:), nyear, False, False)
  probz(i,:,:,:) = ttest(z_ave(0,:,:,:), z_var(0,:,:,:), nyear, z_ave(i+1,:,:,:), z_var(i+1,:,:,:), nyear, False, False)
  diff_u(i,:,:,:) = u_ave(0,:,:,:) - u_ave(i+1,:,:,:)
  diff_v(i,:,:,:) = v_ave(0,:,:,:) - v_ave(i+1,:,:,:)
  diff_z(i,:,:,:) = z_ave(0,:,:,:) - z_ave(i+1,:,:,:)
end do

sig = 0.05
diff_u = mask(diff_u, (probu.lt.sig).or.(probv.lt.sig), True)
diff_v = mask(diff_v, (probu.lt.sig).or.(probv.lt.sig), True)
diff_z = mask(diff_z, probz.lt.sig, True)
printVarSummary(diff_z)

;----------------------------------------------------------------
tpmaskpath = "/home/yangsong3/wuyt/sysu/scripts/TP_CESM_F/model_need_input_data/"
; tpmaskname = "TPmask_1000m_C20190410.nc"
tpmaskname = "TPmask_1500m_C20200713.nc"
ftp = addfile(tpmaskpath+tpmaskname, "r")

tpmask = ftp->tp
printVarSummary(tpmask)
;-------------------------------------------------------------------
;-------------------------------------------------------------------
Month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20200704/"
figname = "UV500_TropoT_Apr_Jul_multi_model"

wks = gsn_open_wks("pdf", figpath+figname)

nplot = 4*(ncase-1)
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
 res_topo@cnLineColor            = "blue"
 res_topo@cnInfoLabelOn          = False
 res_topo@cnLevelSelectionMode   ="ExplicitLevels"
 res_topo@cnLevels               = (/0.5/)
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
    cnres@cnInfoLabelOn = True
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
;res@lbLabelBarOn = False

;res@gsnRightString = "mm/day"

  vcres@vcRefMagnitudeF         = 2.0  ; define vector ref mag
  
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnLevelSpacingF = 10
  res@cnMaxLevelValF = 100
  res@cnMinLevelValF = -100
  res@cnFillPalette = "CBR_coldhot"
  ; cnres@cnLevelSelectionMode = "ExplicitLevels"
  ; cnres@cnLevels = (/12460,12470,12480,12490,12500,15210,12520,12540/);12520,12560
  ; cnres@cnLevelSelectionMode = "ManualLevels"
  ; cnres@cnLevelSpacingF = 10.
  ; cnres@cnMaxLevelValF = 5900.
  ; cnres@cnMinLevelValF = 5880.

  ll = 2
  do j = 0,ncase-2
    res@gsnRightString = case2(j)
    do i = 0,3
      res@gsnLeftString = month(i+ll)
      plot(i*4+j) = gsn_csm_contour_map(wks, diff_z(j,i+ll,:,:), res)
      ; plotc(i*4+j) = gsn_csm_contour(wks, z_ave(0,i+ll,:,:), cnres)
      ; overlay(plot(i*4+j), plotc(i*4+j))
      plotv(i*4+j) = gsn_csm_vector(wks, diff_u(j,i+ll,:,:), diff_v(j,i+ll,:,:), vcres)
      overlay(plot(i*4+j), plotv(i*4+j))
      plot_topo(i*4+j) = gsn_csm_contour(wks, tpmask, res_topo)
      overlay(plot(i*4+j), plot_topo(i*4+j))
    end do
  end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
; pres@gsnPanelMainString = "UV200_H200"
;pres@pmLabelBarOrthogonalPosF = -0.01

gsn_panel(wks, plot, (/4,4/), pres)


end