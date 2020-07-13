
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

casename = (/"CTRL","NOQ1"/) ;"NOTP"TP
ncase = dimsizes(casename)
year = ispan(6, 30, 1)
nyear = dimsizes(year)
;month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
month = (/"Apr","May","Jun","Jul","Aug","Sep"/)
nmonth = dimsizes(month)
LEV = 850
z = new((/ncase,nmonth,96,144/), float)
u = z
v = z

do i = 0, ncase-1
  diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(i)+"_F2000_F19G16_CAM4/pre/"
  fili = systemfunc("ls "+diri+"F2000_"+casename(i)+".cam.h0.0002-0[4-9].nc")
  print(fili)
  f = addfiles(fili, "r")
  C = f[:]->PRECC
  L = f[:]->PRECL
  Z = C 
  Z = (C+L)*1000*86400
  Z@units = "mm/day"
  U = f[:]->U(:,22,:,:)
  V = f[:]->V(:,22,:,:)
  z(i,:,:,:) = Z ;reshape(Z, (/nyear,nmonth,96,144/))
  copy_VarCoords(Z(0,:,:), z(i,0,:,:))
  u(i,:,:,:) = U ;reshape(U, (/nyear,nmonth,96,144/))
  v(i,:,:,:) = V ;reshape(V, (/nyear,nmonth,96,144/))
end do
printVarSummary(z)
; z!0 = "case"
; z&case = casename
; z!1 = "month"
; z&month = (/"Apr","May","Jun","Jul","Aug","Sep"/)
; printVarSummary(z)
; copy_VarCoords(z, u)
; copy_VarCoords(z, v)

;------------ttest--------------------------------------
; u_var = dim_variance_n_Wrap(u, 1)
; v_var = dim_variance_n_Wrap(v, 1)
; z_var = dim_variance_n_Wrap(z, 1)

; u_ave = dim_avg_n_Wrap(u, 1)
; v_ave = dim_avg_n_Wrap(v, 1)
; z_ave = dim_avg_n_Wrap(z, 1)


; probu = ttest(u_ave(0,:,:,:), u_var(0,:,:,:), nyear, u_ave(1,:,:,:), u_var(1,:,:,:), nyear, False, False)
; probv = ttest(v_ave(0,:,:,:), v_var(0,:,:,:), nyear, v_ave(1,:,:,:), v_var(1,:,:,:), nyear, False, False)
; probz = ttest(z_ave(0,:,:,:), z_var(0,:,:,:), nyear, z_ave(1,:,:,:), z_var(1,:,:,:), nyear, False, False)

diff_z = z(0,:,:,:)
diff_z = z(0,:,:,:) - z(1,:,:,:)
diff_u = u(0,:,:,:)
diff_u = u(0,:,:,:) - u(1,:,:,:)
diff_v = v(0,:,:,:)
diff_v = v(0,:,:,:) - v(1,:,:,:)

; sig = 0.05
; diff_u = mask(diff_u, (probu.lt.sig).or.(probv.lt.sig), True)
; diff_v = mask(diff_v, (probu.lt.sig).or.(probv.lt.sig), True)
; diff_z = mask(diff_z, probz.lt.sig, True)
printVarSummary(diff_z)
;----------------------------------------------------------------
tpmaskpath = "/home/yangsong3/wuyt/sysu/scripts/TP_CESM_F/model_need_input_data/"
tpmaskname = "TPmask_1000m_C20190410.nc"
ftp = addfile(tpmaskpath+tpmaskname, "r")

tpmask = ftp->tp
printVarSummary(tpmask)
;-------------------------------------------------------------------
;-------------------------------------------------------------------

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20200704/"
figname = "UV850_Precip_Apr_Jul_noQ1"

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
 res_topo@cnLineDashPattern      = 16
 res_topo@cnLineColor            = "gray"
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

;=====================================================
    cnres = True
    cnres@gsnDraw = False
    cnres@gsnFrame = False

    cnres@cnFillOn = False
    cnres@cnLinesOn = True
    cnres@cnLineLabelsOn = True
    cnres@cnLineLabelFontColor = "red"
    cnres@cnLineColor = "purple" ;goldenrod2 ;dodgerblue2
    cnres@cnLineThicknessF = 5
    cnres@cnInfoLabelOn = True
    ;cnres@cnLineDrawOrder = "PostDraw"
    cnres@cnLabelDrawOrder = "PostDraw"
    ;cnres@cnLevelFlag = "LineAndLabel"
    cnres@cnLineLabelInterval = 1
    cnres@gsnLeftString = ""
    cnres@gsnRightString = ""

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
  res@cnLevelSpacingF = 1 ;0.5
  res@cnMaxLevelValF = 10.
  res@cnMinLevelValF = -10.
  res@cnFillPalette = "CBR_drywet"
  ; cnres@cnLevelSelectionMode = "ExplicitLevels"
  ; cnres@cnLevels = (/12460,12470,12480,12490,12500,15210,12520,12540/);12520,12560
  cnres@cnLevelSelectionMode = "ManualLevels"
  cnres@cnLevelSpacingF = 5.
  cnres@cnMaxLevelValF = 20.
  cnres@cnMinLevelValF = 0.

  do i = 0,5
  res@gsnLeftString = month(i)
  
  plot(i) = gsn_csm_contour_map(wks, diff_z(i,:,:), res)
  plotc(i) = gsn_csm_contour(wks, z(0,i,:,:), cnres)
  ;overlay(plot(i), plotc(i))
  plotv(i) = gsn_csm_vector(wks, diff_u(i,:,:), diff_v(i,:,:), vcres)
  overlay(plot(i), plotv(i))
  plot_topo(i) = gsn_csm_contour(wks, tpmask, res_topo)
  overlay(plot(i), plot_topo(i))
end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
; pres@gsnPanelMainString = "UV200_H200"
;pres@pmLabelBarOrthogonalPosF = -0.01

gsn_panel(wks, plot, (/2,3/), pres)


end