
begin

; lonstr = 50
; lonend = 180
; latstr = -10
; latend = 40
lonstr = 50
lonend = 150
latstr = 0 ;-10
latend = 50

casename = (/"CTRL","NOTPQ1_MAR","NOTPQ1_APR","NOTPQ1_MAY"/) ;"NOTP""NOTPQ1",
ncase = dimsizes(casename)
pentad = ispan(1, 73, 1)
npen = dimsizes(pentad)
lat = fspan(-90, 90, 96)
lat@unit = "degrees_north"
lon = fspan(0, 357.5, 144)
lon@units = "degrees_east"
nlat = dimsizes(lat)
nlon = dimsizes(lon)
year = ispan(6, 30, 1)
nyear = dimsizes(year)

LEV = (/850,500/)

z = new((/ncase,nyear,npen,nlat,nlon/), float)
u = z
v = z
prec = z

do i = 0, ncase-1
  diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(i)+"_F2000_F19G16_CAM4/pro/daily/"
  cname = "F2000_"+casename(i)+".PRECC.daily.06-30.nc"
  lname = "F2000_"+casename(i)+".PRECL.daily.06-30.nc"
  fc = addfile(diri+cname, "r")
  fl = addfile(diri+lname, "r")
  PRECL = fl->PRECL ;(4:13,:,:,:)
  PRECC = fc->PRECC ;(4:13,:,:,:)
  PREC = (PRECC+PRECL)*86400*1000
  copy_VarCoords(PRECC, PREC)
  PREC@units = "mm/day"
  printVarSummary(PREC)


  zname = systemfunc("ls "+diri+"F2000_"+casename(i)+".plevel.Z3.*.daily.nc")
  uname = systemfunc("ls "+diri+"F2000_"+casename(i)+".plevel.U.*.daily.nc")
  vname = systemfunc("ls "+diri+"F2000_"+casename(i)+".plevel.V.*.daily.nc")
  fz = addfiles(zname, "r")
  fu = addfiles(uname, "r")
  fv = addfiles(vname, "r")
  ZX = fz[:]->Z3(:,{LEV(1)},:,:)
  UX = fu[:]->U(:,{LEV(0)},:,:)
  VX = fv[:]->V(:,{LEV(0)},:,:)
  Z = reshape(ZX, (/nyear,365,nlat,nlon/))  
  U = reshape(UX, (/nyear,365,nlat,nlon/))
  V = reshape(VX, (/nyear,365,nlat,nlon/))

  do p = 0, npen-1
    prec(i,:,p,:,:) = dim_sum_n_Wrap(PREC(:,5*p:5*p+4,:,:), 1)/5.0
    z(i,:,p,:,:) = dim_sum_n_Wrap(Z(:,5*p:5*p+4,:,:), 1)/5.0
    u(i,:,p,:,:) = dim_sum_n_Wrap(U(:,5*p:5*p+4,:,:), 1)/5.0
    v(i,:,p,:,:) = dim_sum_n_Wrap(V(:,5*p:5*p+4,:,:), 1)/5.0
  end do
  copy_VarCoords(PRECL(0,0,:,:), z(0,0,0,:,:))
end do

z!0 = "case"
z&case = casename
z!1 = "year"
z&year = year
z!2 = "pentad"
z&pentad = pentad
; z!3 = "lat"
; z&lat = lat
; z!4 = "lon"
; z&lon = lon
printVarSummary(z)
copy_VarCoords(z, u)
copy_VarCoords(z, v)
copy_VarCoords(z, prec)

;------------ttest--------------------------------------
pens = 24
pene = 34
pen = ispan(pens, pene, 1)
ux = u(:,:,pens-1:pene-1,:,:)
vx = v(:,:,pens-1:pene-1,:,:)
zx = z(:,:,pens-1:pene-1,:,:)
precx = prec(:,:,pens-1:pene-1,:,:)
u_var = dim_variance_n_Wrap(ux, 1)
v_var = dim_variance_n_Wrap(vx, 1)
z_var = dim_variance_n_Wrap(zx, 1)
prec_var = dim_variance_n_Wrap(precx, 1)

u_ave = dim_avg_n_Wrap(ux, 1)
v_ave = dim_avg_n_Wrap(vx, 1)
z_ave = dim_avg_n_Wrap(zx, 1)
prec_ave = dim_avg_n_Wrap(precx, 1)

probu = u_ave(0:ncase-1,:,:,:)
probv = probu
probp = probu
diff_u = u_ave(0:ncase-1,:,:,:)
diff_v = diff_u
diff_prec = diff_u
sig = 0.05

do k = 1, ncase-1
  probu(k,:,:,:) = ttest(u_ave(0,:,:,:), u_var(0,:,:,:), nyear, u_ave(k,:,:,:), u_var(k,:,:,:), nyear, False, False)
  probv(k,:,:,:) = ttest(v_ave(0,:,:,:), v_var(0,:,:,:), nyear, v_ave(k,:,:,:), v_var(k,:,:,:), nyear, False, False)
  ; probz = ttest(z_ave(0,:,:,:), z_var(0,:,:,:), nyear, z_ave(1,:,:,:), z_var(1,:,:,:), nyear, False, False)
  probp(k,:,:,:) = ttest(prec_ave(0,:,:,:), prec_var(0,:,:,:), nyear, prec_ave(k,:,:,:), prec_var(k,:,:,:), nyear, False, False)

  ; diff_z = z_ave(0,:,:,:)
  ; diff_z = z_ave(0,:,:,:) - z_ave(k,:,:,:)
  diff_u(k,:,:,:) = u_ave(0,:,:,:) - u_ave(k,:,:,:)
  diff_v(k,:,:,:) = v_ave(0,:,:,:) - v_ave(k,:,:,:)
  diff_prec(k,:,:,:) = prec_ave(0,:,:,:) - prec_ave(k,:,:,:)


  ; diff_u(k,:,:,:) = mask(diff_u(k,:,:,:), (probu(k,:,:,:).lt.sig).or.(probv(k,:,:,:).lt.sig), True)
  ; diff_v(k,:,:,:) = mask(diff_v(k,:,:,:), (probu(k,:,:,:).lt.sig).or.(probv(k,:,:,:).lt.sig), True)
  ; diff_prec(k,:,:,:) = mask(diff_prec(k,:,:,:), probp(k,:,:,:).lt.sig, True)
end do
printVarSummary(diff_u)
;----------------------------------------------------------------
tpmaskpath = "/home/yangsong3/wuyt/sysu/scripts/TP_CESM_F/model_need_input_data/"
; tpmaskname = "TPmask_1000m_C20190410.nc"
tpmaskname = "TPmask_1500m_C20200713.nc"
ftp = addfile(tpmaskpath+tpmaskname, "r")

tpmask = ftp->tp
printVarSummary(tpmask)
;-------------------------------------------------------------------
;-------------------------------------------------------------------

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20200704/"
figname = "Pentad_UV850_Precip_WPSH_diff"

wks = gsn_open_wks("pdf", figpath+figname)

plot = new((/(ncase-1)*5/), graphic)
plotv = plot
plotc1 = plot
plotc2 = plot
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
  cnres@cnLevelSelectionMode = "ExplicitLevels"
  cnres@cnLevels = (/5885/);12520,12560
  ; cnres@cnLevelSelectionMode = "ManualLevels"
  ; cnres@cnLevelSpacingF = 5.
  ; cnres@cnMaxLevelValF = 20.
  ; cnres@cnMinLevelValF = 0.

  do i = 0,ncase-2
  ; res@gsnLeftString = month(i+3)
  
  ; plot(i) = gsn_csm_contour_map(wks, diff_z(i+3,:,:), res)
  ; plotc(i) = gsn_csm_contour(wks, z_ave(0,i+3,:,:), cnres)
  ; ;overlay(plot(i), plotc(i))
  ; plotv(i) = gsn_csm_vector(wks, diff_u(i+3,:,:), diff_v(i+3,:,:), vcres)
  ; overlay(plot(i), plotv(i))
  ; plot_topo(i) = gsn_csm_contour(wks, tpmask, res_topo)
  ; overlay(plot(i), plot_topo(i))

  res@gsnRightString = "Ctrl - "+casename(i+1)
    do j = 0,4
      res@gsnLeftString = "Pentad "+pen(j+1)
      ; plot(j*4+i) = gsn_csm_contour_map(wks, diff_prec(i+1,2*j+1,:,:), res)
      ; cnres@cnLineColor = "red"
      ; plotc1(j*4+i) = gsn_csm_contour(wks, z_ave(0,2*j+1,:,{105:lonend}), cnres)
      ; overlay(plot(j*4+i), plotc1(j*4+i))
      ; cnres@cnLineColor = "blue"
      ; plotc2(j*4+i) = gsn_csm_contour(wks, z_ave(i+1,2*j+1,:,{105:lonend}), cnres)
      ; overlay(plot(j*4+i), plotc2(j*4+i))
      ; plotv(j*4+i) = gsn_csm_vector(wks, diff_u(i+1,2*j+1,:,:), diff_v(i+1,2*j+1,:,:), vcres)
      ; overlay(plot(j*4+i), plotv(j*4+i))
      ; plot_topo(j*4+i) = gsn_csm_contour(wks, tpmask, res_topo)
      ; overlay(plot(j*4+i), plot_topo(j*4+i))
      plot(j+i*5) = gsn_csm_contour_map(wks, diff_prec(i+1,j+1,:,:), res)
      cnres@cnLineColor = "red"
      plotc1(j+i*5) = gsn_csm_contour(wks, z_ave(0,j+1,:,{105:lonend}), cnres)
      overlay(plot(j+i*5), plotc1(j+i*5))
      cnres@cnLineColor = "blue"
      plotc2(j+i*5) = gsn_csm_contour(wks, z_ave(i+1,j+1,:,{105:lonend}), cnres)
      overlay(plot(j+i*5), plotc2(j+i*5))
      plotv(j+i*5) = gsn_csm_vector(wks, diff_u(i+1,j+1,:,:), diff_v(i+1,j+1,:,:), vcres)
      overlay(plot(j+i*5), plotv(j+i*5))
      plot_topo(j+i*5) = gsn_csm_contour(wks, tpmask, res_topo)
      overlay(plot(j+i*5), plot_topo(j+i*5))
    end do
end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
; pres@gsnPanelMainString = "UV200_H200"
pres@pmLabelBarOrthogonalPosF = -0.02
pres@lbLabelFontHeightF = 0.015

gsn_panel(wks, plot, (/3,5/), pres)


end