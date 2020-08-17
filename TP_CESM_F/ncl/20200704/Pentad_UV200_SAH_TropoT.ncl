
begin

; lonstr = 50
; lonend = 180
; latstr = -10
; latend = 40
lonstr = 50
lonend = 150
latstr = 0 ;-10
latend = 50

casename = (/"CTRL","NOTPQ1"/) ;,"NOTPQ1_MAR","NOTPQ1_APR","NOTPQ1_MAY";"NOTP"
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

LEV = (/200/)
levs = 500
levt = 200

z = new((/ncase,nyear,npen,nlat,nlon/), float)
u = z
v = z
t = z

do i = 0, ncase-1
  diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(i)+"_F2000_F19G16_CAM4/pro/daily/"
  cname = "F2000_"+casename(i)+".PS.daily.06-30.nc"
  fc = addfile(diri+cname, "r")
  PS0 = fc->PS
  PS = reshape(PS0, (/nyear*365,nlat,nlon/))
  ; copy_VarCoords(PS0(0,0,:,:), PS(0,:,:))
  ; printVarSummary(PS)

  tname = systemfunc("ls "+diri+"F2000_"+casename(i)+".plevel.T.*.daily.nc")
  zname = systemfunc("ls "+diri+"F2000_"+casename(i)+".plevel.Z3.*.daily.nc")
  uname = systemfunc("ls "+diri+"F2000_"+casename(i)+".plevel.U.*.daily.nc")
  vname = systemfunc("ls "+diri+"F2000_"+casename(i)+".plevel.V.*.daily.nc")
  ft = addfiles(tname, "r")
  fz = addfiles(zname, "r")
  fu = addfiles(uname, "r")
  fv = addfiles(vname, "r")
  TX = ft[:]->T
  ZX = fz[:]->Z3(:,{LEV},:,:)
  UX = fu[:]->U(:,{LEV},:,:)
  VX = fv[:]->V(:,{LEV},:,:)

  copy_VarCoords(UX, PS)

  plev = tofloat(TX&lev)
  plev = plev*100.
  plev@units = "Pa"
  dp = dpres_plevel_Wrap(plev, PS, min(plev), 0)
  Tdp = TX*dp
  copy_VarCoords(TX, Tdp)
  tdp = dim_sum_n_Wrap(Tdp(:,{levs:levt},:,:), 1)/dim_sum_n_Wrap(dp(:,{levs:levt},:,:), 1)
  T = reshape(tdp, (/nyear,365,nlat,nlon/))
  Z = reshape(ZX, (/nyear,365,nlat,nlon/))  
  U = reshape(UX, (/nyear,365,nlat,nlon/))
  V = reshape(VX, (/nyear,365,nlat,nlon/))

  do p = 0, npen-1
    t(i,:,p,:,:) = dim_sum_n_Wrap(T(:,5*p:5*p+4,:,:), 1)/5.0
    z(i,:,p,:,:) = dim_sum_n_Wrap(Z(:,5*p:5*p+4,:,:), 1)/5.0
    u(i,:,p,:,:) = dim_sum_n_Wrap(U(:,5*p:5*p+4,:,:), 1)/5.0
    v(i,:,p,:,:) = dim_sum_n_Wrap(V(:,5*p:5*p+4,:,:), 1)/5.0
  end do
  copy_VarCoords(ZX(0,:,:), z(0,0,0,:,:))
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
copy_VarCoords(z, t)

;------------ttest--------------------------------------
pens = 24
pene = 37
pen = ispan(pens, pene, 1)
ux = u(:,:,pens-1:pene-1,:,:)
vx = v(:,:,pens-1:pene-1,:,:)
zx = z(:,:,pens-1:pene-1,:,:)
tx = t(:,:,pens-1:pene-1,:,:)
u_var = dim_variance_n_Wrap(ux, 1)
v_var = dim_variance_n_Wrap(vx, 1)
z_var = dim_variance_n_Wrap(zx, 1)
t_var = dim_variance_n_Wrap(tx, 1)

u_ave = dim_avg_n_Wrap(ux, 1)
v_ave = dim_avg_n_Wrap(vx, 1)
z_ave = dim_avg_n_Wrap(zx, 1)
t_ave = dim_avg_n_Wrap(tx, 1)

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
figname = "Pentad_UV200_SAH_TropoT"

wks = gsn_open_wks("pdf", figpath+figname)

plot = new((/ncase*6/), graphic)
plotv = plot
plotc = plot
; plotc2 = plot
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

  vcres@vcRefMagnitudeF         = 5.0  ; define vector ref mag

  res@cnLevelSelectionMode = "ManualLevels"
  res@cnLevelSpacingF = 1 ;0.5
  res@cnMaxLevelValF = 15.
  res@cnMinLevelValF = 1.
  res@cnFillPalette = "CBR_wet"
  cnres@cnLevelSelectionMode = "ExplicitLevels"
  cnres@cnLevels = (/5885/);12520,12560
  ; cnres@cnLevelSelectionMode = "ManualLevels"
  ; cnres@cnLevelSpacingF = 5.
  ; cnres@cnMaxLevelValF = 20.
  ; cnres@cnMinLevelValF = 0.

  do i = 0,ncase-1
  res@gsnRightString = casename(i)

  do j = 0,5
      res@gsnLeftString = "Pentad "+pen(2*j+1)
      plot(j*ncase+i) = gsn_csm_contour_map(wks, prec_ave(i,2*j+1,:,:), res)
      ; if (i.eq.0) then
      ;   cnres@cnLevels = (/5885/)
      ; end if
      ; if (i.eq.1) then
      ;   cnres@cnLevels = (/5890/)
      ; end if
      plotc(j*ncase+i) = gsn_csm_contour(wks, z_ave(i,2*j+1,:,{110:lonend}), cnres)
      overlay(plot(j*ncase+i), plotc(j*ncase+i))
      plotv(j*ncase+i) = gsn_csm_vector(wks, u_ave(i,2*j+1,:,:), v_ave(i,2*j+1,:,:), vcres)
      overlay(plot(j*ncase+i), plotv(j*ncase+i))
      plot_topo(j*ncase+i) = gsn_csm_contour(wks, tpmask, res_topo)
      overlay(plot(j*ncase+i), plot_topo(j*ncase+i))

      ; plot(i*5+j) = gsn_csm_contour_map(wks, prec_ave(i,2*j+1,:,:), res)
      ; plotc(i*5+j) = gsn_csm_contour(wks, z_ave(i,2*j+1,:,{110:lonend}), cnres)
      ; overlay(plot(i*5+j), plotc(i*5+j))
      ; plotv(i*5+j) = gsn_csm_vector(wks, u_ave(i,2*j+1,:,:), v_ave(i,2*j+1,:,:), vcres)
      ; overlay(plot(i*5+j), plotv(i*5+j))
      ; plot_topo(i*5+j) = gsn_csm_contour(wks, tpmask, res_topo)
      ; overlay(plot(i*5+j), plot_topo(i*5+j))
    end do
end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
; pres@gsnPanelMainString = "UV200_H200"
;pres@pmLabelBarOrthogonalPosF = -0.01

gsn_panel(wks, plot, (/6,ncase/), pres)


end