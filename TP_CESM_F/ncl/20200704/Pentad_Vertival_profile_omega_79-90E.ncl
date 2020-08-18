
begin

; lonstr = 50
; lonend = 180
; latstr = -10
; latend = 40
lonstr = 50
lonend = 150
latstr = 0 ;-10
latend = 50

casename = (/"CTRL","NOTPQ1","NOTPQ1_MAR","NOTPQ1_APR","NOTPQ1_MAY"/) ;;"NOTP"
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
nlev = 37

lonW = 70
lonE = 90

u = new((/ncase,nyear,npen,nlev,nlat/), float)
v = u

do i = 0, ncase-1
  diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(i)+"_F2000_F19G16_CAM4/pro/daily/"

  uname = systemfunc("ls "+diri+"F2000_"+casename(i)+".plevel.U.*.daily.nc")
  vname = systemfunc("ls "+diri+"F2000_"+casename(i)+".plevel.OMEGA.*.daily.nc")
  fu = addfiles(uname, "r")
  fv = addfiles(vname, "r")
  UX = fu[:]->U ;(:,{LEV},:,:)
  VX = fv[:]->OMEGA ;(:,{LEV},:,:)
  VX = VX*100.*-1.
  printVarSummary(UX)
  
  UX0 = dim_avg_n_Wrap(UX(:,:,:,{lonW:lonE}), 3)
  VX0 = dim_avg_n_Wrap(VX(:,:,:,{lonW:lonE}), 3)

  U = reshape(UX0, (/nyear,365,nlev,nlat/))
  V = reshape(VX0, (/nyear,365,nlev,nlat/))

  do p = 0, npen-1
    u(i,:,p,:,:) = dim_sum_n_Wrap(U(:,5*p:5*p+4,:,:), 1)/5.0
    v(i,:,p,:,:) = dim_sum_n_Wrap(V(:,5*p:5*p+4,:,:), 1)/5.0
  end do
  copy_VarCoords(UX(0,:,:,0), u(0,0,0,:,:))
end do

; latx = tofloat(UX&lat)
; lonx = tofloat(UX&lon )

u!0 = "case"
u&case = casename
u!1 = "year"
u&year = year
u!2 = "pentad"
u&pentad = pentad
; z!3 = "lat"
; z&lat = lat
; z!4 = "lon"
; z&lon = lon
printVarSummary(u)
copy_VarCoords(u, v)

;------------ttest--------------------------------------
pens = 24
pene = 37
pen = ispan(pens, pene, 1)
ux = u(:,:,pens-1:pene-1,:,:)
vx = v(:,:,pens-1:pene-1,:,:)
u_var = dim_variance_n_Wrap(ux, 1)
v_var = dim_variance_n_Wrap(vx, 1)

u_ave = dim_avg_n_Wrap(ux, 1)
v_ave = dim_avg_n_Wrap(vx, 1)

probu = u_ave(0:ncase-1,:,:,:)
probv = probu
diff_u = u_ave(0:ncase-1,:,:,:)
diff_v = diff_u
sig = 0.05

do k = 1, ncase-1
  probu(k,:,:,:) = ttest(u_ave(0,:,:,:), u_var(0,:,:,:), nyear, u_ave(k,:,:,:), u_var(k,:,:,:), nyear, False, False)
  probv(k,:,:,:) = ttest(v_ave(0,:,:,:), v_var(0,:,:,:), nyear, v_ave(k,:,:,:), v_var(k,:,:,:), nyear, False, False)

  ; diff_z(k,:,:,:) = z_ave(0,:,:,:)
  diff_u(k,:,:,:) = u_ave(0,:,:,:) - u_ave(k,:,:,:)
  diff_v(k,:,:,:) = v_ave(0,:,:,:) - v_ave(k,:,:,:)


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
topo_m = dim_avg_n_Wrap(tpmask(:,{lonW:lonE}), 1)
topo_merid = topo_m
topo_merid = 1013.25*(1-topo_m*0.0065/288.15)^5.25145
copy_VarCoords(tpmask(0,:), topo_merid)
printVarSummary(topo_merid)
;-------------------------------------------------------------------
;-------------------------------------------------------------------

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20200704/"
figname = "Pentad_Vertical_profile_omega_70-90E"

wks = gsn_open_wks("pdf", figpath+figname)

plot = new((/(ncase-1)*5/), graphic)
plotv = plot
plotz1 = plot
plotz2 = plot
plot_topo = plot

;=====================================================
    res_topo                   = True

    res_topo@gsnDraw  = False           ;先别直接画出来
    res_topo@gsnFrame = False           ;也不要翻页

    res_topo@gsnLeftString  = " "              ;不要标示出变量名和单位等信息
    res_topo@gsnRightString = " "

    ;res_topo@gsnMaximize        = True
    ;res_topo@trYMaxF            = 1000.0
    res_topo@gsnYRefLine        = 1000 ;(/100,1000/)          ; create a reference line 
    res_topo@gsnBelowYRefLineColor = "black" ;(/"black","black"/) ; above ref line fill grey
    ;res_topo@gsnXYBelowFillColors = "goldenrod2"
    res_topo@gsnAboveYRefLineColor = "blue"
    res_topo@xyLineColors = "black"
    res_topo@xyLineThicknesses = 4.
;=====================================================

;=====================================================

    res = True
    res@gsnDraw  = False
    res@gsnFrame = False
    res@vpHeightF = 0.4
    res@vpWidthF  = 0.7

    res@cnFillOn = True
    res@cnLinesOn = False
    res@cnLineLabelsOn = False
    res@cnInfoLabelOn = False

    res@gsnStringFontHeightF = 0.03
    res@tmXBLabelFontHeightF = 0.02
    res@tmYLLabelFontHeightF = 0.02

    ;res@pmLabelBarOrthogonalPosF = 0.12
    res@lbLabelBarOn = False
    ;res@lbLabelFontHeightF = 0.02
    res@cnFillPalette = "BlueWhiteOrangeRed"
    res@tiYAxisString = ""


 ;------------------------------------------
     res@trYMinF   = 100
     res@trYMaxF   = 1000
      
     res@tmYLMode      = "Explicit"
     res@tmYLValues    = (/"1000","850","750","500","300","200","100"/)
     res@tmYLLabels    = (/"1000","850","750","500","300","200","100"/)

     res@vcRefAnnoOrthogonalPosF =  -0.35 ;-1.0             ; move ref vector up
     ;res@vcGlyphStyle            = "CurlyVector"  ;"LineArrow"    ; turn on curly vectors
     ;res@vcMinDistanceF          = 0.017
     res@vcGlyphStyle            = "LineArrow"
     res@vcMinDistanceF          = 0.03
     res@vcRefLengthF            = 0.04            ; define length of vec ref
     
     res@vcLineArrowThicknessF   = 3.0
     ;res@vcLineArrowHeadMinSizeF = 0.008
     ;res@vcLineArrowHeadMaxSizeF = 0.3
     res@vcVectorDrawOrder = "PostDraw" 
     
     res@vcRefAnnoString2On = False
     res@vcRefAnnoFontHeightF = 0.02

    res@cnLevelSelectionMode = "ManualLevels"
    res@cnLevelSpacingF = 1
    res@cnMaxLevelValF = 10
    res@cnMinLevelValF = -10
    res@cnFillPalette = "BlueWhiteOrangeRed"
    res@vcRefMagnitudeF         = 10.             ; define vector ref mag

    res@gsnRightString = lonW+"~S~o~N~E-"+lonE+"~S~o~N~E"

  n = ncase-1
  do i = 0,ncase-2
  res@gsnRightString = "Ctrl -"+casename(i+1)

  do j = 0,4
      res@gsnLeftString = "Pentad "+pen(2*j+1)
      plot(j*n+i) = gsn_csm_pres_hgt_vector(wks, diff_v(i+1,2*j+1,{1000:100},{0:60}),\
 diff_u(i+1,2*j+1,{1000:100},{0:60}), diff_v(i+1,2*j+1,{1000:100},{0:60}), res)

      ; plotz1(j*ncase+i) = gsn_csm_contour(wks, z_ave(i,2*j+1,:,:), cnres)

      ; plotv(j*ncase+i) = gsn_csm_vector(wks, u_ave(i,2*j+1,:,:), v_ave(i,2*j+1,:,:), vcres)
      ; overlay(plot(j*ncase+i), plotv(j*ncase+i))
      ; plot_topo(j*ncase+i) = gsn_csm_xy(wks, lat(48:79), topo_merid({0:60}), res_topo)
      plot_topo(j*n+i) = gsn_csm_y(wks, topo_merid({0:60}), res_topo)
      overlay(plot(j*n+i), plot_topo(j*n+i))

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

gsn_panel(wks, plot, (/5,ncase-1/), pres)


end