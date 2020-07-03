
begin

lonstr = 40
lonend = 130
latstr = 0
latend = 70

yrStrt = 1979
yrLast = 2010

month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
nm = dimsizes(month)
season = "DJF"

LEV = (/850,200/)
u = new((/2,181,360/), float)
v = u
precp = new((/181,360/), float)
z = precp


diri = "/home/yangsong3/data-observation/ERA-interim/monthly/"
filiu = "uwnd.mon.mean.nc"
filiv = "vwnd.mon.mean.nc"
filiz = "hgt.mon.mean.nc"
filic = "2_metre_temperature.mon.mean.nc"


fu = addfile(diri+filiu, "r")
fv = addfile(diri+filiv, "r")

TIME = fu->time
YYYY = cd_calendar(TIME, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

u1 = short2flt(fu->u(iYYYY,{LEV},:,:))
u_zonal = dim_avg_n_Wrap(u1, 2)
u1 = u1 - conform_dims(dimsizes(u1), u_zonal, (/0,1,3/))
v1 = short2flt(fv->v(iYYYY,{LEV},:,:))
v_zonal = dim_avg_n_Wrap(v1, 2)
v1 = v1 - conform_dims(dimsizes(v1), v_zonal, (/0,1,3/))
printVarSummary(u1)

fz = addfile(diri+filiz, "r")
z1 = short2flt(fz->z(iYYYY,{500},:,:))
z_zonal = dim_avg_n_Wrap(z1, 1)
z1 = z1 - conform_dims(dimsizes(z1), z_zonal, (/0,2/))

fc = addfile(diri+filic, "r")
c1 = short2flt(fc->t2m(iYYYY,:,:))
c_zonal = dim_avg_n_Wrap(c1, 1)
c1 = c1 - conform_dims(dimsizes(c1), c_zonal, (/0,2/))
PRECP = c1

U1 = month_to_season(u1, season)
u = dim_avg_n_Wrap(U1, 0)
V1 = month_to_season(v1, season)
v = dim_avg_n_Wrap(V1, 0)
PRECP1 = month_to_season(PRECP, season)
precp = dim_avg_n_Wrap(PRECP1, 0)
Z1 = month_to_season(z1, season)
z = dim_avg_n_Wrap(Z1, 0)


; ux = u(0:1,:,:,:)
; vx = v(0:1,:,:,:)
; px = precp(0:1,:,:)
; zx = z(0:1,:,:)

; ux(0,:,:,:) = u(0,:,:,:) - u(1,:,:,:)
; ux(1,:,:,:) = u(1,:,:,:) - u(2,:,:,:)
; vx(0,:,:,:) = v(0,:,:,:) - v(1,:,:,:)
; vx(1,:,:,:) = v(1,:,:,:) - v(2,:,:,:)
; px(0,:,:) = precp(0,:,:) - precp(1,:,:)
; px(1,:,:) = precp(1,:,:) - precp(2,:,:)
; zx(0,:,:) = z(0,:,:) - z(1,:,:)
; zx(1,:,:) = z(1,:,:) - z(2,:,:)
; printVarSummary(ux)

;----------------------------------------------------------------
tpmaskpath = "/home/yangsong3/wuyt/data/"
tpmaskname = "topo.sfc.1x1.nc"
ftp = addfile(tpmaskpath+tpmaskname, "r")

tpmask = ftp->topo
printVarSummary(tpmask)
u(0,:,:) = mask(u(0,:,:), tpmask.ge.1500, False)
v(0,:,:) = mask(v(0,:,:), tpmask.ge.1500, False)
;==================================================================
;==================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20191206/4season/"


  figname = "ERAI_circulation_"+season

  wks = gsn_open_wks("pdf", figpath+figname)
  nplot = 2
  plot = new(nplot, graphic)
  plot_topo = plot
  plotv = plot


   res_topo            = True
   res_topo@cnFillOn   = False
   res_topo@cnLinesOn  = True
   res_topo@gsnFrame   = False
   res_topo@gsnDraw    = False
   res_topo@cnLineLabelsOn         = False
   res_topo@cnLineThicknessF       = 5
   res_topo@cnLineColor            = "grey55"
   res_topo@cnLineDashPattern      = 16
   res_topo@cnInfoLabelOn          = False
   res_topo@cnLevelSelectionMode   ="ExplicitLevels"
   res_topo@cnLevels               = (/1500/)
   res_topo@gsnLeftString = ""
   res_topo@gsnRightString = ""


  res = True
  res@gsnDraw  = False
  res@gsnFrame = False

  res@mpOutlineOn = True
  res@mpFillOn = True
  res@mpGeophysicalLineThicknessF = 2
  res@mpGeophysicalLineColor = "grey45" ;"grey25"
  res@mpFillDrawOrder = "PreDraw"

  res@mpMinLatF = latstr
  res@mpMaxLatF = latend
  res@mpMinLonF = lonstr
  res@mpMaxLonF = lonend
  res@mpCenterLonF = (lonstr+lonend)/2

  res@cnFillOn = True
  res@cnLinesOn = False
  res@cnLineLabelsOn = False
  res@cnInfoLabelOn = False

  res@gsnStringFontHeightF = 0.035
  res@tmXBLabelFontHeightF = 0.025
  res@tmYLLabelFontHeightF = 0.025
  res@pmLabelBarOrthogonalPosF = 0.25
  ;------------------------------------------

    vcres                         = True             ; plot mods desired
    vcres@gsnFrame                = False            ; so we can draw time stamp
    vcres@gsnDraw                 = False
    vcres@vcRefAnnoOrthogonalPosF =  -1.0             ; move ref vector up

    ;vcres@vcGlyphStyle            = "CurlyVector"  ;"LineArrow"    ; turn on curly vectors
    ;vcres@vcMinDistanceF          = 0.017
    vcres@vcGlyphStyle            = "LineArrow"
        
    vcres@vcLineArrowThicknessF   = 3.0
    ;vecres@vcLineArrowHeadMinSizeF = 0.008
    ;vecres@vcLineArrowHeadMaxSizeF = 0.3

    vcres@vcRefAnnoString2On = False
    vcres@vcRefAnnoFontHeightF = 0.02
    vcres@gsnLeftString = ""
    vcres@gsnRightString = ""
  ;------------------------------------------
  
    ;res@lbLabelBarOn = False
    ;res@cnLevelSelectionMode = "ManualLevels"
    res@gsnLeftString = "UV"+LEV(0)+"_Precip"
 
  res@gsnRightString = "ERAI"
  vcres@vcRefMagnitudeF         = 5.0  ; define vector ref mag
  vcres@vcMinDistanceF          = 0.04
  vcres@vcRefLengthF            = 0.05            ; define length of vec ref

  res@cnLevelSpacingF = 2.
  res@cnMaxLevelValF = 20.
  res@cnMinLevelValF = 0.
  res@cnFillPalette = "hotcold_18lev"
  plot(0) = gsn_csm_contour_map(wks, precp, res)
  plot_topo(0) = gsn_csm_contour(wks, tpmask, res_topo)
  plotv(0) = gsn_csm_vector(wks, u({LEV(0)},:,:), v({LEV(0)},:,:), vcres)
  overlay(plot(0), plot_topo(0))
  overlay(plot(0), plotv(0))

;-------------------------------------------------------------------
  
  res@gsnLeftString = "UV"+LEV(1)+"_H500"

  res@gsnRightString = "ERAI"
  vcres@vcRefMagnitudeF         = 20.0  ; define vector ref mag
  res@cnLevelSpacingF = 5.
  res@cnMaxLevelValF = 420.
  res@cnMinLevelValF = 380.
  res@cnFillPalette = "WhiteYellowOrangeRed"
  plot(1) = gsn_csm_contour_map(wks, z, res)
  plot_topo(1) = gsn_csm_contour(wks, tpmask, res_topo)
  plotv(1) = gsn_csm_vector(wks, u({LEV(1)},:,:), v({LEV(1)},:,:), vcres)
  overlay(plot(1), plot_topo(1))
  overlay(plot(1), plotv(1))

 
  pres = True
  pres@gsnMaximize = True
  pres@gsnPanelMainString = season
  ;pres@gsnPanelLabelBar = True
  ;pres@pmLabelBarOrthogonalPosF = -0.01

  gsn_panel(wks, plot, (/2,1/), pres)



end