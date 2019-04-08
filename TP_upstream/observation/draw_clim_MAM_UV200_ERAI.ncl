
begin

yrStrt = 1979
yrLast = 2016
year = ispan(yrStrt, yrLast, 1)

diri = "/home/yangsong3/data-observation/ERA-interim/monthly/"
filiu = "uwnd.mon.mean.nc"
filiv = "vwnd.mon.mean.nc"

fu = addfile(diri+filiu, "r")
fv = addfile(diri+filiv, "r")

TIME = fu->time
YYYY = cd_calendar(TIME, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

U = short2flt(fu->u(iYYYY,{200},::-1,:)) ; -90 to 90
V = short2flt(fv->v(iYYYY,{200},::-1,:)) ; -90 to 90
printVarSummary(U)

u = clmMonTLL(U)
v = clmMonTLL(V)
printVarSummary(u)

uu = month_to_season(u, "MAM")
vv = month_to_season(v, "MAM")
printVarSummary(uu)

u0 = dim_avg_n_Wrap(uu, 0)
v0 = dim_avg_n_Wrap(vv, 0)
speed = wind_speed(u0, v0)
printVarSummary(speed)

;==================================================================

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_upstream/observation/"
figname = "draw_clim_MAM_UV200_ERAI"
latstr = 0
latend = 70
lonstr = -30
lonend = 150

wks = gsn_open_wks("pdf", figpath+figname)
nplot = 1
plot = new(nplot, graphic)
vcid = plot
plot_topo = plot

;------------------------------------------------------------
; TOPO
;------------------------------------------------------------
 
 f2 = "~/wuyt/data/TOPO_1500.nc"
 topo_in = addfile(f2,"r")
 topo = topo_in->TP

 xlat = topo_in->LATIXY(:,0)
 xlon = topo_in->LONGXY(0,:)
 topo!0 = "lsmlat"
 topo!1 = "lsmlon"
 topo&lsmlat = xlat
 topo&lsmlon = xlon
 
 res_topo            = True
 res_topo@cnFillOn   = False
 res_topo@cnLinesOn  = True
 res_topo@gsnFrame   = False
 res_topo@gsnDraw    = False
 res_topo@cnLineLabelsOn         = False
 res_topo@cnLineThicknessF       = 3
 res_topo@cnLineColor            = "gray30"
 res_topo@cnInfoLabelOn          = False
 res_topo@cnLevelSelectionMode   ="ExplicitLevels"
 res_topo@cnLevels               = (/1500/)
;------------------------------------------------------------
 res = True
 res@gsnDraw  = False
 res@gsnFrame = False
 res@gsnAddCyclic = True
 res@gsnMaximize = True

 res@mpFillOn = True
 res@mpMinLatF = latstr
 res@mpMaxLatF = latend
 res@mpMinLonF = lonstr
 res@mpMaxLonF = lonend
 res@mpGeophysicalLineThicknessF =2
 res@mpGeophysicalLineColor = "grey40"

 res@cnFillOn = True
 res@cnLinesOn = False
 res@cnLineLabelsOn = False
 ;res@lbLabelBarOn = False
 res@cnFillPalette = "BlueWhiteOrangeRed"

 res@tmYLLabelFontHeightF = 0.02
 res@tmXBLabelFontHeightF = 0.02
 res@gsnStringFontHeightF = 0.03
 res@gsnRightString = "m/s"

 res@cnLevelSelectionMode = "ManualLevels"  ;ExplicitLevels
 res@cnMaxLevelValF = 60.
 res@cnMinLevelValF = 10.
 res@cnLevelSpacingF = 5.

;************************************************
; setting : vector
;************************************************

  vecres                         = True            ; vector only resources
  vecres@gsnDraw                 = False           ; don't draw
  vecres@gsnFrame                = False           ; don't advance frame
  vecres@vcGlyphStyle            = "LineArrow" ;"CurlyVector"   ; curly vectors
  vecres@vcRefMagnitudeF         = 30 ;1.3               ; define vector ref mag
  vecres@vcRefLengthF            = 0.045           ; define length of vec ref
;  vecres@vcRefAnnoOn             = False
  vecres@gsnRightString          = " "             ; turn off right string
  vecres@gsnLeftString           = " "             ; turn off left string
  vecres@tiXAxisString           = " "             ; turn off axis label
  vecres@vcRefAnnoString2On      = False
  vecres@vcRefAnnoOrthogonalPosF = -1.    ; move ref vector into plot
  vecres@vcRefAnnoFontHeightF    = .025    ; move ref vector into plot
  vecres@vcLineArrowHeadMaxSizeF      = 0.01
  vecres@vcLineArrowHeadMinSizeF      = 0.008
  vecres@vcLineArrowThicknessF      = 3

  vecres@vcMonoLineArrowColor    = True

  res@gsnLeftString = "Clim MAM UV200 ERA-Interm"
  plot(0) = gsn_csm_contour_map(wks, speed, res)
  vcid(0) = gsn_csm_vector(wks, u0(::8,::8), v0(::8,::8), vecres)
  plot_topo(0) = gsn_csm_contour(wks, topo, res_topo)
  overlay(plot(0), vcid(0))
  overlay(plot(0), plot_topo(0))

draw(plot)
end