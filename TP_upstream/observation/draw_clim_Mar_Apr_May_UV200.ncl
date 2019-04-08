
begin

yrStrt = 1979
yrLast = 2013
year = ispan(yrStrt, yrLast, 1)

diri = "/home/yangsong3/data-observation/NCEP_DOE/monthly/"
filiu = "uwnd.mon.mean.nc"
filiv = "vwnd.mon.mean.nc"

fu = addfile(diri+filiu, "r")
fv = addfile(diri+filiv, "r")

TIME = fu->time
YYYY = cd_calendar(TIME, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

U = short2flt(fu->uwnd(iYYYY,{200},::-1,:)) ; -90 to 90
V = short2flt(fv->vwnd(iYYYY,{200},::-1,:)) ; -90 to 90
printVarSummary(U)

u = clmMonTLL(U)
v = clmMonTLL(V)
printVarSummary(u)

u3 = u(2,:,:) ;March
v3 = v(2,:,:)

u4 = u(3,:,:) ;April
v4 = v(3,:,:)

u5 = u(4,:,:) ;May
v5 = v(4,:,:)
printVarSummary(u3)

speed3 = wind_speed(u3, v3)
speed4 = wind_speed(u4, v4)
speed5 = wind_speed(u5, v5)
printVarSummary(speed3)

;==================================================================

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_upstream/climatology/"
figname = "draw_clim_Mar_Apr_May_UV200"
latstr = 0
latend = 70
lonstr = -30
lonend = 150

wks = gsn_open_wks("pdf", figpath+figname)
nplot = 3
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
 res@lbLabelBarOn = False
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

  res@gsnLeftString = "March"
  plot(0) = gsn_csm_contour_map(wks, speed3, res)
  vcid(0) = gsn_csm_vector(wks, u3(::4,::4), v3(::4,::4), vecres)
  plot_topo(0) = gsn_csm_contour(wks, topo, res_topo)
  overlay(plot(0), vcid(0))
  overlay(plot(0), plot_topo(0))

  res@gsnLeftString = "April"
  plot(1) = gsn_csm_contour_map(wks, speed4, res)
  vcid(1) = gsn_csm_vector(wks, u4(::4,::4), v4(::4,::4), vecres)
  plot_topo(1) = gsn_csm_contour(wks, topo, res_topo)
  overlay(plot(1), vcid(1))
  overlay(plot(1), plot_topo(1))

  res@gsnLeftString = "May"
  plot(2) = gsn_csm_contour_map(wks, speed5, res)
  vcid(2) = gsn_csm_vector(wks, u5(::4,::4), v5(::4,::4), vecres)
  plot_topo(2) = gsn_csm_contour(wks, topo, res_topo)
  overlay(plot(2), vcid(2))
  overlay(plot(2), plot_topo(2))


pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = "Clim 3-4-5 UV200 NCEP_DOE"
gsn_panel(wks, plot, (/3,1/), pres)
end