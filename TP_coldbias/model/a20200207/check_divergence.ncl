
begin

yrStrt = 1979
yrLast = 2016
lonstr = 60 ;10
lonend = 110 ;150
latstr = 20 ;-10
latend = 50 ;70

diri  = "/home/yangsong3/data-observation/ERA-interim/monthly/"
uname = "uwnd.mon.mean.nc"
vname = "vwnd.mon.mean.nc"

fu = addfile(diri+uname, "r")
fv = addfile(diri+vname, "r")

TIME = fu->time
YYYY = cd_calendar(TIME, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

u1 = short2flt(fu->u(iYYYY,:,{latstr:latend},{lonstr:lonend}))
v1 = short2flt(fv->v(iYYYY,:,{latstr:latend},{lonstr:lonend}))
printVarSummary(u1)

lat = u1&latitude
lon = u1&longitude

dv = uv2dv_cfd(u1, v1, lat, lon, 0)
copy_VarCoords(u1, dv)

season = "JJA"
dv0 = month_to_season(dv, season)
dv_JJA = dim_avg_n_Wrap(dv0, 0)

 ;----------------------------------------------------------
 wks = gsn_open_wks("pdf", "check_divergence")
 plot = new(1, graphic)

  res   = True
  res@gsnDraw  = False
  res@gsnFrame = False
  res@gsnAddCyclic = False

  res@cnFillOn      = True
  res@cnLinesOn     = False
  res@cnInfoLabelOn = False
  res@cnLineLabelsOn = False

  res@mpOutlineOn   = True
  res@mpFillOn      = False
  res@mpMinLatF     = latstr
  res@mpMaxLatF     = latend
  res@mpMinLonF     = lonstr
  res@mpMaxLonF     = lonend
  res@mpGeophysicalLineThicknessF = 4
  res@gsnRightStringFontHeightF = 0.035

  res@gsnStringFontHeightF = 0.04

  res@tmXBLabelFontHeightF = 0.025
  res@tmYLLabelFontHeightF = 0.025

  res@lbLabelFontHeightF   = 0.015
  res@pmLabelBarHeightF    = 0.1
  res@pmLabelBarWidthF     = 0.75
  res@pmLabelBarOrthogonalPosF = 0.12
  ;res@o=pmLabelBarParallelPosF = 

  res@gsnSpreadColors = True
  res@gsnSpreadColorEnd   = 0
  res@gsnSpreadColorStart = 127
  res@cnFillPalette = "MPL_BrBG" ;CBR_drywet
  ; res@cnLevelSelectionMode = "ManualLevels"
  ; res@gsnRightString = "g/m~S~2~N~/s"
  ; res@cnMaxLevelValF  = 0.05
  ; res@cnMinLevelValF  = -0.05
  ; res@cnLevelSpacingF = 0.005

  plot(0) = gsn_csm_contour_map(wks, dv_JJA({500},:,:), res)

draw(plot)


end