
begin

lonstr = -20
lonend = 160
latstr = 0
latend = 70

yrStrt = 1979
yrLast = 2019
year = ispan(yrStrt, yrLast, 1)
nyear = dimsizes(year)
season = "JJA"

LEVs = 1000
LEVe = 300

latstr_east = 30
latend_east = 45
lonstr_east = 105
lonend_east = 140

latstr_west = 30
latend_west = 55
lonstr_west = -10
lonend_west = 60

latstr_south = 10
latend_south = 28
lonstr_south = 70
lonend_south = 90

latstr_north = 45
latend_north = 60
lonstr_north = 85
lonend_north = 120 

latstr_tp = 25 ;28
latend_tp = 40 ;38
lonstr_tp = 75 ;80
lonend_tp = 105

region = (/"West","East","South","North","TP"/)
;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)
printVarSummary(topo)


;===================================================================
diri = "/home/yangsong3/data-observation/ERA5-monthly/pressure/"
fili1 = "cloud_ice_water_content_mon_1x1_1979_2019.nc" ;"temp.monthly.197901-201812.nc"
fili2 = "cloud_liquid_water_content_mon_1x1_1979_2019.nc"
f1 = addfile(diri+fili1, "r")
f2 = addfile(diri+fili2, "r")

Time = f1->time
YYYY = cd_calendar(Time, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

var1 = short2flt(f1->ciwc(iYYYY,::-1,::-1,:))
var2 = short2flt(f2->clwc(iYYYY,::-1,::-1,:))
t = var1
t = var1+var2
printVarSummary(t)

lat = t&latitude
lon = t&longitude
nlat = dimsizes(lat)
nlon = dimsizes(lon)

filiq = "q_mon_1x1_1979_2019.nc" ;"temp.monthly.197901-201812.nc"
filiu = "u_mon_1x1_1979_2020.nc"
filiv = "v_mon_1x1_1979_2020.nc"

fq = addfile(diri+filiq, "r")
fu = addfile(diri+filiu, "r")
fv = addfile(diri+filiv, "r")
q = short2flt(fq->q(iYYYY,::-1,::-1,:))
u = short2flt(fu->u(iYYYY,::-1,::-1,:))
v = short2flt(fv->v(iYYYY,::-1,::-1,:))

filips = "/home/yangsong3/data-observation/ERA5-monthly/surface/sur_pressure_mon_1x1_1979_2019.nc"
fps = addfile(filips, "r")
ps = short2flt(fps->sp(iYYYY,::-1,:))
printVarSummary(ps)


plev = tofloat(t&level)
plev = plev*100
plev!0 = "lev"
plev&lev = t&level
plev@units = "Pa"
; print(plev)

dp = dpres_plevel_Wrap(plev, ps, min(plev), 0)
printVarSummary(dp)

tdp = t*dp
copy_VarCoords(t, tdp)
printVarSummary(tdp)
qudp = q*u*dp
qvdp = q*v*dp
copy_VarCoords(q, qudp)
copy_VarCoords(q, qvdp)
printVarSummary(qudp)

t_tropo = dim_sum_n_Wrap(tdp(:,{LEVs:LEVe},:,:), 1)
t_tropo = t_tropo/9.81
copy_VarCoords(t(:,0,:,:), t_tropo)
t2m_clm = month_to_season(t_tropo, season)
printVarSummary(t2m_clm)

qu = dim_sum_n_Wrap(qudp(:,{LEVs:LEVe},:,:), 1)
qv = dim_sum_n_Wrap(qvdp(:,{LEVs:LEVe},:,:), 1)
qu = qu/9.81
qv = qv/9.81
copy_VarCoords(q(:,0,:,:), qu)
copy_VarCoords(q(:,0,:,:), qv)
qdiv = uv2dv_cfd(qu, qv, lat, lon, 1)
qdiv = qdiv*1.e+6
copy_VarCoords(t(:,0,:,:), qdiv)

qu_clm = month_to_season(qu, season)
qv_clm = month_to_season(qv, season)
qdiv_clm = month_to_season(qdiv, season)
printVarSummary(qdiv_clm)


rc = regCoef_n(year, t2m_clm, 0, 0)
copy_VarCoords(t2m_clm(0,:,:), rc)
rc_qu = regCoef_n(year, qu_clm, 0, 0)
rc_qv = regCoef_n(year, qv_clm, 0, 0)
rc_qdiv = regCoef_n(year, qdiv_clm, 0, 0)
copy_VarCoords(qu(0,:,:), rc_qu)
copy_VarCoords(qv(0,:,:), rc_qv)
copy_VarCoords(qdiv(0,:,:), rc_qdiv)

tval = onedtond(rc@tval, dimsizes(rc))
df   = onedtond(rc@nptxy, dimsizes(rc))
b    = tval
b    = 0.5
prob = betainc(df/(df+tval^2), df/2.0, b)
copy_VarCoords(rc, prob)
prob_qu = betainc(df/(df+tval^2), df/2.0, b)
copy_VarCoords(rc_qu, prob_qu)
prob_qv = betainc(df/(df+tval^2), df/2.0, b)
copy_VarCoords(rc_qv, prob_qv)
prob_qdiv = betainc(df/(df+tval^2), df/2.0, b)
copy_VarCoords(rc_qdiv, prob_qdiv)


rc@long_name = "regression coefficient"
prob@long_name = "probability"

trend = rc*10
copy_VarCoords(rc, trend)
printVarSummary(trend)

; rc_qu = rc_qu*10
; rc_qv = rc_qv*10
; rc_qdiv = rc_qdiv*10

sig = 0.05
rc_qu = mask(rc_qu, prob_qu.le.sig .and. prob_qv.le.sig, True)
rc_qv = mask(rc_qv, prob_qu.le.sig .and. prob_qv.le.sig, True)

;===============================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_surrounding/20210705/"
figname = "Fig5"
wks = gsn_open_wks("eps", figpath+figname)

plot = new(2, graphic)
plot_topo = plot
plotv = plot
plot_prob = new(2, graphic)
dum1 = plot
dum2 = plot
dum3 = plot
dum4 = plot

 res_topo            = True
 res_topo@cnFillOn   = False
 res_topo@cnLinesOn  = True
 res_topo@gsnFrame   = False
 res_topo@gsnDraw    = False
 res_topo@cnLineLabelsOn         = False
 res_topo@cnLineThicknessF       = 5
 res_topo@cnLineColor            = "forestgreen"
 res_topo@cnInfoLabelOn          = False
 res_topo@cnLevelSelectionMode   ="ExplicitLevels"
 res_topo@cnLevels               = (/2000/)
 res_topo@gsnLeftString = ""
 res_topo@gsnRightString = ""


;----------------------------------------------
; draw prob
;----------------------------------------------
 res3  = True
 res3@gsnDraw   = False
 res3@gsnFrame  = False

 ;res3@cnFillOn  = False 
 res3@cnLinesOn = False
 res3@cnLineLabelsOn = False
 res3@cnInfoLabelOn  = False
 ;res3@lbLabelBarOn   = False

 res3@gsnRightString = ""

 res3@cnLevelSelectionMode = "ManualLevels"
 res3@cnMinLevelValF = 0.00
 res3@cnMaxLevelValF = 0.10
 res3@cnLevelSpacingF = 0.01

 ; res3@cnFillPattern = 17 ;(/17/)
 res3@cnFillColor  =  "black"
 res3@gsnLeftString = ""
 res3@cnFillDotSizeF = 0.003
 res3@cnFillScaleF  = 1.2
 ;----------------------------------------------
 ;----------------------------------------------

res = True
res@gsnDraw = False
res@gsnFrame = False
res@gsnAddCyclic = True
res@gsnMaximize = True

res@mpOutlineOn = True
res@mpFillOn = True
res@mpGeophysicalLineThicknessF = 5
res@mpGeophysicalLineColor = "grey25" ;"red" 
; res@mpFillDrawOrder = "PostDraw"
; res@mpCenterLonF = 180
res@mpMinLatF = latstr
res@mpMaxLatF = latend
res@mpMinLonF = lonstr
res@mpMaxLonF = lonend
res@mpCenterLonF = (lonstr+lonend)/2


res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False
; res@cnFillMode = "CellFill"

res@gsnStringFontHeightF = 0.025
res@tmXBLabelFontHeightF = 0.02
res@tmYLLabelFontHeightF = 0.02

res@cnFillPalette = "MPL_RdBu"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -0.01
res@cnMaxLevelValF = 0.01
res@cnLevelSpacingF = 0.002
; res@lbLabelBarOn = False

res@pmLabelBarHeightF = 0.1
res@pmLabelBarWidthF = 0.7
res@lbLabelFontHeightF = 0.02
res@pmLabelBarOrthogonalPosF = 0.2

res@gsnRightString = ""


;-----------------------------------
gres = True
gres@gsLineColor = "blue"
gres@gsLineThicknessF = 5.0
gres@gsLineDashPattern = 16

 glon1 = (/lonend_north,lonend_north,lonstr_north,lonstr_north,lonend_north/)
 glat1 = (/latstr_north,latend_north,latend_north,latstr_north,latstr_north/) ; north
 glon2 = (/lonend_east,lonend_east,lonstr_east,lonstr_east,lonend_east/)
 glat2 = (/latstr_east,latend_east,latend_east,latstr_east,latstr_east/) ; east
 glon3 = (/lonend_south,lonend_south,lonstr_south,lonstr_south,lonend_south/)
 glat3 = (/latstr_south,latend_south,latend_south,latstr_south,latstr_south/) ; south
 glon4 = (/lonend_west,lonend_west,lonstr_west,lonstr_west,lonend_west/)
 glat4 = (/latstr_west,latend_west,latend_west,latstr_west,latstr_west/) ; west
;-----------------------------------

res@gsnLeftString = "(a) Column-integrated Cloud Ice and Liq Water"
res@gsnRightString = "10~S~-1~N~"
plot(0) = gsn_csm_contour_map(wks, trend, res)
plot_topo(0) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(0), plot_topo(0))

plot_prob(0) = gsn_csm_contour(wks,prob,res3)
plot_prob(0) = ShadeLtContour(plot_prob(0),sig,17)
overlay(plot(0), plot_prob(0))

dum1 = gsn_add_polyline(wks, plot(0), glon1, glat1, gres)
dum2 = gsn_add_polyline(wks, plot(0), glon2, glat2, gres)
dum3 = gsn_add_polyline(wks, plot(0), glon3, glat3, gres)
dum4 = gsn_add_polyline(wks, plot(0), glon4, glat4, gres)



;------------------------------------------

  vcres                         = True             ; plot mods desired
  vcres@gsnFrame                = False            ; so we can draw time stamp
  vcres@gsnDraw                 = False
  vcres@vcRefAnnoOrthogonalPosF = -0.9    ; -0.35          ; move ref vector up

  vcres@vcGlyphStyle            = "LineArrow"  ;"CurlyVector"    ; turn on curly vectors
  ;vcres@vcMinDistanceF          = 0.017
  ;vcres@vcGlyphStyle            = "LineArrow"
  vcres@vcMinDistanceF          = 0.03
  vcres@vcRefLengthF            = 0.045            ; define length of vec ref
  ; vcres@vcLineArrowColor        = "red"
  
  vcres@vcLineArrowThicknessF   = 5.0
  ;vecres@vcLineArrowHeadMinSizeF = 0.008
  ;vecres@vcLineArrowHeadMaxSizeF = 0.3

  vcres@vcRefAnnoString2On = False
  vcres@vcRefAnnoFontHeightF = 0.02
  vcres@gsnLeftString = ""
  vcres@gsnRightString = ""
  
  vcres@vcRefMagnitudeF         = 1.0 
;------------------------------------------
res@cnFillPalette = "CBR_drywet"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -1
res@cnMaxLevelValF = 1
res@cnLevelSpacingF = 0.2

res@gsnSpreadColors = True
res@gsnSpreadColorEnd   = 0
res@gsnSpreadColorStart = 10

res@gsnLeftString = "(b) WV Flux and Divergence"
res@gsnRightString = "10~S~-6~N"
plot(1) = gsn_csm_contour_map(wks, rc_qdiv, res)
plot_topo(1) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(1), plot_topo(1))


plotv(1) = gsn_csm_vector(wks, rc_qu, rc_qv, vcres)
overlay(plot(1), plotv(1))

plot_prob(1) = gsn_csm_contour(wks,prob_qdiv,res3)
plot_prob(1) = ShadeLtContour(plot_prob(1),sig,17)
overlay(plot(1), plot_prob(1))

dum11 = gsn_add_polyline(wks, plot(1), glon1, glat1, gres)
dum22 = gsn_add_polyline(wks, plot(1), glon2, glat2, gres)
dum33 = gsn_add_polyline(wks, plot(1), glon3, glat3, gres)
dum44 = gsn_add_polyline(wks, plot(1), glon4, glat4, gres)

pres = True
pres@gsnMaximize = True
gsn_panel(wks, plot, (/2,1/), pres)
end