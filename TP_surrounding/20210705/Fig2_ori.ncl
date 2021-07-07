
begin


yrStrt = 1979
yrLast = 2019
year = ispan(yrStrt, yrLast, 1)
nyear = dimsizes(year)
season = "JJA"

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

region = (/"West","North","East","TP","South"/)
;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo ;(::-1,:)
printVarSummary(topo)


fa      = addfile ("/home/yangsong3/wuyt/data/lsm_1279l4_0.1x0.1.grb_v4_unpack.nc","r")
; lsdata = doubletofloat(a->lsm(0,:,:))
lsdata = fa->lsm(0,:,:)
printVarSummary(lsdata)
;=============================================================================================================
diri = "/home/yangsong3/data-observation/ERA5-monthly/surface/"
fili1 = "surface_thermal_downwards_mon_1x1_1979_2019.nc"
fili2 = "surface_net_thermal_mon_1x1_1979_2019.nc"
fili3 = "t2m_mon_1x1_1979_2019.nc"

f1 = addfile(diri+fili1, "r")
Time = f1->time
YYYY = cd_calendar(Time, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

f2 = addfile(diri+fili2, "r")
f3 = addfile(diri+fili3, "r")

STRD = short2flt(f1->strd(iYYYY,:,:))
STRD = STRD/86400.
STR  = short2flt(f2->str(iYYYY,:,:))
STR  = STR/86400.
T2M  = short2flt(f3->t2m(iYYYY,:,:))

STRU = STRD - STR
copy_VarCoords(STRD, STRU)
printVarSummary(STRU)

stru = month_to_season(STRU, season)
t2m = month_to_season(T2M, season)

lat = STRD&latitude
lon = STRD&longitude
lsm = linint2_Wrap(lsdata&longitude, lsdata&latitude, lsdata, True, lon, lat, 0)
printVarSummary(lsm)
delete(lsdata)
lsm1 = conform_dims(dimsizes(stru), lsm, (/1,2/))

t2m  = where(lsm1.ge.0.8, t2m, t2m@_FillValue)
stru = where(lsm1.ge.0.8, stru, stru@_FillValue)

stru_west = wgt_areaave_Wrap(stru(:,{latstr_west:latend_west},{lonstr_west:lonend_west}), 1., 1., 0)
stru_east = wgt_areaave_Wrap(stru(:,{latstr_east:latend_east},{lonstr_east:lonend_east}), 1., 1., 0)
stru_south = wgt_areaave_Wrap(stru(:,{latstr_south:latend_south},{lonstr_south:lonend_south}), 1., 1., 0)
stru_north = wgt_areaave_Wrap(stru(:,{latstr_north:latend_north},{lonstr_north:lonend_north}), 1., 1., 0)

t2m_west = wgt_areaave_Wrap(t2m(:,{latstr_west:latend_west},{lonstr_west:lonend_west}), 1., 1., 0)
t2m_east = wgt_areaave_Wrap(t2m(:,{latstr_east:latend_east},{lonstr_east:lonend_east}), 1., 1., 0)
t2m_south = wgt_areaave_Wrap(t2m(:,{latstr_south:latend_south},{lonstr_south:lonend_south}), 1., 1., 0)
t2m_north = wgt_areaave_Wrap(t2m(:,{latstr_north:latend_north},{lonstr_north:lonend_north}), 1., 1., 0)

stru = mask(stru, topo.ge.1500, True)
t2m  = mask(t2m, topo.ge.1500, True)
stru_tp = wgt_areaave_Wrap(stru(:,{latstr_tp:latend_tp},{lonstr_tp:lonend_tp}), 1., 1., 0)
t2m_tp  = wgt_areaave_Wrap(t2m(:,{latstr_tp:latend_tp},{lonstr_tp:lonend_tp}), 1., 1., 0)
printVarSummary(stru_tp)

lines = new((/3,5,nyear/), float)
lines(2,0,:) = t2m_west ;- t2m_west(0)
lines(2,1,:) = t2m_north ;- t2m_north(0)
lines(2,2,:) = t2m_east ;- t2m_east(0
lines(2,3,:) = t2m_tp ;- t2m_tp(0)
lines(2,4,:) = t2m_south ;- t2m_south(0))

lines(0,0,:) = stru_west ;- stru_west(0)
lines(0,1,:) = stru_north ;- stru_north(0)
lines(0,2,:) = stru_east ;- stru_east(0)
lines(0,3,:) = stru_tp ;- stru_tp(0)
lines(0,4,:) = stru_south ;- stru_south(0)

lines_runave = runave_n(lines(0,:,:), 10, 0, 1)
lines(1,:,:) = (/lines_runave/)
; printVarSummary(lines_runave)

;===============================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_surrounding/20210705/"
figname = "Fig2"
wks = gsn_open_wks("eps", figpath+figname)

Left = (/"(a)","(b)","(c)","(d)","(e)"/)

plot = new(6, graphic)

res = True
res@gsnDraw  = False
res@gsnFrame = False
; res@gsnMaximize = True

res@vpWidthF = 0.8
res@vpHeightF = 0.4

res@trXMinF = yrStrt
res@trXMaxF = yrLast
; res@tmXTOn = False

res@gsnYRefLine = 0
res@gsnYRefLineColor = "gray25"
res@gsnYRefLineDashPattern = 16
res@gsnYRefLineThicknessF = 8


; res@tmYLMode = "Explicit"
; res@tmYLValues = ispan(-20, 20, 5)
; res@tmYLLabels = ispan(-20, 20, 5)
; res@tmYLPrecision = 4
res@tiYAxisString = ""
res@gsnStringFontHeightF = 0.035
res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025
res@tmYLPrecision = 4

Colors = (/"black","purple"/)
Patterns = (/0,16/)

res@xyLineThicknessF = 5
res@xyLineColors = Colors
res@xyDashPatterns = Patterns

; ; res@gsnRightString = "trend: "+sprintf("%4.2f", trend(j))+" K/decade"
; res@gsnLeftString = season+" T2m trend ("+yrStrt+"-"+yrLast+")"
; plot(0) = gsn_csm_xy(wks, year, pp, res)

res2 = True
res2@xyLineColor = "firebrick2"
res2@xyDashPattern = 0
res2@xyLineThicknessF = 5
res2@tiYAxisString = ""
res2@tmYRPrecision = 4
res2@tmYRLabelFontHeightF = 0.025


res@gsnLeftString = "(a) West region ("+lonstr_west+"-"+lonend_west+"E, "+latstr_west+"-"+latend_west+"N)"
res@trYMinF = 420
res@trYMaxF = 456
res2@trYMinF = 294
res2@trYMaxF = 300
plot(0) = gsn_csm_xy2(wks, year, lines(0:1,0,:), lines(2,0,:), res, res2)

res@gsnLeftString = "(b) North region ("+lonstr_north+"-"+lonend_north+"E, "+latstr_north+"-"+latend_north+"N)"
res@trYMinF = 380
res@trYMaxF = 410
res2@trYMinF = 287
res2@trYMaxF = 295
plot(1) = gsn_csm_xy2(wks, year, lines(0:1,1,:), lines(2,1,:), res, res2)

res@gsnLeftString = "(c) East region ("+lonstr_east+"-"+lonend_east+"E, "+latstr_east+"-"+latend_east+"N)"
res@trYMinF = 415
res@trYMaxF = 445
res2@trYMinF = 294
res2@trYMaxF = 300
plot(2) = gsn_csm_xy2(wks, year, lines(0:1,2,:), lines(2,2,:), res, res2)

res@gsnLeftString = "(d) TP region ("+lonstr_tp+"-"+lonend_tp+"E, "+latstr_tp+"-"+latend_tp+"N)"
res@trYMinF = 360
res@trYMaxF = 380
res2@trYMinF = 281
res2@trYMaxF = 286
plot(3) = gsn_csm_xy2(wks, year, lines(0:1,3,:), lines(2,3,:), res, res2)

res@gsnLeftString =  "(e) South region ("+lonstr_south+"-"+lonend_south+"E, "+latstr_south+"-"+latend_south+"N)"
res@trYMinF = 445
res@trYMaxF = 485
res2@trYMinF = 298
res2@trYMaxF = 308
plot(4) = gsn_csm_xy2(wks, year, lines(0:1,4,:), lines(2,4,:), res, res2)

;===========================================================
; legend
;===========================================================
lb_colors = (/"black","purple","firebrick2"/)
lb_pattern = (/0,16,0/)
xpos = (/0.72,0.76/)
ypos = (/0.45,0.42,0.39/)
gsres = True
gsres@gsLineThicknessF = 5

do i = 0, 2
    gsres@gsLineColor = lb_colors(i)
    gsres@gsLineDashPattern = lb_pattern(i)
    gsn_polyline_ndc(wks, xpos, (/ypos(i),ypos(i)/), gsres)
end do

Text = (/"Surface upward LW","10-yr running average upward LW","T2m"/)
txres = True
txres@txFontHeightF = 0.01
txres@txJust             = "CenterLeft"
tx_xpos = 0.77
tx_ypos = ypos
do i = 0, 2
    gsn_text_ndc(wks, Text(i), tx_xpos, tx_ypos(i), txres)
end do

pres = True
pres@gsnMaximize = True
gsn_panel(wks, plot, (/2,3/), pres)

end