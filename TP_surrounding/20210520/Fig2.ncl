
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

region = (/"TP","South","West","East","North"/)
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
fili1 = "surface_solar_downwards_mon_1x1_1979_2019.nc"
fili2 = "surface_thermal_downwards_mon_1x1_1979_2019.nc"
fili3 = "surface_net_solar_mon_1x1_1979_2019.nc"
fili4 = "surface_net_thermal_mon_1x1_1979_2019.nc"

f1 = addfile(diri+fili1, "r")
Time = f1->time
YYYY = cd_calendar(Time, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

f2 = addfile(diri+fili2, "r")
f3 = addfile(diri+fili3, "r")
f4 = addfile(diri+fili4, "r")

SSRD = short2flt(f1->ssrd(iYYYY,:,:))
SSRD = SSRD/86400.
STRD = short2flt(f2->strd(iYYYY,:,:))
STRD = STRD/86400.
SSR  = short2flt(f3->ssr(iYYYY,:,:))
SSR  = SSR/86400.
STR  = short2flt(f4->str(iYYYY,:,:))
STR  = STR/86400.

SSRU = SSRD - SSR
SSRU = SSRU*(-1)
copy_VarCoords(SSRD, SSRU)
STRU = STRD - STR
copy_VarCoords(STRD, STRU)
printVarSummary(SSRU)


ssrd = month_to_season(SSRD, season)
ssru = month_to_season(SSRU, season)
strd = month_to_season(STRD, season)
stru = month_to_season(STRU, season)

lat = SSRD&latitude
lon = SSRD&longitude
lsm = linint2_Wrap(lsdata&longitude, lsdata&latitude, lsdata, True, lon, lat, 0)
printVarSummary(lsm)
delete(lsdata)
lsm1 = conform_dims(dimsizes(ssrd), lsm, (/1,2/))

ssrd = where(lsm1.ge.0.8, ssrd, ssrd@_FillValue)
ssru = where(lsm1.ge.0.8, ssru, ssru@_FillValue)
strd = where(lsm1.ge.0.8, strd, strd@_FillValue)
stru = where(lsm1.ge.0.8, stru, stru@_FillValue)

stru_west = wgt_areaave_Wrap(stru(:,{latstr_west:latend_west},{lonstr_west:lonend_west}), 1., 1., 0)
stru_east = wgt_areaave_Wrap(stru(:,{latstr_east:latend_east},{lonstr_east:lonend_east}), 1., 1., 0)
stru_south = wgt_areaave_Wrap(stru(:,{latstr_south:latend_south},{lonstr_south:lonend_south}), 1., 1., 0)
stru_north = wgt_areaave_Wrap(stru(:,{latstr_north:latend_north},{lonstr_north:lonend_north}), 1., 1., 0)

strd_west = wgt_areaave_Wrap(strd(:,{latstr_west:latend_west},{lonstr_west:lonend_west}), 1., 1., 0)
strd_east = wgt_areaave_Wrap(strd(:,{latstr_east:latend_east},{lonstr_east:lonend_east}), 1., 1., 0)
strd_south = wgt_areaave_Wrap(strd(:,{latstr_south:latend_south},{lonstr_south:lonend_south}), 1., 1., 0)
strd_north = wgt_areaave_Wrap(strd(:,{latstr_north:latend_north},{lonstr_north:lonend_north}), 1., 1., 0)

ssrd_west = wgt_areaave_Wrap(ssrd(:,{latstr_west:latend_west},{lonstr_west:lonend_west}), 1., 1., 0)
ssrd_east = wgt_areaave_Wrap(ssrd(:,{latstr_east:latend_east},{lonstr_east:lonend_east}), 1., 1., 0)
ssru_south = wgt_areaave_Wrap(ssru(:,{latstr_south:latend_south},{lonstr_south:lonend_south}), 1., 1., 0)
ssrd_north = wgt_areaave_Wrap(ssrd(:,{latstr_north:latend_north},{lonstr_north:lonend_north}), 1., 1., 0)

stru = mask(stru, topo.ge.1500, True)
strd = mask(strd, topo.ge.1500, True)
ssru = mask(ssru, topo.ge.1500, True)
stru_tp = wgt_areaave_Wrap(stru(:,{latstr_tp:latend_tp},{lonstr_tp:lonend_tp}), 1., 1., 0)
strd_tp = wgt_areaave_Wrap(strd(:,{latstr_tp:latend_tp},{lonstr_tp:lonend_tp}), 1., 1., 0)
ssru_tp = wgt_areaave_Wrap(ssru(:,{latstr_tp:latend_tp},{lonstr_tp:lonend_tp}), 1., 1., 0)
printVarSummary(ssru_tp)

lines = new((/3,5,nyear/), float)
lines(0,0,:) = stru_tp
lines(0,1,:) = stru_south
lines(0,2,:) = stru_west
lines(0,3,:) = stru_north
lines(0,4,:) = stru_east

lines(1,0,:) = strd_tp
lines(1,1,:) = strd_south
lines(1,2,:) = strd_west
lines(1,3,:) = strd_north
lines(1,4,:) = strd_east

lines(2,0,:) = ssru_tp
lines(2,1,:) = ssru_south
lines(2,2,:) = ssrd_west
lines(2,3,:) = ssrd_north
lines(2,4,:) = ssrd_east

; lines(3,:,:) = lines(2,:,:) + lines(1,:,:)

lines_runave = runave_n(lines, 10, 0, 2)
printVarSummary(lines_runave)

line = new((/2,3,5,nyear/), float)
line(0,:,:,:) = lines
line(1,:,:,:) = lines_runave
printVarSummary(line)

tt = line(0,0,:,:)
tt(0,:) = line(1,1,0,:) - line(1,2,0,:)
tt(1,:) = line(1,1,1,:) - line(1,2,1,:)
tt(2,:) = line(1,1,2,:) + line(1,2,2,:)
tt(3,:) = line(1,1,3,:) + line(1,2,3,:)
tt(4,:) = line(1,1,4,:) + line(1,2,4,:)
corr = escorc_n(line(1,0,:,:), tt, 1, 1)
print(sprintf("%4.2f", corr))
;===============================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_surrounding/20210520/"
figname = "Fig2"
wks = gsn_open_wks("pdf", figpath+figname)

Left = (/"(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)","(j)","(k)","(l)","(m)","(n)","(o)"/)
Right = (/"upR","upR","upR","upR","upR","downR","downR","downR","downR","downR","upS","upS","downS","downS","downS"/)

plot = new(15, graphic)

res = True
res@gsnDraw  = False
res@gsnFrame = False
; res@gsnMaximize = True

res@vpWidthF = 0.6
res@vpHeightF = 0.4

res@trXMinF = yrStrt
res@trXMaxF = yrLast
res@tmXTOn = False

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

Colors = (/"black","purple","green4"/)
Patterns = (/0,16,16/)

res@xyLineThicknessF = 10
res@xyLineColors = Colors
res@xyDashPatterns = Patterns

; ; res@gsnRightString = "trend: "+sprintf("%4.2f", trend(j))+" K/decade"
; res@gsnLeftString = season+" T2m trend ("+yrStrt+"-"+yrLast+")"
; plot(0) = gsn_csm_xy(wks, year, pp, res)

do i = 0, 2
    do j = 0, 4
        res@gsnLeftString = Left(j+i*5)+" "+Right(j+i*5)
        res@gsnRightString = region(j)
        if (i.eq.0) then
            res@gsnCenterString = sprintf("%4.2f", corr(j))
        end if
        plot(j+i*5) = gsn_csm_xy(wks, year, line(:,i,j,:), res)
        res@gsnCenterString = ""
    end do
end do

pres = True
pres@gsnMaximize = True
gsn_panel(wks, plot, (/3,5/), pres)

end