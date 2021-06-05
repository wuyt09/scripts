
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

latstr_east = 30
latend_east = 45
lonstr_east = 105
lonend_east = 140

latstr_west = 30
latend_west = 55
lonstr_west = -10
lonend_west = 60

latstr_south = 10
latend_south = 30 ;25
lonstr_south = 70 ;50
lonend_south = 90 ;100

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
topo = ftopo->topo(:,:)
printVarSummary(topo)

fa      = addfile ("/home/yangsong3/wuyt/data/lsm_1279l4_0.1x0.1.grb_v4_unpack.nc","r")
; lsdata = doubletofloat(a->lsm(0,:,:))
lsdata1 = fa->lsm(0,:,:)
lsdata = tofloat(lsdata1)
copy_VarCoords(lsdata1, lsdata)
printVarSummary(lsdata)
;=============================================================================================================
diri = "/home/yangsong3/data-observation/ERA5-monthly/surface/"
fili1 = "surface_solar_downwards_mon_1x1_1979_2019.nc"
fili2 = "surface_thermal_downwards_mon_1x1_1979_2019.nc"
fili3 = "surface_net_solar_mon_1x1_1979_2019.nc"
fili4 = "surface_net_thermal_mon_1x1_1979_2019.nc"
fili5 = "surface_latent_heat_flux_mon_1x1_1979_2019.nc"
fili6 = "surface_sensible_heat_flux_mon_1x1_1979_2019.nc"

f1 = addfile(diri+fili1, "r")
Time = f1->time
YYYY = cd_calendar(Time, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

f2 = addfile(diri+fili2, "r")
f3 = addfile(diri+fili3, "r")
f4 = addfile(diri+fili4, "r")
f5 = addfile(diri+fili5, "r")
f6 = addfile(diri+fili6, "r")


SSRD = short2flt(f1->ssrd(iYYYY,:,:))
SSRD = SSRD/86400.
STRD = short2flt(f2->strd(iYYYY,:,:))
STRD = STRD/86400.
SSR  = short2flt(f3->ssr(iYYYY,:,:))
SSR  = SSR/86400.
STR  = short2flt(f4->str(iYYYY,:,:))
STR  = STR/86400.
SLHF = short2flt(f5->slhf(iYYYY,:,:))
SLHF = SLHF/86400.
; SLHF = SLHF*(-1)
SSHF = short2flt(f6->sshf(iYYYY,:,:))
SSHF = SSHF/86400.
; SSHF = SSHF*(-1)


G0 = STR+SSR+SLHF+SSHF
copy_VarCoords(STR,G0)
t2m_clm = month_to_season(G0, season)

lsm = linint2_Wrap(lsdata&longitude, lsdata&latitude, lsdata, True, SSR&longitude, SSR&latitude, 0)
printVarSummary(lsm)
delete(lsdata)
lsm1 = conform_dims(dimsizes(t2m_clm), lsm, (/1,2/))
x_clm = where(lsm1.ge.0.8, t2m_clm, SSR@_FillValue)
copy_VarCoords(t2m_clm, x_clm)
printMinMax(x_clm, False)
printMinMax(t2m_clm, False)
; exit
lines_west = wgt_areaave_Wrap(x_clm(:,{latstr_west:latend_west},{lonstr_west:lonend_west}), 1., 1., 0)
lines_east = wgt_areaave_Wrap(x_clm(:,{latstr_east:latend_east},{lonstr_east:lonend_east}), 1., 1., 0)
lines_north = wgt_areaave_Wrap(x_clm(:,{latstr_north:latend_north},{lonstr_north:lonend_north}), 1., 1., 0)
lines_south = wgt_areaave_Wrap(x_clm(:,{latstr_south:latend_south},{lonstr_south:lonend_south}), 1., 1., 0)

x_mask = mask(t2m_clm, topo.ge.1500, True)
copy_VarCoords(t2m_clm, x_mask)
lines_tp = wgt_areaave_Wrap(x_mask(:,{latstr_tp:latend_tp},{lonstr_tp:lonend_tp}), 1., 1., 0)

lines = new((/5,nyear/), float)
lines(0,:) = lines_west - lines_west(0)
lines(1,:) = lines_east - lines_east(0)
lines(2,:) = lines_south - lines_south(0)
lines(3,:) = lines_north - lines_north(0)
lines(4,:) = lines_tp - lines_tp(0)
copy_VarCoords(t2m_clm(:,0,0), lines(0,:))
lines!0 = "region"
lines&region = region
printVarSummary(lines)

print(sprintf("%4.2f", lines_west)+"    "+sprintf("%4.2f", lines_east)+"    "+sprintf("%4.2f", lines_tp))
print(sprintf("%4.2f", lines_north)+"   "+sprintf("%4.2f", lines_south))
; exit()
;=======================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_surrounding/20210520/"
figname = "check_ERA5_surf_energy_budget"

Left = (/"(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)","(j)","(k)","(l)"/)

wks = gsn_open_wks("pdf", figpath+figname)
gsn_define_colormap(wks, "temp_diff_18lev")
plot = new(1, graphic)


res = True
res@gsnDraw  = False
res@gsnFrame = False
res@gsnMaximize = True

res@vpWidthF = 0.8
res@vpHeightF = 0.4

res@trXMinF = yrStrt
res@trXMaxF = yrLast
; res@tmXBMode = "Explicit"
; res@tmXBValues = year
; res@tmXBLabels = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
res@tmXTOn = False
; res@tiXAxisString = "Month"

res@gsnYRefLine = 0
res@gsnYRefLineColor = "gray25"
res@gsnYRefLineDashPattern = 16
res@gsnYRefLineThicknessF = 8


; res@tmYLMode = "Explicit"
; res@tmYLValues = ispan(-20, 20, 5)
; res@tmYLLabels = ispan(-20, 20, 5)
; res@tmYLPrecision = 4
res@tmYROn = False
res@tiYAxisString = ""
res@gsnStringFontHeightF = 0.025

;Colors = (/"deeppink2", "blue", "purple","green4","orange2", "black" /);
Colors = (/"forestgreen","red","orange2","blue","black","forestgreen","red","orange2","blue","black"/)
; Patterns = (/0,16,16,0,0/)
Patterns = (/0,0,0,0,0,16,16,16,16,16/)

res@xyLineThicknessF = 10
res@xyLineColors = Colors
res@xyDashPatterns = Patterns

; res@gsnRightString = "trend: "+sprintf("%4.2f", trend(j))+" K/decade"
res@gsnLeftString = season+" G0 trend ("+yrStrt+"-"+yrLast+")"
plot(0) = gsn_csm_xy(wks, year, lines, res)

; res = True
; res@gsnDraw = False
; res@gsnFrame = False
; res@gsnMaximize = True

; res@mpOutlineOn = True
; res@mpFillOn = True
; res@mpGeophysicalLineThicknessF = 5
; res@mpGeophysicalLineColor = "grey25" ;"red" 
; ; res@mpFillDrawOrder = "PostDraw"
; ; res@mpCenterLonF = 180
; res@mpMinLatF = latstr
; res@mpMaxLatF = latend
; res@mpMinLonF = lonstr
; res@mpMaxLonF = lonend


; res@cnFillOn = True
; res@cnLinesOn = False
; res@cnLineLabelsOn = False
; res@cnInfoLabelOn = False
; ; res@cnFillMode = "CellFill"

; res@gsnStringFontHeightF = 0.025
; res@tmXBLabelFontHeightF = 0.02
; res@tmYLLabelFontHeightF = 0.02

; res@cnFillPalette = "CBR_coldhot"
; res@cnLevelSelectionMode = "ManualLevels"
; res@cnMinLevelValF = -1
; res@cnMaxLevelValF = 1
; res@cnLevelSpacingF = 0.2

; plot(0) = gsn_csm_contour_map(wks, lsm, res)
draw(plot(0))
end