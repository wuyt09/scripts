
begin

yrStrt = 1979
yrLast = 2019
year = ispan(yrStrt, yrLast, 1)
nyear = dimsizes(year)
season = "JJA"

latstr_west = 30
latend_west = 55
lonstr_west = -10
lonend_west = 60

;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo ;(::-1,:)
printVarSummary(topo)

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

SSRU = SSRD - SSR
SSRU = SSRU*(-1)
copy_VarCoords(SSRD, SSRU)
STRU = STRD - STR
copy_VarCoords(STRD, STRU)
printVarSummary(SSRU)


nlat = dimsizes(SSR&latitude)
nlon = dimsizes(SSR&longitude)
var = (/"upR","downR","upS","downS","SH","LH"/)
x = new((/6,nyear*12,nlat,nlon/), float)
x!0 = "var"
x&var = var
x(0,:,:,:) = STRU
x(1,:,:,:) = STRD
x(2,:,:,:) = SSRU
x(3,:,:,:) = SSRD
x(4,:,:,:) = SSHF
x(5,:,:,:) = SLHF
copy_VarCoords(SSR, x(0,:,:,:))

x_jja = new((/6,nyear,nlat,nlon/), float)
x_jja!0 = "var"
x_jja&var = var
do i = 0, 5
    x_jja(i,:,:,:) = month_to_season(x(i,:,:,:), season)
end do
x_jja&time = year
printVarSummary(x_jja)

x_region = wgt_areaave_Wrap(x_jja(:,:,{latstr_west:latend_west},{lonstr_west:lonend_west}), 1., 1., 0)
printVarSummary(x_region)

x_region_runave = runave_n_Wrap(x_region, 10, 0, 1)
printVarSummary(x_region_runave)

lines = new((/2,6,nyear/), float)
lines(0,:,:) = x_region
lines(1,:,:) = x_region_runave
copy_VarCoords(x_region, lines(0,:,:))
printVarSummary(lines)

corr = new(6, float)
corr(0) = 1
do i = 0, 4
    corr(i+1) = escorc_n(lines(1,0,:), lines(1,i+1,:), 0, 0)
end do
;=======================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_surrounding/20210520/"
figname = "ERA5_surf_energy_budget_WEST"

Left = (/"(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)","(j)","(k)","(l)"/)

wks = gsn_open_wks("pdf", figpath+figname)
gsn_define_colormap(wks, "temp_diff_18lev")
plot = new(6, graphic)


res = True
res@gsnDraw  = False
res@gsnFrame = False

res@vpHeightF = 0.4
res@vpWidthF  = 0.7

res@tiYAxisString = ""
res@gsnStringFontHeightF = 0.03
res@tmXBLabelFontHeightF = 0.02
res@tmYLLabelFontHeightF = 0.02
res@gsnRightString = "W/m~S~2~N~"

res@trXMinF = yrStrt
res@trXMaxF = yrLast

Colors = (/"black","purple"/)
Patterns = (/0,16/)
res@xyLineThicknessF = 10
res@xyLineColors = Colors
res@xyDashPatterns = Patterns

do i = 0, 5
    res@gsnLeftString = Left(i)+" "+var(i)
    res@gsnCenterString = sprintf("%4.2f", corr(i))
    plot(i) = gsn_csm_xy(wks, year ,lines(:,i,:), res)
end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelMainString = "Surface energy bugdet: West region"
gsn_panel(wks, plot, (/2,3/), pres)
end