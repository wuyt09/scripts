
begin

yrStrt = 1979
yrLast = 2016
year = ispan(yrStrt, yrLast, 1)
nyear = dimsizes(year)
season = "JJA"

latstr_tp = 25 ;28
latend_tp = 40 ;38
lonstr_tp = 75 ;80
lonend_tp = 105

;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo ;(::-1,:)
printVarSummary(topo)

;=============================================================================================================
diri1 = "/home/yangsong3/data-observation/ERA5-monthly/surface/"
fili1 = "surface_sensible_heat_flux_mon_1x1_1979_2019.nc"

diri2 = "/home/yangsong3/data-observation/wuyt/ERAI/monthly/1x1/"
fili2 = "sshf.mon.mean.7916.r1.nc"


f1 = addfile(diri1+fili1, "r")
Time = f1->time
YYYY = cd_calendar(Time, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

f2 = addfile(diri2+fili2, "r")

SSHF1 = short2flt(f1->sshf(iYYYY,:,:))
SSHF1 = SSHF1/86400.
SSHF2 = short2flt(f2->sshf(iYYYY,:,:))
; SSHF2 = SSHF1/86400.
; SSHF = SSHF*(-1)
printVarSummary(SSHF1)

nlat = dimsizes(SSHF1&latitude)
nlon = dimsizes(SSHF1&longitude)
var = (/"SH1","SH2"/)
x = new((/2,nyear*12,nlat,nlon/), float)
x!0 = "var"
x&var = var
x(0,:,:,:) = SSHF1
x(1,:,:,:) = SSHF2
copy_VarCoords(SSHF1, x(0,:,:,:))

x = mask(x, topo.ge.1500, True)
printVarSummary(x)

x_jja = new((/2,nyear,nlat,nlon/), float)
x_jja!0 = "var"
x_jja&var = var
do i = 0, 1
    x_jja(i,:,:,:) = month_to_season(x(i,:,:,:), season)
end do
x_jja&time = year
printVarSummary(x_jja)

x_region = wgt_areaave_Wrap(x_jja(:,:,{latstr_tp:latend_tp},{lonstr_tp:lonend_tp}), 1., 1., 0)
printVarSummary(x_region)

; x_region_runave = runave_n_Wrap(x_region, 10, 0, 1)
; printVarSummary(x_region_runave)

; lines = new((/2,6,nyear/), float)
; lines(0,:,:) = x_region
; lines(1,:,:) = x_region_runave
; copy_VarCoords(x_region, lines(0,:,:))
; printVarSummary(lines)

; corr = new(6, float)
; corr(0) = 1
; do i = 0, 4
;     corr(i+1) = escorc_n(lines(1,0,:), lines(1,i+1,:), 0, 0)
; end do
;=======================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_surrounding/20210520/"
figname = "check_SH_ERAI_ERA5_TP"

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

; do i = 0, 5
    ; res@gsnLeftString = Left(i)+" "+var(i)
    ; res@gsnCenterString = sprintf("%4.2f", corr(i))
;     plot(i) = gsn_csm_xy(wks, year ,x_region(:,:), res)
; end do

plot(0) = gsn_csm_xy(wks, year, x_region, res)
draw(plot(0))

; pres = True
; pres@gsnMaximize = True
; pres@gsnPanelMainString = "Surface energy bugdet: TP region"
; gsn_panel(wks, plot, (/2,3/), pres)
end