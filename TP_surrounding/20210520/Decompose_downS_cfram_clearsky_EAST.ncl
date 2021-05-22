
begin

yrStrt = 1979
yrLast = 2019
year = ispan(yrStrt, yrLast, 1)
nyear = dimsizes(year)
season = "JJA"

latstr = 35
latend = 45
lonstr = 105
lonend = 140

;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo ;(::-1,:)
printVarSummary(topo)

;=============================================================================================================
diri = "/home/yangsong3/data-observation/ERA5-monthly/surface/"
fili1 = "surface_solar_downwards_mon_1x1_1979_2019.nc"
fili2 = "surface_solar_downwards_clear_sky_mon_1x1_1979_2019.nc"
; fili3 = "surface_net_thermal_mon_1x1_1979_2019.nc"
; fili4 = "surface_net_thermal_clear_sky_mon_1x1_1979_2019.nc"

f1 = addfile(diri+fili1, "r")
Time = f1->time
YYYY = cd_calendar(Time, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

f2 = addfile(diri+fili2, "r")
; f3 = addfile(diri+fili3, "r")
; f4 = addfile(diri+fili4, "r")

STRD = short2flt(f1->ssrd(iYYYY,:,:))
STRD = STRD/86400.
STRDC = short2flt(f2->ssrdc(iYYYY,:,:))
STRDC = STRDC/86400.
; STR  = short2flt(f3->str(iYYYY,:,:))
; STR  = STR/86400.
; STRC  = short2flt(f4->strc(iYYYY,:,:))
; STRC  = STRC/86400.

STRD_CRE = STRD - STRDC
copy_VarCoords(STRD, STRD_CRE)
printVarSummary(STRD_CRE)

strd = month_to_season(STRD, season)
strd&time = year
strd_cre = month_to_season(STRD_CRE, season)
strd_cre&time = year
printVarSummary(strd_cre)

obs = runave_n_Wrap(strd, 10, 0, 0)
printVarSummary(obs)
obs_region = wgt_areaave_Wrap(obs(:,{latstr:latend},{lonstr:lonend}), 1., 1., 0)
printVarSummary(obs_region)

cloud = runave_n_Wrap(strd_cre, 10, 0, 0)
printVarSummary(cloud)
cloud_region = wgt_areaave_Wrap(cloud(:,{latstr:latend},{lonstr:lonend}), 1., 1., 0)
printVarSummary(cloud_region)
; print(cloud_region)

;================================================================================
nn = 31
nx = 4
lat = cloud&latitude
lon = cloud&longitude
nlat = dimsizes(lat)
nlon = dimsizes(lon)

diri2 = "/home/yangsong3/wuyt/sysu/scripts/TP_surrounding/cfram20210513_clearsky/"
nvname = (/"rad_conv","lwu","lwd","swu","swd","lwu_toa","lwd_toa","swu_toa","swd_toa"/)  ; ,"error"
nvars = dimsizes(nvname)


UNDEF = -999.
x = new((/nn,6,nvars,nlat,nlon/),float,UNDEF)
x!0 = "time"
x!1 = "process"
x!2 = "var"
x!3 = "lat"
x!4 = "lon"
x&var = nvname
x&lat = lat
x&lon = lon
x&lat@units = "degrees_north"
x&lon@units = "degrees_east"
x@long_name = "forcing"
x@units     = "W/m2"


do nnn = 0, nn-1
    n = nnn+1
    fili1 = diri2+"baseline_radsfc_ranc_1.grd"
    fili2 = diri2+"solar_radsfc_ranc_"+n+".grd"
    fili3 = diri2+"o3_radsfc_ranc_"+n+".grd"
    fili4 = diri2+"wv_radsfc_ranc_"+n+".grd"
    fili5 = diri2+"albedo_radsfc_ranc_"+n+".grd"
    fili6 = diri2+"warm_radsfc_ranc_"+n+".grd"

    do nv = 0, nvars-1
        x(nnn,0,nv,:,:) = fbindirread(fili1,nv,(/nlat,nlon/),"float")
        x(nnn,1,nv,:,:) = fbindirread(fili2,nv,(/nlat,nlon/),"float")
        x(nnn,2,nv,:,:) = fbindirread(fili3,nv,(/nlat,nlon/),"float")
        x(nnn,3,nv,:,:) = fbindirread(fili4,nv,(/nlat,nlon/),"float")
        x(nnn,4,nv,:,:) = fbindirread(fili5,nv,(/nlat,nlon/),"float")
        x(nnn,5,nv,:,:) = fbindirread(fili6,nv,(/nlat,nlon/),"float")
    end do
end do

printVarSummary(x)


x_new = x(:,:,nx,:,:)

; cfram_wb = x_new(:,5,:,:) - x_new(:,0,:,:)
error = x_new(:,5,:,:) - dim_sum_n_Wrap(x_new(:,1:4,:,:), 1) + 3*x_new(:,0,:,:)
copy_VarCoords(x_new(:,0,:,:), error)
; error_TP = error(29,{latstr:latend},{lonstr:lonend}) ;wgt_areaave_Wrap(error(:,{latstr:latend},{lonstr:lonend}), 1, 1, 0)
; print(error_TP)
; exit

x_new(:,5,:,:) = where(abs(error).ge.5, x_new(:,5,:,:) - error, x_new(:,5,:,:))

; error = x_new(:,5,:,:) - dim_sum_n_Wrap(x_new(:,1:4,:,:), 1) + 3*x_new(:,0,:,:)
; copy_VarCoords(x_new(:,0,:,:), error)
; error_TP = error(29,{latstr:latend},{lonstr:lonend}) ;wgt_areaave_Wrap(error(:,{latstr:latend},{lonstr:lonend}), 1, 1, 0)
; print(error_TP)
; exit

x_region = wgt_areaave_Wrap(x_new(:,:,{latstr:latend},{lonstr:lonend}), 1., 1., 0)
printVarSummary(x_region)

yrs = 1984
yre = 2014
year_new = ispan(yrs, yre, 1)
process = (/"Solar","O3","WV","Albedo","CRE","Sum of a-e"/)
np = dimsizes(process)
lines = new((/nn,np/), float)
lines!0 = "time"
lines!1 = "process"
lines&time = year_new
lines&process = process

do i = 0, 3
    lines(:,i) = x_region(:,i+1) - x_region(:,0)
end do
lines(:,4) = cloud_region(5:35)
lines(:,4) = lines(:,4) - cloud_region(4)

cfram_total = dim_sum_n_Wrap(lines(:,0:4), 1)
; lines(:,5) = x_region(:,5) - x_region(:,0) - cfram_total ;
lines(:,5) = cfram_total
obs_total = obs_region(5:35)
obs_total = obs_total - obs_region(4)
; lines(:,5) = obs_total - cfram_total
printVarSummary(lines)


;=======================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_surrounding/20210520/"
figname = "Decompose_downS_cfram_clearsky_EAST"

Left = (/"(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)","(j)","(k)","(l)"/)

wks = gsn_open_wks("pdf", figpath+figname)
gsn_define_colormap(wks, "temp_diff_18lev")
plot = new(6, graphic)
plotc = new(1, graphic)

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

res@trXMinF = yrs-1
res@trXMaxF = yre+1
res@trYMinF = -1
res@trYMaxF = 12

res@gsnXYBarChart = True
res@gsnYRefLine = 0
res@gsnAboveYRefLineBarColors = "orangered"
res@gsnBelowYRefLineBarColors = "dodgerblue"
; res@gsnXYBarChartColors = "orangered"
; res@gsnXYBarChartColors2 = "dodgerblue"

; Colors = (/"black","purple"/)
; Patterns = (/0,16/)
; res@xyLineThicknessF = 10
; res@xyLineColors = Colors
; res@xyDashPatterns = Patterns

do i = 0, 5
    res@gsnLeftString = Left(i)+" "+process(i)
    plot(i) = gsn_csm_xy(wks, year_new, lines(:,i), res)
end do

delete(res@gsnAboveYRefLineBarColors)
delete(res@gsnBelowYRefLineBarColors)
; res@xyCurveDrawOrder = "PostDraw"
res@gsnLeftString = ""
res@gsnRightString = ""
res@gsnXYBarChartColors = "black"
res@gsnXYBarChartPatterns = 17
res@gsnXYBarChartFillScaleF = 0.7
res@gsnXYBarChartFillDotSizeF = 0.0025
plotc(0) = gsn_csm_xy(wks, year_new, obs_total, res)
overlay(plot(5), plotc(0))

pres = True
pres@gsnMaximize = True
pres@gsnPanelMainString = "Decompose surface downward SW: East region"
gsn_panel(wks, plot, (/2,3/), pres)
end