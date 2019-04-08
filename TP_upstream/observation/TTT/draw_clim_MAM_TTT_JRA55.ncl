
begin

yrStrt = 1979
yrLast = 2016
year = ispan(yrStrt, yrLast, 1)

diri = "/home/yangsong3/wangzq/DATA/jra/"
fili = "air.jra55.mon.mean.nc"

f = addfile(diri+fili, "r")
TIME = f->time
YYYY = cd_calendar(TIME, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

T = short2flt(f->air(iYYYY,:,:,:))
printVarSummary(T)

;===========integration========================================
lev = T&level
plev = lev*100.
psfc = max(plev)
ptop = min(plev)

dP = dpres_plevel_Wrap(plev, psfc, ptop, 0)
dP!0 = "lev"
dP&lev = lev
printVarSummary(dP)

dp = conform_dims(dimsizes(T), dP, 1)
Tdp = T*dp
copy_VarCoords(T, Tdp)
TTT = dim_sum_n_Wrap(Tdp(:,{400:200},:,:), 1)/sum(dP({400:200}))
copy_VarCoords(T(:,0,:,:), TTT)
printVarSummary(TTT)
;===============================================================

t = month_to_season(TTT, "MAM")
printVarSummary(t)
print(t(:,{25},{75}))

;wgt = year
;wgt = 1
;ttt = dim_avg_wgt_n_Wrap(t, wgt, 1, 0)
ttt = dim_avg_n_Wrap(t, 0)
print(ttt({25},{75}))

;=================================================================

latstr = 20
latend = 50
lonstr = 60
lonend = 120

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_upstream/observation/TTT/"
figname = "draw_clim_MAM_TTT_JRA55"

wks = gsn_open_wks("pdf", figpath+figname)
nplot = 1
plot = new(nplot, graphic)
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
 res_topo@cnLineColor            = "gray20"
 res_topo@cnInfoLabelOn          = False
 res_topo@cnLevelSelectionMode   ="ExplicitLevels"
 res_topo@cnLevels               = (/1500/)



res = True
res@gsnDraw  = False
res@gsnFrame = False
res@gsnMaximize = True

res@mpOutlineOn = True
res@mpFillOn = True
res@mpGeophysicalLineThicknessF = 2

res@mpMinLatF = latstr
res@mpMaxLatF = latend
res@mpMinLonF = lonstr
res@mpMaxLonF = lonend

res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False

res@gsnStringFontHeightF = 0.025
res@tmXBLabelFontHeightF = 0.02
res@tmYLLabelFontHeightF = 0.02

;res@lbLabelBarOn = False
res@cnFillPalette = "temp_19lev"
;res@cnLevelSelectionMode = "ManualLevels"
;res@cnLevelSpacingF = 0.1
;res@cnMaxLevelValF = 1.
;res@cnMinLevelValF = -1.

res@gsnRightString = "K"
res@gsnLeftString  = "MAM tropo T"

plot(0) = gsn_csm_contour_map(wks, ttt, res)
plot_topo(0) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(0), plot_topo(0))

draw(plot)
end