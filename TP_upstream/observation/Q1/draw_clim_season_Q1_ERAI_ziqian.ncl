
begin

yrStrt = 1979
yrLast = 2014
year = ispan(yrStrt, yrLast, 1)

diri = "/home/yangsong3/wangzq/DATA/"
fili = "ERAI.mon.Q1.1979-2014.computation.nc" ;  

f = addfile(diri+fili, "r")
TIME = f->time
YYYY = cd_calendar(TIME, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

Q1 = f->Q1(iYYYY,:,:)    ; units W/m^2
printVarSummary(Q1)

;---------------------------------------------------------

Season = (/"MAM","JJA","SON","DJF"/)

q1Sea = month_to_seasonN(Q1, Season)
printVarSummary(q1Sea)
q1sea = dim_avg_n_Wrap(q1Sea, 1)
;q1sea = q1Sea(:,22,:,:)
printVarSummary(q1sea)    

;========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_upstream/climatology/"
figname = "draw_clim_season_Q1_ERAI_ziqian"

wks = gsn_open_wks("pdf", figpath+figname)
nplot = 4
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

res@mpOutlineOn = True
res@mpFillOn = True
res@mpGeophysicalLineThicknessF = 2

res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False

res@gsnStringFontHeightF = 0.025
res@tmXBLabelFontHeightF = 0.02
res@tmYLLabelFontHeightF = 0.02

res@lbLabelBarOn = False
res@cnFillPalette = "cmp_b2r"
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 50.
res@cnMaxLevelValF = 400.
res@cnMinLevelValF = -400.

res@gsnRightString = "W/m~S~2~N~"
do i = 0, dimsizes(Season)-1
    res@gsnLeftString = Season(i)
    plot(i) = gsn_csm_contour_map(wks, q1sea(i,:,:), res)
    plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(i), plot_topo(i))
end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelMainString = "Q1 climatology"
pres@gsnStringFontHeightF = 0.035
pres@gsnPanelLabelBar = True

gsn_panel(wks, plot, (/2,2/), pres)

latstr = 20
latend = 50
lonstr = 60
lonend = 120
res@mpMinLatF = latstr
res@mpMaxLatF = latend
res@mpMinLonF = lonstr
res@mpMaxLonF = lonend

do i = 0, dimsizes(Season)-1
    res@gsnLeftString = Season(i)
    plot(i) = gsn_csm_contour_map(wks, q1sea(i,:,:), res)
    plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(i), plot_topo(i))
end do
gsn_panel(wks, plot, (/2,2/), pres)

end