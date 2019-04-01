
begin

diri = "./"
fili = "TP_NOTOPO.cam.h0.0250-*.nc"
filis = systemfunc("ls "+diri+fili)
print(filis)

f = addfiles(filis, "r")

x1 = f[:]->DTCOND
x2 = f[:]->DTV
hyam = f[0]->hyam
hybm = f[0]->hybm
ps = f[:]->PS
printVarSummary(ps)

x = x1
x = x1+x2
x = x*86400
printVarSummary(x)

p0  = 1000.
p0@units = "hPa"
plevo = (/1, 2, 3, 5, 7, 10, 20, 30, 50, 70, 100, 125, 150, 175, 200,\
         225, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750,\
         775, 800, 825, 850, 875, 900, 925, 950, 975, 1000/)
plevo@units = "hPa"
xp = vinth2p(x, hyam, hybm, plevo, ps, 1, p0, 1, False)
printVarSummary(xp)

xp_profile = dim_avg_n_Wrap(xp(:,:,{25:40},{75:105}), (/2,3/))
printVarSummary(xp_profile)

;=================================================================================

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_B/NOTOPO/check/"
figname = "check_DTCOND_DTV"

wks = gsn_open_wks("ps", figpath+figname)

;==================================================================
; spatial distribution at specfic levels
;==================================================================

nplot = 12
plot = new(nplot,graphic)
plot_topo = plot
gsn_define_colormap(wks,"WhiteYellowOrangeRed")

 ;------------------------------------------------------------
 ; TOPO
 ;------------------------------------------------------------
 
 ff = "~/wuyt/data/TOPO_1500.nc"
 topo_in = addfile(ff,"r")
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
 res_topo@cnLineThicknessF       = 2
 res_topo@cnLineColor            = "grey4"
 res_topo@cnInfoLabelOn          = False
 res_topo@cnLevelSelectionMode   ="ExplicitLevels"
 res_topo@cnLevels               = (/1500/)

;---------------------------------------

res                 =  True
res@gsnDraw         = False
res@gsnFrame        = False
res@gsnAddCyclic    = True
;res@gsnSpreadColors = True

res@mpOutlineOn  = True
res@mpFillOn     = True
;res@mpCenterLonF = 0.
;res@mpMinLatF    = latstr
;res@mpMaxLatF    = latend
;res@mpMinLonF    = lonstr
;res@mpMaxLonF    = lonend
;res@mpCenterLonF = 180
res@mpGeophysicalLineThicknessF = 2

res@cnFillOn     = True
res@cnLinesOn    = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn  = False
res@lbLabelBarOn   = False

res@cnLineLabelsOn         = False
res@cnLineThicknessF       = 5
res@cnLineColor            = "black"
res@cnInfoLabelOn          = False
;res@cnLevelSelectionMode   ="ExplicitLevels"
;res@cnLevels               = (/0.9/)
res@cnLevelSelectionMode = "ManualLevels"
res@cnMaxLevelValF = 7
res@cnMinLevelValF = 0
res@cnLevelSpacingF = 1

res@gsnStringFontHeightF = 0.035

do i = 0,11
    res@gsnLeftString = "MONTH = "+(i+1)
    plot(i) = gsn_csm_contour_map(wks, xp(i,{850},:,:), res)
    plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(i), plot_topo(i))
end do
pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = "850hPa"
gsn_panel(wks, plot, (/3,4/), pres)

do i = 0,11
    res@gsnLeftString = "MONTH = "+(i+1)
    plot(i) = gsn_csm_contour_map(wks, xp(i,{500},:,:), res)
    plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(i), plot_topo(i))
end do
pres@gsnPanelMainString = "500hPa"
gsn_panel(wks, plot, (/3,4/), pres)

;==================================================================
; profile
;==================================================================

plots = new(nplot, graphic)

resxy = True
resxy@gsnDraw  = False
resxy@gsnFrame = False

resxy@trYReverse = True
resxy@xyYStyle = "Log"
resxy@tmYLMode = "Explicit"
resxy@tmYLValues = (/ 1000, 700,500,400,300,200,100,50,30,10/)
resxy@tmYLLabels = ""+resxy@tmYLValues

resxy@xyLineThicknessF = 3.0
resxy@xyLineColors = "red"
resxy@xyDashPatterns = (/0,0/) ;(/0,0,1,1/)
resxy@trYMinF   = 10
resxy@trYMaxF   = 1000
resxy@trXMinF   = -2
resxy@trXMaxF   = 5
resxy@gsnXRefLine = 0
resxy@gsnXRefLineThicknessF = 3
resxy@gsnXRefLineColor = "grey"
resxy@gsnXRefLineDashPattern = 1
resxy@tiYAxisString = "hPa"

do i = 0, 11
    resxy@gsnLeftString = "MONTH = "+(i+1)
    plots(i) = gsn_csm_xy(wks, xp_profile(i,:), xp_profile&lev_p, resxy) 
end do

Pres = True
Pres@gsnMaximize = True
Pres@gsnPanelMainString = "sum of heating (entire TP)"
gsn_panel(wks, plots, (/3,4/), Pres)


end

