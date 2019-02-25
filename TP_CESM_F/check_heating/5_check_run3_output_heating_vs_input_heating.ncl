
begin

diri = "/home/yangsong3/data-model/wuyt/TP_NoTOPO_heat_F2000_F19G16_CAM4/pro/"
cname = "TP_noTopo_heat.DTCOND.mon.0601-3012.nc"
dname = "TP_noTopo_heat.DTV.mon.0601-3012.nc"

aname = "TP_noTopo_heat.hyam.mon.0601-3012.nc"
bname = "TP_noTopo_heat.hybm.mon.0601-3012.nc"

pname = "TP_noTopo_heat.PS.mon.0601-3012.nc"

fa = addfile(diri+aname, "r")
fb = addfile(diri+bname, "r")
fc = addfile(diri+cname, "r")
fd = addfile(diri+dname, "r")
fp = addfile(diri+pname, "r")

hyam = fa->hyam
hybm = fb->hybm
h1 = fc->DTCOND
h2 = fd->DTV
ps = fp->PS

h = h1
h = h1+h2
printVarSummary(h)


plevo = (/1000, 975, 950, 925, 900, 875, 850, 825, 800, 775, 750, 700,\
         650, 600, 550, 500, 450, 400, 350, 300, 250, 225, 200, 175,\
         150, 125, 100, 70, 50, 30, 20, 10, 7, 5, 3, 2, 1/)
plevo@units = "hPa"
p0  = 1000.
p0@units = "hPa"
hp = vinth2p(h, hyam, hybm, plevo, ps, 1, p0, 1, False)
hp = hp*86400

hp1 = clmMonTLLL(hp)
printVarSummary(hp)

;---------------------------------------

diri2 = "/home/yangsong3/data-model/wuyt/TP_NoTOPO_heat_F2000_F19G16_CAM4/pro/"
fili = "heat_float_input.nc"
f = addfile(diri2+fili, "r")
hh = f->heat
hh = hh*86400
printVarSummary(hh)

h = h*86400

;=======plot=================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/check_heating/"
figname = "5_check_run3_output_heating_vs_input_heating"

wks = gsn_open_wks("pdf", figpath+figname)
res  = True
res@gsnFrame  = False
res@gsnDraw   = False

nplot = 6
plot = new(nplot,graphic)
plot_topo = plot
gsn_define_colormap(wks,"cmp_flux")

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
 res_topo@cnLineThicknessF       = 5
 res_topo@cnLineColor            = "black"
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

res@cnLevelSelectionMode = "ManualLevels"
res@cnMaxLevelValF  = 5  ;5e-5
res@cnMinLevelValF  = -5  ;-5e-5
res@cnLevelSpacingF = 0.5  ;5e-6

res@gsnStringFontHeightF = 0.035

res@gsnLeftString = "850hPa_run3_output"
plot(0) = gsn_csm_contour_map(wks, h(2,{850},:,:), res)
plot_topo(0) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(0), plot_topo(0))
res@gsnLeftString = "850hPa_input"
plot(1) = gsn_csm_contour_map(wks, hh(2,{850},:,:), res)
plot_topo(1) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(1), plot_topo(1))


res@gsnLeftString = "500hPa"
plot(2) = gsn_csm_contour_map(wks, h(2,{500},:,:), res)
plot_topo(2) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(2), plot_topo(2))
plot(3) = gsn_csm_contour_map(wks, hh(2,{500},:,:), res)
plot_topo(3) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(3), plot_topo(3))


res@gsnLeftString = "200hPa"
plot(4) = gsn_csm_contour_map(wks, h(2,{200},:,:), res)
plot_topo(4) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(4), plot_topo(4))
plot(5) = gsn_csm_contour_map(wks, hh(2,{200},:,:), res)
plot_topo(5) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(5), plot_topo(5))

pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
gsn_panel(wks, plot, (/3,2/), pres)



;========================================================================================
;======================================================
hp1@_FillValue = -999.
p = h(0:2,:,{35},{90})
p(0,:) = h(2,:,{35},{90})
p(1,:) = hh(2,:,{35},{90})
p(2,:) = hh(2,:,{35},{90}) - h(2,:,{35},{90})
printVarSummary(p)


Colors = (/"red","blue","black"/) ;,"purple"
Texts  = (/"Ctrl_run","No_topo","No_topo-Ctrl_run"/)  ;,"cltr-run-wgt","No-TOPO-wgt"
dum5 = new(3,graphic)
dum6 = dum5

resxy = True
resxy@gsnDraw  = False
resxy@gsnFrame = False

resxy@trYReverse = True
resxy@xyYStyle = "Log"
resxy@tmYLMode          = "Explicit"                ; explicit labels
resxy@tmYLValues        = (/1000,700,500,400,300,200,100,50/)
resxy@tmYLLabels        = ""+resxy@tmYLValues         ; make strings

resxy@xyLineThicknessF = 3.0
resxy@xyLineColors = Colors
resxy@xyDashPatterns = (/0,0/) ;(/0,0,1,1/)

resxy@trXMinF   = -2
resxy@trXMaxF   = 5
resxy@gsnXRefLine = 0
resxy@gsnXRefLineThicknessF = 3
resxy@gsnXRefLineColor = "grey"
resxy@gsnXRefLineDashPattern = 1
resxy@tiMainString = "sum of heating (35N, 90E)"
resxy@tiYAxisString = "hPa"

plots = gsn_csm_xy(wks, p, p&lev , resxy)  ; p(:,{1000:50}), plevo(0:28)
; lev/100: Pa converts to hPa
; *86400: K/s converts to K/day


    res_lines                   = True                  ; polyline mods desired
    res_lines@gsLineDashPattern = 0                     ; solid line
    res_lines@gsLineThicknessF  = 3                     ; line thicker

    ;res_lines@gsMarkerColor     = "blue"
    ;res_lines@gsMarkerIndex     = 16
    ;res_lines@gsMarkerSizeF     = 0.01
    
    res_text                    = True                  ; text mods desired
    res_text@txFontHeightF      = 0.018                 ; change text size
    res_text@txJust             = "CenterLeft"          ; text justification

        ;yy = (/100,100/)
    xx = (/2.0,2.5/)
    text_pos = 2.6
    do i = 0,2
        yy = (/100+i*50,100+i*50/)
        res_lines@gsLineColor       = Colors(i)
        dum5(i) = gsn_add_polyline(wks, plots, xx, yy, res_lines)
        dum6(i) = gsn_add_text(wks, plots, Texts(i), text_pos, yy(0), res_text)
    end do

draw(plots)
frame(wks)

end