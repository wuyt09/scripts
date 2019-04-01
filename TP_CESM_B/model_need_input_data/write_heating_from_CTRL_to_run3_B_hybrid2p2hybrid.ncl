begin


;===============================================================================================
; read data, choose 250-300 year for calculation
;===============================================================================================
diri1 = "/home/yangsong3/data-model/wuyt/TP_CTRL_B2000_F19G16_CAM4/pro/"
cname = "B2000_f19g16_CP_CTRL.cam.h0.DTCOND.200-349.nc"  ; control run
dname = "B2000_f19g16_CP_CTRL.cam.h0.DTV.200-349.nc"  ; control run
caname   = "B2000_f19g16_CP_CTRL.cam.h0.hyam.30001.nc"
cbname   = "B2000_f19g16_CP_CTRL.cam.h0.hybm.30001.nc"
cpsname  = "B2000_f19g16_CP_CTRL.cam.h0.PS.200-349.nc"

fc = addfile(diri1+cname, "r")
fd = addfile(diri1+dname, "r")
xc = fc->DTCOND(600:1199,:,:,:)  ; 25001-30012
xd = fd->DTV(600:1199,:,:,:)
printVarSummary(xc)
printVarSummary(xd)
x = xc ; copy meta data
x = xc + xd
printVarSummary(x)

;time = fc->time(600:1199)
;TIME = ut_calendar(time, -1)/100
;print(TIME)

fca = addfile(diri1+caname, "r")
fcb = addfile(diri1+cbname, "r")
fcsp = addfile(diri1+cpsname, "r")
cnhyam = fca->hyam
cnhybm = fcb->hybm
cnps   = fcsp->PS(600:1199,:,:)
printVarSummary(cnps)
printVarSummary(cnhyam)

;--------------------------------------------------------------------------
diri2 = "/home/yangsong3/data-model/wuyt/TP_NoTOPO_B2000_F19G16_CAM4/pro/"
raname   = "TP_noTOPO_hybrid.hyam.mon.30001.nc"
rbname   = "TP_noTOPO_hybrid.hybm.mon.30001.nc"
rpsname  = "TP_noTOPO_hybrid.PS.mon.25001-35001.nc"
fra = addfile(diri2+raname, "r")
frb = addfile(diri2+rbname, "r") 
frsp = addfile(diri2+rpsname, "r")

crhyam = fra->hyam
crhybm = frb->hybm
crps   = frsp->PS(0:599,:,:) ; 25001-30012

;time2 = frsp->time(0:599)
;TIME2 = ut_calendar(time2, -1)
;print(TIME2)


p0  = 1000.
p0@units = "hPa"

;plevo = (/1000, 975, 950, 925, 900, 875, 850, 825, 800, 775, 750, 700,\
;         650, 600, 550, 500, 450, 400, 350, 300, 250, 225, 200, 175,\
;         150, 125, 100, 70, 50, 30, 20, 10, 7, 5, 3, 2, 1/)
plevo = (/1, 2, 3, 5, 7, 10, 20, 30, 50, 70, 100, 125, 150, 175, 200,\
         225, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750,\
         775, 800, 825, 850, 875, 900, 925, 950, 975, 1000/)
plevo@units = "hPa"

xp = vinth2p(x, cnhyam, cnhybm, plevo, cnps, 1, p0, 1, False)
xp!1 = "lev"
xp&lev = plevo
xp@_FillValue = -999.  ; ;0.
printVarSummary(xp)

z = pres2hybrid_Wrap(plevo*100, crps, p0*100, xp, crhyam, crhybm, 0)
; careful for the plevo, it requires top-to-bottom lev
printVarSummary(z)

zp = vinth2p(z, crhyam, crhybm, plevo, crps, 1, p0, 1, False)
zp!1 = "lev"
zp&lev = plevo
printVarSummary(zp)

;===============================================================
; write data
;===========================================================
z1 = clmMonTLLL(z)
z1@long_name = "monthly sum of heating "
z1@_FillValue = -999.
z1 = where(abs(z1).gt.0.002, z1@_FillValue, z1)
z1@_FillValue = -999.
printVarSummary(z1)

dataname = "heating_B_H2P2H_C20190304.nc"
datapath = "./"
system("rm -f "+datapath+dataname)
output = addfile(datapath+dataname, "c")
output->heat  = z1
printVarSummary(z1)
exit

;===========check if interpolation works=========================
x = x*86400
z = z*86400
xp = xp*86400
zp = zp*86400

xtest = x(50,:,{27:40},{80:110})
ztest = z(50,:,{27:40},{80:110})
testt = ztest - xtest
print("clm control heating = "+xtest+" K/day")
print("clm run1111 heating = "+ztest+" K/day")
print("diff heating = "+testt+" K/day")

p = xp(0:2,:,{35},{90})
p(0,:) = xp(50,:,{35},{90})
p(1,:) = zp(50,:,{35},{90})
p(2,:) = p(0,:) - p(1,:)
printVarSummary(p)

pp = xp(0:2,:,{35},{90})
pp(0,:) = dim_avg_n_Wrap(xp(50,:,{25:40},{75:105}),(/1,2/))
pp(1,:) = dim_avg_n_Wrap(zp(50,:,{25:40},{75:105}),(/1,2/))
pp(2,:) = pp(0,:) - pp(1,:)
printVarSummary(pp)
delete(xc)
delete(xd)
delete(x)
delete(z)


;===============================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_B/model_need_input_data/"
figname = "write_heating_from_CTRL_to_run3_B_hybrid2hybrid"
wks = gsn_open_wks("ps", figpath+figname)
print("wks")
;==================================================================
; spatial distribution at specfic levels
;==================================================================
res  = True
res@gsnFrame  = False
res@gsnDraw   = False

nplot = 6
plot = new(nplot,graphic)
plot_topo = plot
gsn_define_colormap(wks,"cmp_flux")
print("nplot")

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
print("res_topo")

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

print("res")
printMinMax(xp(50,{850},:,:), 1)
printMinMax(zp(50,{850},:,:), 1)

res@gsnLeftString = "850hPa_before_interpolation"
plot(0) = gsn_csm_contour_map(wks, xp(50,{850},:,:), res)
plot_topo(0) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(0), plot_topo(0))
res@gsnLeftString = "850hPa_after_interpolation"
plot(1) = gsn_csm_contour_map(wks, zp(50,{850},:,:), res)
plot_topo(1) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(1), plot_topo(1))
print("end 1-2")


printMinMax(xp(50,{500},:,:), 1)
printMinMax(zp(50,{500},:,:), 1)
res@gsnLeftString = "500hPa"
plot(2) = gsn_csm_contour_map(wks, xp(50,{500},:,:), res)
plot_topo(2) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(2), plot_topo(2))
plot(3) = gsn_csm_contour_map(wks, zp(50,{500},:,:), res)
plot_topo(3) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(3), plot_topo(3))
print("end 3-4")


printMinMax(xp(50,{200},:,:), 1)
printMinMax(zp(50,{200},:,:), 1)
res@gsnLeftString = "200hPa"
plot(4) = gsn_csm_contour_map(wks, xp(50,{200},:,:), res)
plot_topo(4) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(4), plot_topo(4))
plot(5) = gsn_csm_contour_map(wks, zp(50,{200},:,:), res)
plot_topo(5) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(5), plot_topo(5))
print("end 5-6")

pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
gsn_panel(wks, plot, (/3,2/), pres)

;==================================================================
; profile
;==================================================================

plots = new(2, graphic)
Colors = (/"red","blue","black"/) ;,"purple"
Texts  = (/"Ctrl_run","No_topo","No_topo-Ctrl_run"/)  ;,"cltr-run-wgt","No-TOPO-wgt"
dum5 = new(3,graphic)
dum6 = dum5

resxy = True
resxy@gsnDraw  = False
resxy@gsnFrame = False

resxy@trYReverse = True
resxy@xyYStyle = "Log"
resxy@tmYLMode = "Explicit"
resxy@tmYLValues = (/ 1000, 700,500,400,300,200,100,50,30,10/)
resxy@tmYLLabels = ""+resxy@tmYLValues

resxy@xyLineThicknessF = 3.0
resxy@xyLineColors = Colors
resxy@xyDashPatterns = (/0,0/) ;(/0,0,1,1/)
resxy@trYMinF   = 10
resxy@trYMaxF   = 1000
resxy@trXMinF   = -2
resxy@trXMaxF   = 5
resxy@gsnXRefLine = 0
resxy@gsnXRefLineThicknessF = 3
resxy@gsnXRefLineColor = "grey"
resxy@gsnXRefLineDashPattern = 1
resxy@tiMainString = "sum of heating (35N, 90E)"
resxy@tiYAxisString = "hPa"

plots(0) = gsn_csm_xy(wks, p, p&lev, resxy) 

resxy@tiMainString = "sum of heating (entire TP)"
plots(1) = gsn_csm_xy(wks, pp, p&lev, resxy)


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
        dum5(i) = gsn_add_polyline(wks, plots(0), xx, yy, res_lines)
        dum6(i) = gsn_add_text(wks, plots(0), Texts(i), text_pos, yy(0), res_text)
    end do

Pres = True
Pres@gsnMaximize = True
gsn_panel(wks, plots, (/2,1/), Pres)

end
