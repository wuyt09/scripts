begin


;===============================================================================================
; read data, choose 250-300 year for calculation
;===============================================================================================
diri1 = "/home/yangsong3/data-model/wuyt/TP_CTRL_B2000_F19G16_CAM4/pro/"
cname = "B2000_f19g16_CP_CTRL.cam.h0.DTCOND.200-349.nc"  ; control run
dname = "B2000_f19g16_CP_CTRL.cam.h0.DTV.200-349.nc"  ; control run
caname   = "B2000_f19g16_CP_CTRL.cam.h0.hyam.200-349.nc"
cbname   = "B2000_f19g16_CP_CTRL.cam.h0.hybm.200-349.nc"
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
cnhyam = fca->hyam(600)
cnhybm = fcb->hybm(600)
cnps   = fcsp->PS(600:1199,:,:)
printVarSummary(cnps)
printVarSummary(cnhyam)

;--------------------------------------------------------------------------
diri2 = "/home/yangsong3/data-model/wuyt/TP_NoTOPO_B2000_F19G16_CAM4/pro/"
raname   = "TP_noTOPO_hybrid.hyam.mon.25001-35001.nc"
rbname   = "TP_noTOPO_hybrid.hybm.mon.25001-35001.nc"
rpsname  = "TP_noTOPO_hybrid.PS.mon.25001-35001.nc"
fra = addfile(diri2+raname, "r")
frb = addfile(diri2+rbname, "r") 
frsp = addfile(diri2+rpsname, "r")

crhyam = fra->hyam(0)
crhybm = frb->hybm(0)
crps   = frsp->PS(0:599,:,:) ; 25001-30012

;time2 = frsp->time(0:599)
;TIME2 = ut_calendar(time2, -1)
;print(TIME2)


p0  = 100000.
p0@units = "Pa"
z = hyi2hyo_Wrap(p0, cnhyam, cnhybm, cnps, x, crhyam, crhybm, 0)
printVarSummary(z)

z1 = clmMonTLLL(z)
;z0 = todouble(z1)
z1@long_name = "monthly sum of heating "
;copy_VarMeta(z1, z0)
printVarSummary(z1)

print("Z = "+z(6,:,{35},{90})*86400)
print("X = "+x(6,:,{35},{90})*86400)
exit
;===============================================================
;file1 = "heating.mon.nc"
;system("rm -f "+diri+file1)
;output = addfile(diri+file1, "c")
;;output->month = month
;;output->lev   = (/z1&lev/)
;;output->lat   = (/z1&lat/)
;;output->lon   = (/z1&lon/)
;output->heat  = z0
;printVarSummary(z0)

;=================================================================
xtest = x(50,:,{35},{90})*86400 ;dim_avg_n_Wrap(x, 0)
ztest = z(50,:,{35},{90})*86400 ;dim_avg_n_Wrap(z, 0)
testt = ztest - xtest
print("clm control heating = "+xtest+" K/day")
print("clm run1111 heating = "+ztest+" K/day")
print("diff heating = "+testt+" K/day")
;print(xc&lev)
;print(xd&lev)

;plevo = (/1000,950,900,850,800,750,700,650,600,550,500,450,400,350,300,250,200,150,100,50,10,5,1/)
plevo = (/1,5,10,50,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000/)
plevo@units = "hPa"
p0 = 1000.
p0@units = "hPa"
y1 = vinth2p(x, cnhyam, cnhybm, plevo, cnps, 1, p0, 1, False)
y2 = vinth2p(z, crhyam, crhybm, plevo, crps, 1, p0, 1, False)
printVarSummary(y1)
printVarSummary(y2)

;p = new((/3,dimsizes(plevo),dimsizes(y1&lat),dimsizes(y1&lon)/), float)
p = y1(0:2,:,{35},{90})
p(0,:) = y1(50,:,{35},{90})
p(1,:) = y2(50,:,{35},{90})
p(2,:) = y2(50,:,{35},{90}) - y1(50,:,{35},{90})
;copy_VarMeta(y1, p(0,:,:,:))
;p!1 = "level"
;p&level = plevo
;p&level@units = "hPa"
printVarSummary(p)
print(p*86400)

;===============================================================
wks = gsn_open_wks("ps", "cal_monthly_diffusion_condensation_heating_hybrid2hybrid")

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
res@gsnAddCyclic    = False
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

;res@cnLevelSelectionMode = "ManualLevels"
;res@cnMaxLevelValF  = 5  ;5e-5
;res@cnMinLevelValF  = -5  ;-5e-5
;res@cnLevelSpacingF = 0.5  ;5e-6

res@gsnStringFontHeightF = 0.035

res@gsnLeftString = "850hPa"
plot(0) = gsn_csm_contour_map(wks, y1(50,{850},:,:), res)
plot_topo(0) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(0), plot_topo(0))
plot(1) = gsn_csm_contour_map(wks, y2(50,{850},:,:), res)
plot_topo(1) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(1), plot_topo(1))


res@gsnLeftString = "500hPa"
plot(2) = gsn_csm_contour_map(wks, y1(50,{500},:,:), res)
plot_topo(2) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(2), plot_topo(2))
plot(3) = gsn_csm_contour_map(wks, y2(50,{500},:,:), res)
plot_topo(3) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(3), plot_topo(3))


res@gsnLeftString = "200hPa"
plot(4) = gsn_csm_contour_map(wks, y1(50,{200},:,:), res)
plot_topo(4) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(4), plot_topo(4))
plot(5) = gsn_csm_contour_map(wks, y2(50,{200},:,:), res)
plot_topo(5) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(5), plot_topo(5))


pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
gsn_panel(wks, plot, (/3,2/), pres)

;========================================================================================
;p = new((/3,dimsizes(plevo)/), float)
;y1test = dim_avg_n_Wrap(y1, 0)
;y2test = dim_avg_n_Wrap(y2, 0)
;diff0 = y2test(:,{35},{90}) - y1test(:,{35},{90})
;p(0,:) = y1test(:,{35},{90}) ;*86400
;p(1,:) = y2test(:,{35},{90}) ;*86400
;p(2,:) = diff0
;print(p)


Colors = (/"red","blue","black"/) ;,"purple"
Texts  = (/"Ctrl_run","No_topo","No_topo-Ctrl_run"/)  ;,"cltr-run-wgt","No-TOPO-wgt"
dum5 = new(3,graphic)
dum6 = dum5

resxy = True
resxy@gsnDraw  = False
resxy@gsnFrame = False

resxy@trYReverse = True
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

plots = gsn_csm_xy(wks, p*86400, p&lev_p, resxy)  ;/100, 
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
