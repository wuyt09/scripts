
begin

lonstr = 60
lonend = 120
latstr = 10
latend = 60


diri = "/home/yangsong3/wuyt/sysu/scripts/TP_CESM_B/model_need_input_data/"
fili = "heating_B_H2P2H_C20190304.nc"

f = addfile(diri+fili, "r")
heating = f->heat
heating = heating*86400
heating@units = "K/day"
printVarSummary(heating)
;print(heating&lat)
;print(heating&lon)

lon = heating&lon
lat = heating&lat
lev = heating&lev
printVarSummary(lon)
month = ispan(1, 12, 1)

system("rm -f ./heating_input_hybrid.txt")
write_table("./heating_input_hybrid.txt", "w", [/"mon","lev","lat","lon","heat"/], "%8s %8s %8s %8s %8s")

do i = 0, 11
    do k = 61,68  ;25.57-38.84
        do l = 30,42 ;75-105
            do j = 0,25
                write_table("./heating_input_hybrid.txt", "a", [/month(i),lev(j),lat(k),lon(l),heating(i,j,k,l)/],\
                "%8.3i %8.3f %8.3f %8.3f %8.3f")
                ;asciiwrite("./heating_input_hybrid.txt", heating(i,j,k,l))
            end do
        end do
    end do
end do


diri2 = "/home/yangsong3/data-model/wuyt/TP_NoTOPO_B2000_F19G16_CAM4/pro/"
raname   = "TP_NOTOPO.hyam.mon.30001.nc"
rbname   = "TP_NOTOPO.hybm.mon.30001.nc"
rpsname  = "TP_NOTOPO.PS.mon.26001-30508.nc"
fra = addfile(diri2+raname, "r")
frb = addfile(diri2+rbname, "r") 
frsp = addfile(diri2+rpsname, "r")

hyam = fra->hyam
hybm = frb->hybm
ps0   = frsp->PS(0:539,:,:)
ps = clmMonTLL(ps0)
printVarSummary(ps)

;plevo = (/1000, 975, 950, 925, 900, 875, 850, 825, 800, 775, 750, 700,\
;         650, 600, 550, 500, 450, 400, 350, 300, 250, 225, 200, 175,\
;         150, 125, 100, 70, 50, 30, 20, 10, 7, 5, 3, 2, 1/)
plevo = (/1, 2, 3, 5, 7, 10, 20, 30, 50, 70, 100, 125, 150, 175, 200,\
         225, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750,\
         775, 800, 825, 850, 875, 900, 925, 950, 975, 1000/)
plevo@units = "hPa"
p0 = 1000.
p0@units = "hPa"

;heating@_FillValue = 0.
heat = vinth2p(heating, hyam, hybm, plevo, ps, 1, p0, 1, False)
printVarSummary(heat)

;=======================================================
maskpath = "/home/yangsong3/wuyt/sysu/scripts/TP_CESM_B/model_need_input_data/"
maskname = "TPmask_C20190304.nc"
fmask = addfile(maskpath+maskname, "r")
tp = fmask->tp
tp0 = conform_dims(dimsizes(heat), tp, (/2,3/))
printVarSummary(tp0)
heat = mask(heat, tp0.ge.0.9, True)
print(heat(0,4,{35},{90}))
print(heat(0,4,{20},{90}))

;h_profile = dim_avg_n_Wrap(heat(:,:,{25:40},{75:105}), (/2,3/))
h_profile = wgt_areaave_Wrap(heat(:,:,{20:45},{65:110}), 1, 1, 1)
printVarSummary(h_profile)

;wgt1 = heat(0,0,0,{75:105})
;wgt2 = heat(0,0,{25:40},0)
;wgt1 = 1
;wgt2 = 1 
;h_p1 = dim_avg_wgt_n_Wrap(heat(:,:,:,{75:105}), wgt1, 1, 3)
;h_p2 = dim_avg_wgt_n_Wrap(h_p1(:,:,{25:40}), wgt2, 1, 2)
;h_profile = h_p2
;printVarSummary(h_profile)

;===============================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_B/model_need_input_data/"
figname = "check_heating_input_spatial_profile_12month"
wks = gsn_open_wks("ps", figpath+figname)

;==================================================================
; spatial distribution at specfic levels
;==================================================================
res  = True
res@gsnFrame  = False
res@gsnDraw   = False

nplot = 12
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
res@mpMinLatF    = latstr
res@mpMaxLatF    = latend
res@mpMinLonF    = lonstr
res@mpMaxLonF    = lonend
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

do i = 0,11
    res@gsnLeftString = "MONTH = "+(i+1)
    plot(i) = gsn_csm_contour_map(wks, heat(i,{500},:,:), res)
    plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(i), plot_topo(i))
end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
pres@txString = "500hPa"
gsn_panel(wks, plot, (/3,4/), pres)

delete(plot)
delete(plot_topo)
plot = new(nplot,graphic)
plot_topo = plot

do i = 0,11
    res@gsnLeftString = "MONTH = "+(i+1)
    plot(i) = gsn_csm_contour_map(wks, heat(i,{400},:,:), res)
    plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(i), plot_topo(i))
end do
pres@txString = "400hPa"
gsn_panel(wks, plot, (/3,4/), pres)


;==================================================================
; profile
;==================================================================
delete(res)
plots = new(nplot, graphic)

 res = True
 res@gsnDraw  = False
 res@gsnFrame = False

 res@gsnXRefLine = 0
 res@gsnXRefLineDashPattern = 2
 res@gsnXRefLineColor = "grey"
 res@gsnXRefLineThicknessF = 3
 res@trXMinF = -6.
 res@trXMaxF = 6
 res@trYMinF   = 100
 res@trYMaxF   = 1000
 res@trYLog = True
 
 res@tmYLMode      = "Explicit"
 res@tmYLValues    = (/"1000","850","750","500","300","100"/)
 res@tmYLLabels    = (/"1000","850","750","500","300","100"/)
 res@xyLineColors = "red"
 res@xyLineThicknessF = 3
 res@xyDashPattern = 0

 res@trYReverse = True
 
 res@gsnRightString = "K"
 res@gsnStringFontHeightF = 0.035

do i = 0, 11
    res@gsnLeftString = "MONTH = "+(i+1)
    plots(i) = gsn_csm_xy(wks, h_profile(i,:), h_profile&lev_p, res) 
end do

Pres = True
Pres@gsnMaximize = True
gsn_panel(wks, plots, (/3,4/), Pres)

end