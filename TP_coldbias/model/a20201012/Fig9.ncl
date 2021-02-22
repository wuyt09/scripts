
begin

lonstr = 60
lonend = 150
latstr = 15
latend = 50


yrStrt = 198102
yrLast = 200601
nyear = 25

casename = (/"Ctrl_FAMIP","Hist_FAMIP"/) ;,"Hist_IO","Hist_PO","Hist_AO","Hist_PIO","Hist_PAO","Hist_IAO"
ncase = dimsizes(casename)

season = "DJF"
nlat = 96
nlon = 144

LEV = 300

z = new((/ncase,nyear,nlat,nlon/), float)
z!0 = "case"
z&case = casename
v = z
u = z
h = z

do i = 0, ncase-1
    diri = "/home/yangsong3/data-model/wuyt/TPbias/TPbias_"+casename(i)+"/a20191206/"

    fili1 = casename(i)+".cam.h0.Z3.1979-2005.nc"
    fili2 = casename(i)+".cam.h0.U.1979-2005.nc"
    fili3 = casename(i)+".cam.h0.V.1979-2005.nc"
    f1 = addfile(diri+fili1, "r")
    f2 = addfile(diri+fili2, "r")
    f3 = addfile(diri+fili3, "r")

    TIME = f1->time
    YYYY = cd_calendar(TIME, -1)
    iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

    ZX = f1->Z3(iYYYY,{300},:,:)
    Z_zonal = dim_avg_n_Wrap(ZX, 2)
    Z = ZX - conform_dims(dimsizes(ZX), Z_zonal, (/0,1/))
    copy_VarCoords(ZX, Z)

    U = f2->U(iYYYY,{200},:,:)
    V = f3->V(iYYYY,{200},:,:)
    H = f1->Z3(iYYYY,{500},:,:)

    z(i,:,:,:) = month_to_season(Z, season)
    u(i,:,:,:) = month_to_season(U, season)
    v(i,:,:,:) = month_to_season(V, season)
    h(i,:,:,:) = month_to_season(H, season)
end do

printVarSummary(z)

;=============================================================================
; ttest
;=============================================================================
z_ave = dim_avg_n_Wrap(z, 1)
z_var = dim_variance_n_Wrap(z, 1)

u_ave = dim_avg_n_Wrap(u, 1)
u_var = dim_variance_n_Wrap(u, 1)

v_ave = dim_avg_n_Wrap(v, 1)
v_var = dim_variance_n_Wrap(v, 1)

h_ave = dim_avg_n_Wrap(h, 1)
; h_var = dim_variance_n_Wrap(h, 1)

probz = z_ave(0,:,:)
probz = ttest(z_ave(0,:,:), z_var(0,:,:), nyear, z_ave(1,:,:), z_var(1,:,:), nyear, False, False)

probu = u_ave(0,:,:)
probu = ttest(u_ave(0,:,:), u_var(0,:,:), nyear, u_ave(1,:,:), u_var(1,:,:), nyear, False, False)

probv = v_ave(0,:,:)
probv = ttest(v_ave(0,:,:), v_var(0,:,:), nyear, v_ave(1,:,:), v_var(1,:,:), nyear, False, False)
printVarSummary(probv)

diff_z = z_ave(0,:,:)
diff_z = z_ave(1,:,:) - z_ave(0,:,:)
diff_u = u_ave(0,:,:)
diff_u = u_ave(1,:,:) - u_ave(0,:,:)
diff_v = v_ave(0,:,:)
diff_v = v_ave(1,:,:) - v_ave(0,:,:)

diff_h = h_ave(0,:,:)
diff_h = h_ave(1,:,:) - h_ave(0,:,:)
;=========================================================
; outline of tp
;=========================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)
topo = smth9_Wrap(topo, 0.5, 0.25, True)
;=========================================================
; plot
;=========================================================

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/a20201012/"
figname = "Fig9"


wks = gsn_open_wks("pdf", figpath+figname)
plot = new(9, graphic)
plotc = plot
plot_topo = plot

 res_topo            = True
 res_topo@cnFillOn   = False
 res_topo@cnLinesOn  = True
 res_topo@gsnFrame   = False
 res_topo@gsnDraw    = False
 res_topo@cnLineLabelsOn         = False
 res_topo@cnLineThicknessF       = 3
 res_topo@cnLineColor            = "gray40"
 res_topo@cnInfoLabelOn          = False
 res_topo@cnLevelSelectionMode   ="ExplicitLevels"
 res_topo@cnLevels               = (/1500/)
 res_topo@gsnLeftString = ""
 res_topo@gsnRightString = ""

res                       = True     ; plot mods desired
res@gsnDraw               = False
res@gsnFrame              = False

res@gsnMaximize           = True     ; Make plot fit the frame
res@cnFillOn              = True     ; turn on color fill
res@cnLinesOn             = False    ; turn of contour lines
res@cnLineLabelsOn        = False
res@cnInfoLabelOn         = False
res@cnFillPalette         = "temp_diff_18lev"
; res@lbLabelBarOn          = False
res@pmLabelBarOrthogonalPosF = 0.2

res@mpMinLatF            = latstr      ; range to zoom in on
res@mpMaxLatF            = latend
res@mpMinLonF            = lonstr
res@mpMaxLonF            = lonend
; res@mpCenterLonF         = 180.
res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025

res@gsnStringFontHeightF = 0.035

; ====================================================================================
cnres = True
cnres@gsnDraw  = False
cnres@gsnFrame = False

cnres@cnFillOn = False
cnres@cnLinesOn = True
cnres@cnLineLabelsOn = True
cnres@cnLineLabelFontColor = "red"
cnres@cnLineLabelInterval = 1
; cnres@cnLineLabelDensityF = 2
cnres@cnLineLabelPlacementMode = "Constant"
cnres@cnLineDashSegLenF = 0.2

cnres@cnLineLabelFontHeightF = 0.02
cnres@cnLineColor = "purple"
cnres@cnLineThicknessF = 4
cnres@cnInfoLabelOn = True
cnres@cnInfoLabelOrthogonalPosF = -0.65
cnres@cnInfoLabelFontHeightF = 0.02
cnres@cnLabelDrawOrder = "PostDraw"
cnres@gsnContourNegLineDashPattern = 2
cnres@gsnContourPosLineDashPattern = 0
cnres@gsnLeftString = ""
cnres@gsnRightString = ""
;====================================================================================

; res@gsnRightString = ""
res@gsnLeftString = "(a) CTRL: U200"
res@cnFillPalette = "WhiteBlueGreenYellowRed" ;
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = 30
res@cnMaxLevelValF = 60
res@cnLevelSpacingF = 10.
plot(0) = gsn_csm_contour_map(wks, u_ave(0,:,:), res)
plot_topo(0) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(0), plot_topo(0))

res@gsnLeftString = "(b) SEN-ALL: U200"
plot(1) = gsn_csm_contour_map(wks, u_ave(1,:,:), res)
plot_topo(1) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(1), plot_topo(1))

res@gsnLeftString = "(c) b-a"
res@cnFillPalette = "temp_diff_18lev"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -5
res@cnMaxLevelValF = 5
res@cnLevelSpacingF = 0.5
plot(2) = gsn_csm_contour_map(wks, diff_u, res)
plot_topo(2) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(2), plot_topo(2))

res@gsnLeftString = "(d) CTRL: V200"
res@cnFillPalette = "temp_diff_18lev" ;WhiteBlueGreenYellowRed
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -10 ;30
res@cnMaxLevelValF = 10 ;60
res@cnLevelSpacingF = 2. ;10.
plot(3) = gsn_csm_contour_map(wks, v_ave(0,:,:), res)
plot_topo(3) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(3), plot_topo(3))

res@gsnLeftString = "(e) SEN-ALL: V200"
plot(4) = gsn_csm_contour_map(wks, v_ave(1,:,:), res)
plot_topo(4) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(4), plot_topo(4))

res@gsnLeftString = "(f) e-d"
res@cnFillPalette = "temp_diff_18lev"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -3
res@cnMaxLevelValF = 3
res@cnLevelSpacingF = 0.3
plot(5) = gsn_csm_contour_map(wks, diff_v, res)
plot_topo(5) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(5), plot_topo(5))

res@gsnRightString = "gpm"
res@gsnLeftString = "(g) CTRL: H300"
res@cnFillPalette = "temp_diff_18lev"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -300
res@cnMaxLevelValF = 300
res@cnLevelSpacingF = 30.
plot(6) = gsn_csm_contour_map(wks, z_ave(0,:,:), res)
plot_topo(6) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(6), plot_topo(6))

cnres@cnLevelSelectionMode = "ManualLevels"
cnres@cnMinLevelValF = 5200
cnres@cnMaxLevelValF = 5800
cnres@cnLevelSpacingF = 100
plotc(6) = gsn_csm_contour(wks, h_ave(0,:,:), cnres)
overlay(plot(6), plotc(6))

res@gsnLeftString = "(h) SEN-ALL: H300"
plot(7) = gsn_csm_contour_map(wks, z_ave(1,:,:), res)
plot_topo(7) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(7), plot_topo(7))
plotc(7) = gsn_csm_contour(wks, h_ave(1,:,:), cnres)
overlay(plot(7), plotc(7))


res@gsnLeftString = "(i) h-g"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -50
res@cnMaxLevelValF = 50
res@cnLevelSpacingF = 5.
plot(8) = gsn_csm_contour_map(wks, diff_z, res)
plot_topo(8) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(8), plot_topo(8))

cnres@gsnContourNegLineDashPattern = 16
cnres@gsnContourPosLineDashPattern = 0
cnres@cnLevelSelectionMode = "ManualLevels"
cnres@cnMinLevelValF = -50
cnres@cnMaxLevelValF = 50
cnres@cnLevelSpacingF = 10
plotc(8) = gsn_csm_contour(wks, diff_h, cnres)
overlay(plot(8), plotc(8))


pres = True
pres@gsnMaximize = True
;pres@gsnPanelLabelBar = True
; pres@gsnPanelMainString = season+" mean heat flux: "+LEV+" hPa "

gsn_panel(wks, plot, (/3,3/), pres)

end