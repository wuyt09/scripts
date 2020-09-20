
begin

lonstr = 60
lonend = 150
latstr = 15
latend = 50


yrStrt = 198102
yrLast = 200601


casename = (/"Ctrl_FAMIP","Hist_FAMIP"/) ;,"Hist_IO","Hist_PO","Hist_AO","Hist_PIO","Hist_PAO","Hist_IAO"
ncase = dimsizes(casename)

season = "DJF"
nlat = 96
nlon = 144

LEV = 300

VT_bar = new((/ncase,nlat,nlon/), float)
VT_bar!0 = "case"
VT_bar&case = casename

VT_ano = VT_bar
VT_tot = VT_bar

do i = 0, ncase-1
    diri = "/home/yangsong3/data-model/wuyt/TPbias/TPbias_"+casename(i)+"/a20191206/"

    fili1 = casename(i)+".cam.h0.T.1979-2005.nc"
    fili2 = casename(i)+".cam.h0.U.1979-2005.nc"
    f1 = addfile(diri+fili1, "r")
    f2 = addfile(diri+fili2, "r")

    TIME = f1->time
    YYYY = cd_calendar(TIME, -1)
    iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

    T = f1->T(iYYYY,{LEV},:,:)
    V = f2->U(iYYYY,{LEV},:,:)

    T_season = month_to_season(T, season)
    V_season = month_to_season(V, season)

    T_bar = dim_avg_n_Wrap(T_season, 0)
    V_bar = dim_avg_n_Wrap(V_season, 0)
    VT_bar(i,:,:) = V_bar*T_bar
    copy_VarCoords(V_bar, VT_bar(i,:,:))

    T_ano = T_season - conform_dims(dimsizes(T_season), T_bar, (/1,2/))
    V_ano = V_season - conform_dims(dimsizes(V_season), V_bar, (/1,2/))
    VT_ano_med = T_ano*V_ano
    VT_ano(i,:,:) = (/dim_avg_n_Wrap(VT_ano_med, 0)/)
    copy_VarCoords(V_bar, VT_ano(i,:,:))

    VT_tot(i,:,:) = VT_bar(i,:,:)+VT_ano(i,:,:)
    copy_VarCoords(V_bar, VT_tot(i,:,:))

end do

printVarSummary(VT_bar)
printVarSummary(VT_ano)
printVarSummary(VT_tot)

diff_bar = VT_bar(0,:,:)
diff_ano = VT_ano(0,:,:)
diff_tot = VT_tot(0,:,:)
diff_bar = VT_bar(1,:,:) - VT_bar(0,:,:)
diff_ano = VT_ano(1,:,:) - VT_ano(0,:,:)
diff_tot = VT_tot(1,:,:) - VT_tot(0,:,:)

;=========================================================
; outline of tp
;=========================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)

;=========================================================
; plot
;=========================================================

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/a20200910/"
figname = season+"_Ctrl_diff_heat_flux_UT"


wks = gsn_open_wks("pdf", figpath+figname)
plot = new(6, graphic)
plot_topo = plot

 res_topo            = True
 res_topo@cnFillOn   = False
 res_topo@cnLinesOn  = True
 res_topo@gsnFrame   = False
 res_topo@gsnDraw    = False
 res_topo@cnLineLabelsOn         = False
 res_topo@cnLineThicknessF       = 3
 res_topo@cnLineColor            = "purple"
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
;res@lbLabelBarOn          = False
res@pmLabelBarOrthogonalPosF = 0.2

res@mpMinLatF            = latstr      ; range to zoom in on
res@mpMaxLatF            = latend
res@mpMinLonF            = lonstr
res@mpMaxLonF            = lonend
; res@mpCenterLonF         = 180.
res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025

res@gsnStringFontHeightF = 0.035


res@gsnRightString = "Ctrl"
res@gsnLeftString = "(a) <U>*<T>"
res@cnFillPalette = "WhiteBlueGreenYellowRed"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = 0
res@cnMaxLevelValF = 20000
res@cnLevelSpacingF = 1000.
plot(0) = gsn_csm_contour_map(wks, VT_bar(0,:,:), res)
plot_topo(0) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(0), plot_topo(0))

res@gsnLeftString = "(b) <U'*T'>"
res@cnFillPalette = "temp_diff_18lev"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -5
res@cnMaxLevelValF = 5
res@cnLevelSpacingF = 0.5
plot(1) = gsn_csm_contour_map(wks, VT_ano(0,:,:), res)
plot_topo(1) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(1), plot_topo(1))

res@gsnLeftString = "(c) <U*T>: a+b"
res@cnFillPalette = "WhiteBlueGreenYellowRed"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = 0
res@cnMaxLevelValF = 20000
res@cnLevelSpacingF = 1000.
plot(2) = gsn_csm_contour_map(wks, VT_tot(0,:,:), res)
plot_topo(2) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(2), plot_topo(2))

res@gsnRightString = "Diff"
res@gsnLeftString = "(d) <U>*<T>"
res@cnFillPalette = "temp_diff_18lev"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -1000
res@cnMaxLevelValF = 1000
res@cnLevelSpacingF = 100.
plot(3) = gsn_csm_contour_map(wks, diff_bar, res)
plot_topo(3) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(3), plot_topo(3))

res@gsnLeftString = "(e) <U'*T'>"
res@cnFillPalette = "temp_diff_18lev"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -5
res@cnMaxLevelValF = 5
res@cnLevelSpacingF = 0.5
plot(4) = gsn_csm_contour_map(wks, diff_ano, res)
plot_topo(4) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(4), plot_topo(4))

res@gsnLeftString = "(f) <U*T>: d+e"
res@cnFillPalette = "temp_diff_18lev"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -1000
res@cnMaxLevelValF = 1000
res@cnLevelSpacingF = 100.
plot(5) = gsn_csm_contour_map(wks, diff_bar, res)
plot_topo(5) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(5), plot_topo(5))

pres = True
pres@gsnMaximize = True
;pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = season+" mean heat flux: "+LEV+" hPa "

gsn_panel(wks, plot, (/2,3/), pres)

end