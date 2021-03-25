
begin

lonstr = 60
lonend = 150
latstr = 15
latend = 50


yrStrt = 198102
yrLast = 200601


casename = (/"Ctrl_FAMIP","Hist_FAMIP"/) ;,"Hist_PO","Hist_FAMIP_CMIP6","Hist_IO","Hist_AO","Hist_PIO","Hist_PAO","Hist_IAO"
ncase = dimsizes(casename)

season = "DJF"
nlat = 96
nlon = 144

LEV = 300

V_season = new((/ncase,nlat,nlon/), float)
V_season!0 = "case"
V_season&case = casename
U_season = V_season
T_season = V_season

do i = 0, ncase-1
    diri = "/home/yangsong3/data-model/wuyt/TPbias/TPbias_"+casename(i)+"/a20191206/"

    fili1 = casename(i)+".cam.h0.T.1979-2005.nc"
    fili2 = casename(i)+".cam.h0.V.1979-2005.nc"    
    fili3 = casename(i)+".cam.h0.U.1979-2005.nc"
    f1 = addfile(diri+fili1, "r")
    f2 = addfile(diri+fili2, "r")
    f3 = addfile(diri+fili3, "r")

    TIME = f1->time
    YYYY = cd_calendar(TIME, -1)
    iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

    T = f1->T(iYYYY,{LEV},:,:)
    V = f2->V(iYYYY,{LEV},:,:)
    U = f3->U(iYYYY,{LEV},:,:)

    U_season_med = month_to_season(U, season)
    V_season_med = month_to_season(V, season)
    T_season_med = month_to_season(T, season)

    U_season(i,:,:) = dim_avg_n_Wrap(U_season_med, 0)
    V_season(i,:,:) = dim_avg_n_Wrap(V_season_med, 0)
    T_season(i,:,:) = dim_avg_n_Wrap(T_season_med, 0)

end do
printVarSummary(T_season)

diff_v = V_season(0:ncase-2,:,:)
diff_u = U_season(0:ncase-2,:,:)
diff_t = T_season(0:ncase-2,:,:)
mean_v = V_season(0:ncase-2,:,:)
mean_u = U_season(0:ncase-2,:,:)
mean_t = T_season(0:ncase-2,:,:)
do i = 0, ncase-2
    diff_v(i,:,:) = V_season(i+1,:,:) - V_season(0,:,:)
    diff_u(i,:,:) = U_season(i+1,:,:) - U_season(0,:,:)
    diff_t(i,:,:) = T_season(i+1,:,:) - T_season(0,:,:)
    mean_v(i,:,:) = (V_season(i+1,:,:) + V_season(0,:,:))/2.
    mean_u(i,:,:) = (U_season(i+1,:,:) + U_season(0,:,:))/2.
    mean_t(i,:,:) = (T_season(i+1,:,:) + T_season(0,:,:))/2.
end do

term = new((/4,nlat,nlon/), float)
term(0,:,:) = diff_u*mean_t
term(1,:,:) = diff_t*mean_u
term(2,:,:) = diff_v*mean_t
term(3,:,:) = diff_t*mean_v
copy_VarCoords(T_season(0,:,:), term(0,:,:))
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
figname = "Fig9_v2"


title = (/"(a) deltaU*avgT","(b) deltaT*avgU","(c) deltaV*avgT","(d) deltaT*avgV"/)

wks = gsn_open_wks("pdf", figpath+figname)
plot = new(4, graphic)
plotc = plot
plot_topo = plot
plotc = plot

 res_topo            = True
 res_topo@cnFillOn   = False
 res_topo@cnLinesOn  = True
 res_topo@gsnFrame   = False
 res_topo@gsnDraw    = False
 res_topo@cnLineLabelsOn         = False
 res_topo@cnLineThicknessF       = 3
 res_topo@cnLineColor            = "grey40"
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
; res@pmLabelBarOrthogonalPosF = 0.2

res@mpMinLatF            = latstr      ; range to zoom in on
res@mpMaxLatF            = latend
res@mpMinLonF            = lonstr
res@mpMaxLonF            = lonend
; res@mpCenterLonF         = 180.
res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025

res@gsnStringFontHeightF = 0.035


res@cnFillPalette = "temp_diff_18lev"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -1000
res@cnMaxLevelValF = 1000
res@cnLevelSpacingF = 100
; res@gsnRightString = "*10~S~-3~N~"

res@gsnLeftString = title(0)
plot(0) = gsn_csm_contour_map(wks, term(0,:,:), res)
plot_topo(0) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(0), plot_topo(0))

res@cnMinLevelValF = -50
res@cnMaxLevelValF = 50
res@cnLevelSpacingF = 5
res@gsnLeftString = title(1)
plot(1) = gsn_csm_contour_map(wks, term(1,:,:), res)
plot_topo(1) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(1), plot_topo(1))

res@cnMinLevelValF = -500
res@cnMaxLevelValF = 500
res@cnLevelSpacingF = 50
res@gsnLeftString = title(2)
plot(2) = gsn_csm_contour_map(wks, term(2,:,:), res)
plot_topo(2) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(2), plot_topo(2))

res@cnMinLevelValF = -10
res@cnMaxLevelValF = 10
res@cnLevelSpacingF = 1
res@gsnLeftString = title(3)
plot(3) = gsn_csm_contour_map(wks, term(3,:,:), res)
plot_topo(3) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(3), plot_topo(3))


pres = True
pres@gsnMaximize = True
; pres@gsnPanelLabelBar = True
; pres@pmLabelBarOrthogonalPosF = -0.02
; pres@lbLabelFontHeightF = 0.015
; pres@gsnPanelMainString = season+" mean: "+LEV+" hPa "

gsn_panel(wks, plot, (/2,2/), pres)

; pres = True
; ; pres@gsnMaximize = True
; pres@gsnDraw  = False
; pres@gsnFrame = False

; height = 0.2
; pres@gsnDraw = True
; pres@gsnPanelLabelBar = True
; pres@lbLabelFontHeightF = 0.012
; pres@gsnPanelTop = 0.7
; pres@gsnPanelBottom = pres@gsnPanelTop - height
; gsn_panel(wks, plot(0:2), (/1,3/), pres)


; pres@gsnPanelLabelBar = True
; pres@gsnPanelTop = pres@gsnPanelBottom
; pres@gsnPanelBottom = pres@gsnPanelTop - height
; pres@gsnPanelLabelBar = True
; pres@lbLabelFontHeightF = 0.012
; gsn_panel(wks, plot(3:5), (/1,3/), pres)

end