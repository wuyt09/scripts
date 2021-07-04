begin

; lonstr = 75
; lonend = 105
; latstr = 25
; latend = 40
lonstr = -20
lonend = 160
latstr = 0
latend = 70

yrStrt = 1979
yrLast = 2017
year = ispan(yrStrt, yrLast, 1)
nyear = dimsizes(year)
season = "JJA"


fili = "/home/yangsong3/wuyt/sysu/scripts/monsoon_onset_pentad/LBQ_onset_pentad_1979_2017.txt"

x = asciiread(fili, (/nyear,4/), "integer")
print(x(nyear-1,:))

monsoon = (/"BOB","SCS","ISM"/)
lines = x(:,1:3)
;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)
printVarSummary(topo)

;===================================================================
diri = "/home/yangsong3/data-observation/ERA5-monthly/pressure/thermal_equation/"
fili = "q1.monmean.1979.2019.nc"
f = addfile(diri+fili, "r")

TIME = f->time
YYYY = cd_calendar(TIME, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

q1 = short2flt(f->q1(iYYYY,:,::-1,:))
q1 = q1*1004.7/9.8
printVarSummary(q1)

filips = "/home/yangsong3/data-observation/ERA5-monthly/surface/sur_pressure_mon_1x1_1979_2019.nc"
fps = addfile(filips, "r")
ps = short2flt(fps->sp(iYYYY,::-1,:))

plev = q1&lev
plev = plev*100
plev!0 = "lev"
plev&lev = q1&lev
plev@units = "Pa"

dp = dpres_plevel_Wrap(plev, ps, min(plev), 0)
printVarSummary(dp)

qdp = q1*dp
copy_VarCoords(q1, qdp)
printVarSummary(qdp)

Q1 = dim_sum_n_Wrap(qdp, 1)
copy_VarCoords(q1(:,0,:,:), Q1)
printVarSummary(Q1)

; t2m = mask(T2M, topo.ge.1500., True)
; TP_t2m = wgt_areaave_Wrap(t2m(:,{latstr:latend},{lonstr:lonend}), 1., 1., 0)

t2m_apr = Q1(3::12,:,:)
t2m_may = Q1(4::12,:,:)
t2m_jja = Q1(5::12,:,:) ;month_to_season(Q1, season)
printVarSummary(t2m_apr)

;========================================================================

rc_apr = regCoef_n(lines, t2m_apr, 0, 0)
rc_may = regCoef_n(lines, t2m_may, 0, 0)
rc_jja = regCoef_n(lines, t2m_jja, 0, 0)
copy_VarCoords(t2m_apr(0,:,:), rc_apr(0,:,:))
copy_VarCoords(t2m_may(0,:,:), rc_may(0,:,:))
copy_VarCoords(t2m_may(0,:,:), rc_jja(0,:,:))
printVarSummary(rc_apr)

;====================================================================
; plot
;====================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_surf_temp_climate_effects/20210702/"
figname = "reg_Apr_May_JJA_TP_Q1_ERA5_on_ISM_BOB_SCS_onset_time"
wks = gsn_open_wks("pdf", figpath+figname)
plot = new(9, graphic)
plot_topo = plot

;----------------------------------------------------
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
 res_topo@cnLevels               = (/2000/)
 res_topo@gsnLeftString = ""
 res_topo@gsnRightString = ""
;----------------------------------------------------
res = True
res@gsnDraw = False
res@gsnFrame = False
res@gsnMaximize = True

res@mpOutlineOn = True
res@mpFillOn = True
res@mpGeophysicalLineThicknessF = 4
res@mpGeophysicalLineColor = "grey25" ;"red" 
; res@mpFillDrawOrder = "PostDraw"
; res@mpCenterLonF = 180
res@mpMinLatF = latstr
res@mpMaxLatF = latend
res@mpMinLonF = lonstr
res@mpMaxLonF = lonend


res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False
; res@cnFillMode = "CellFill"

res@gsnStringFontHeightF = 0.025
res@tmXBLabelFontHeightF = 0.02
res@tmYLLabelFontHeightF = 0.02

res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 5
res@cnMaxLevelValF = 20
res@cnMinLevelValF = -20
res@cnFillPalette = "BlueWhiteOrangeRed"
res@lbLabelBarOn = False

do i = 0, 2
    res@gsnLeftString = monsoon(i)
    res@gsnRightString = "Apr Q1"
    plot(i) = gsn_csm_contour_map(wks, rc_apr(i,:,:), res)
    plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(i), plot_topo(i))
    res@gsnRightString = "May Q1"
    plot(i+3) = gsn_csm_contour_map(wks, rc_may(i,:,:), res)
    plot_topo(i+3) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(i+3), plot_topo(i+3))
    res@gsnRightString = "June Q1"
    plot(i+6) = gsn_csm_contour_map(wks, rc_jja(i,:,:), res)
    plot_topo(i+6) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(i+6), plot_topo(i+6))
end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
pres@pmLabelBarOrthogonalPosF = -0.012
pres@lbLabelFontHeightF = 0.012
gsn_panel(wks, plot, (/3,3/), pres)
end