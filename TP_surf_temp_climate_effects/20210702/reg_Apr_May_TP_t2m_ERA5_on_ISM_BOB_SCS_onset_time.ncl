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
diri = "/home/yangsong3/data-observation/ERA5-monthly/surface/"
fili = "t2m_mon_1x1_1979_2019.nc"
f = addfile(diri+fili, "r")

TIME = f->time
YYYY = cd_calendar(TIME, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

T2M = short2flt(f->t2m(iYYYY,::-1,:))
printVarSummary(T2M)

; t2m = mask(T2M, topo.ge.1500., True)
; TP_t2m = wgt_areaave_Wrap(t2m(:,{latstr:latend},{lonstr:lonend}), 1., 1., 0)

t2m_apr = T2M(3::12,:,:)
t2m_may = T2M(4::12,:,:)
printVarSummary(t2m_apr)

;========================================================================

rc_apr = regCoef_n(lines, t2m_apr, 0, 0)
rc_may = regCoef_n(lines, t2m_may, 0, 0)
copy_VarCoords(t2m_apr(0,:,:), rc_apr(0,:,:))
copy_VarCoords(t2m_may(0,:,:), rc_may(0,:,:))
printVarSummary(rc_apr)

;====================================================================
; plot
;====================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_surf_temp_climate_effects/20210702/"
figname = "reg_Apr_May_TP_t2m_ERA5_on_ISM_BOB_SCS_onset_time"
wks = gsn_open_wks("pdf", figpath+figname)
plot = new(6, graphic)
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
res@cnLevelSpacingF = 0.1
res@cnMaxLevelValF = 0.5
res@cnMinLevelValF = -0.5
res@cnFillPalette = "BlueWhiteOrangeRed"
res@lbLabelBarOn = False

do i = 0, 2
    res@gsnLeftString = monsoon(i)
    res@gsnRightString = "Apr t2m"
    plot(i) = gsn_csm_contour_map(wks, rc_apr(i,:,:), res)
    plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(i), plot_topo(i))
    res@gsnRightString = "May t2m"
    plot(i+3) = gsn_csm_contour_map(wks, rc_may(i,:,:), res)
    plot_topo(i+3) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(i+3), plot_topo(i+3))
end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
pres@pmLabelBarOrthogonalPosF = -0.012
pres@lbLabelFontHeightF = 0.012
gsn_panel(wks, plot, (/2,3/), pres)
end