
begin

lonstr = 75
lonend = 105
latstr = 25
latend = 40
yrStrt = 198102
yrLast = 200601


casename = (/"Ctrl_FAMIP","Hist_FAMIP"/) ;,"Hist_IO","Hist_PO","Hist_AO","Hist_PIO","Hist_PAO","Hist_IAO"

ncase = dimsizes(casename)
season = "DJF"

nyear = 25
t = new((/ncase,nyear,26,96,144/), float)

do i = 0, ncase-1
    diri = "/home/yangsong3/data-model/wuyt/TPbias/TPbias_"+casename(i)+"/a20191206/"
    fili = casename(i)+".cam.h0.T.NoInterpolation.1979-2005.nc"
    f = addfile(diri+fili, "r")

    TIME = f->time
    YYYY = cd_calendar(TIME, -1)
    iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

    T = f->T(iYYYY,::-1,:,:)
    t(i,:,:,:,:) = month_to_season(T, season)
end do
printVarSummary(t)

lev = t&lev
print(lev)


;=======================================================================
; ttest
;=======================================================================

t_ave = dim_avg_n_Wrap(t, 1)
t_var = dim_variance_n_Wrap(t, 1)

prob = ttest(t_ave(0,:,:,:), t_var(0,:,:,:), nyear, t_ave(1,:,:,:), t_var(1,:,:,:), nyear, False, False)
copy_VarCoords(t_ave(0,:,:,:), prob)
printVarSummary(prob)

diff = t_ave(0,:,:,:)
diff = t_ave(1,:,:,:) - t_ave(0,:,:,:)
printVarSummary(diff)

;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)

;=========================================================
; plot 
;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/a20200723/"
figname = season+"_T_multilevel_ttest"

; lonstr = 0
; lonend = 360
; latstr = -30
; latend = 60

wks = gsn_open_wks("pdf", figpath+figname)
plot = new(8, graphic)
plot_topo = plot
plot_prob = plot

 res_topo            = True
 res_topo@cnFillOn   = False
 res_topo@cnLinesOn  = True
 res_topo@gsnFrame   = False
 res_topo@gsnDraw    = False
 res_topo@cnLineLabelsOn         = False
 res_topo@cnLineThicknessF       = 3
 res_topo@cnLineColor            = "grey"
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
res@lbLabelBarOn          = False

res@mpMinLatF            = latstr      ; range to zoom in on
res@mpMaxLatF            = latend
res@mpMinLonF            = lonstr
res@mpMaxLonF            = lonend
; res@mpCenterLonF         = 180.
res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025

res@gsnStringFontHeightF = 0.035
;res@gsnRightString       = "K"

res@pmLabelBarOrthogonalPosF = 0.12
res@pmLabelBarWidthF = 0.7
res@lbLabelFontHeightF = 0.02

res@cnFillPalette  = "CBR_coldhot"

;----------------------------------------------
; draw prob
;----------------------------------------------
 res3  = True
 res3@gsnDraw   = False
 res3@gsnFrame  = False

 ;res3@cnFillOn  = False 
 res3@cnLinesOn = False
 res3@cnLineLabelsOn = False
 res3@cnInfoLabelOn  = False
 ;res3@lbLabelBarOn   = False

 res3@gsnRightString = ""

 res3@cnLevelSelectionMode = "ManualLevels"
 res3@cnMinLevelValF = 0.00
 res3@cnMaxLevelValF = 0.10
 res3@cnLevelSpacingF = 0.01

 ; res3@cnFillPattern = 17 ;(/17/)
 res3@cnFillColor  =  "red"
 res3@gsnLeftString = ""
 res3@cnFillDotSizeF = 0.003
 res3@cnFillScaleF  = 1.2
 ;----------------------------------------------
 ;----------------------------------------------

; do i = 0,3
; res@gsnLeftString  = sprintf("%4.2f", lev(3*i+1))+"hPa"
; res@gsnRightString = "Ctrl"
; res@cnLevelSelectionMode = "ExplicitLevels"
; res@cnMinLevelValF = 220
; res@cnMaxLevelValF = 280
; res@cnLevelSpacingF = 5.
; plot(i) = gsn_csm_contour_map(wks, t_ave(0,3*i+1,:,:), res)
; plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
; overlay(plot(i), plot_topo(i))

; res@gsnRightString = "Diff"
; res@cnLevelSelectionMode = "ManualLevels"
; res@cnMinLevelValF = -2
; res@cnMaxLevelValF = 2
; res@cnLevelSpacingF = 0.2
; plot(i+4) = gsn_csm_contour_map(wks, diff(3*i+1,:,:), res)
; plot_topo(i+4) = gsn_csm_contour(wks, topo, res_topo)
; overlay(plot(i+4), plot_topo(i+4))
; plot_prob(i) = gsn_csm_contour(wks,prob(i,:,:),res3)
; plot_prob(i) = ShadeLtContour(plot_prob(i),0.055,17)
; overlay(plot(i+4), plot_prob(i))
; end do

do i = 0, 7
    res@gsnLeftString  = sprintf("%4.2f", lev(i))+"hPa"
    res@gsnRightString = "Diff"
    res@cnLevelSelectionMode = "ManualLevels"
    res@cnMinLevelValF = -2
    res@cnMaxLevelValF = 2
    res@cnLevelSpacingF = 0.2
    plot(i) = gsn_csm_contour_map(wks, diff(i,:,:), res)
    plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(i), plot_topo(i))
    plot_prob(i) = gsn_csm_contour(wks,prob(i,:,:),res3)
    plot_prob(i) = ShadeLtContour(plot_prob(i),0.055,17)
    overlay(plot(i), plot_prob(i))
    
end do


pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar  = True
pres@gsnPanelMainString = season+" mean: T at diff levels"
gsn_panel(wks, plot, (/2,4/), pres)

end