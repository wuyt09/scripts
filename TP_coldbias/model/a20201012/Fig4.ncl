
begin

lonstr = 75
lonend = 105
latstr = 25
latend = 40
yrStrt = 198102
yrLast = 200601

;,"Hist_PIO","Hist_PAO","Hist_IAO"
casename = (/"Ctrl_FAMIP","Hist_FAMIP","Hist_FAMIP_CMIP6","Hist_IO","Hist_PO","Hist_AO"/)
ncase = dimsizes(casename)

nyear = 25
t = new((/ncase,nyear,12,96,144/), float)

do i = 0, ncase-1
    diri = "/home/yangsong3/data-model/wuyt/TPbias/TPbias_"+casename(i)+"/a20191206/"
    fili1 = casename(i)+".cam.h0.TREFHT.1979-2005.nc"

    f1 = addfile(diri+fili1, "r")

    TIME = f1->time
    YYYY = cd_calendar(TIME, -1)
    iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

    T1 = f1->TREFHT(iYYYY,:,:)
    t(i,:,:,:,:) = reshape(T1, (/nyear,12,96,144/))
    copy_VarCoords(T1(0,:,:), t(i,0,0,:,:))
    ; t(i,:,:,:) = clmMonTLL(T1)
end do
printVarSummary(t)

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
; printMinMax(OMEGA(0,:,:,:), False)
; printMinMax(diff, False)
; exit
; diff = mask(diff, prob.le.0.1, True)
printVarSummary(diff)

Month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
Num = (/"(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)","(j)","(k)","(l)"/)

;=============================================================================
; annul cycle
;
;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)
topo = smth9_Wrap(topo, 0.5, 0.25, True)

TS = t_ave
TOPO = linint2_Wrap(topo&lon, topo&lat, topo, True, TS&lon, TS&lat, 0)
TS = mask(TS, TOPO.ge.1500, True)
ts = TS(0:ncase-2,:,:,:)
do i = 0, ncase-2
    ts(i,:,:,:) = TS(i+1,:,:,:) - TS(0,:,:,:)
end do
printVarSummary(ts)

tdiff = wgt_areaave_Wrap(ts(:,:,{25:40},{75:105}), 1, 1, 0)
printVarSummary(tdiff)

ave5_w  = (tdiff(:,0)+tdiff(:,1)+tdiff(:,11))/3

ave5_s  = (tdiff(:,5)+tdiff(:,6)+tdiff(:,7))/3
print(sprintf("%4.2f", ave5_w))
print(sprintf("%4.2f", ave5_s))

;=========================================================
; plot 
;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/a20201012/"
figname = "Fig4"

pwks = "png"
pwks@wkHeight              = 3048
pwks@wkWidth               = 6048

wks = gsn_open_wks("pdf", figpath+figname) ;
plot = new(12, graphic)
plot_topo = plot
plot_prob = plot

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

res@vpXF = 0.1
res@vpYF = 0.8

res@gsnMaximize           = True     ; Make plot fit the frame
res@cnFillOn              = True     ; turn on color fill
res@cnLinesOn             = False    ; turn of contour lines
res@cnLineLabelsOn        = False
res@cnInfoLabelOn         = False
;res@lbLabelBarOn          = False

res@mpMinLatF            = latstr      ; range to zoom in on
res@mpMaxLatF            = latend
res@mpMinLonF            = lonstr
res@mpMaxLonF            = lonend
res@mpCenterLonF         = 180.
res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025

res@gsnStringFontHeightF = 0.045
res@gsnRightString       = "K"

res@pmLabelBarOrthogonalPosF = 0.12
res@pmLabelBarWidthF = 0.7
res@lbLabelFontHeightF = 0.02
res@cnFillPalette = "CBR_coldhot"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMaxLevelValF = 3
res@cnMinLevelValF = -3
res@cnLevelSpacingF = 0.2
res@lbLabelBarOn = False

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
 res3@cnFillColor  =  "black"
 res3@gsnLeftString = ""
 res3@cnFillDotSizeF = 0.003
 res3@cnFillScaleF  = 1.2
 ;----------------------------------------------
 ;----------------------------------------------
do i = 0, 11
    res@gsnLeftString = Num(i)+" "+Month(i)
    plot(i) = gsn_csm_contour_map(wks, diff(i,:,:), res)
    plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(i), plot_topo(i))
    plot_prob(i) = gsn_csm_contour(wks,prob(i,:,:),res3)
    plot_prob(i) = ShadeLtContour(plot_prob(i),0.055,6)
    overlay(plot(i), plot_prob(i))
end do

; do i = 10, 11
;     res@gsnLeftString = Num(i)+" "+Month(i)
;     plot(i+2) = gsn_csm_contour_map(wks, diff(i,:,:), res)
;     plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
;     overlay(plot(i+2), plot_topo(i))
;     plot_prob(i) = gsn_csm_contour(wks,prob(i,:,:),res3)
;     plot_prob(i) = ShadeLtContour(plot_prob(i),0.055,17)
;     overlay(plot(i+2), plot_prob(i))
; end do

pres = True
; pres@gsnDraw = False
pres@gsnFrame = False
pres@gsnPanelLabelBar = True
pres@lbLabelFontHeightF = 0.012
pres@gsnPanelRowSpec = True
pres@gsnPanelCenter = False
gsn_panel(wks, plot, (/4,4,2,2/), pres)

;================================================================
; annual cycle
;================================================================
delete(plot)
plot = new(1, graphic)
delete(res)
res = True
res@gsnDraw  = False
res@gsnFrame = False
; res@gsnMaximize = True

res@vpWidthF = 0.38
res@vpHeightF = 0.2

res@vpXF = 0.6
res@vpYF = 0.48

res@trXMinF = 1
res@trXMaxF = 12
res@tmXBMode = "Explicit"
res@tmXBValues = ispan(1, 12, 1)
res@tmXBLabels = Month
res@tmXTOn = False
;res@tiXAxisString = "Month"

res@trYMinF = -2
res@trYMaxF = 1
; res@tmYLMode = "Explicit"
; res@tmYLValues = ispan(-, 0, 1)
; res@tmYLLabels = ispan(-5, 0, 1)
res@tmYROn = False
res@tiYAxisString = "Temperature (K)" ;(~S~o~N~C)
res@gsnYRefLine = 0
res@gsnYRefLineDashPattern = 16
res@gsnYRefLineColor = "gray20"
res@gsnYRefLineThicknessF = 5

Colors = (/"black","red","blue","forestgreen","purple","chocolate1","goldenrod3","lightblue"/)
Patterns = (/0,0,0,0,0,16,16,16/)
Thick = (/12,12,4,12,4,12,12,4/)
res@xyLineThicknesses = Thick
res@xyLineColors = Colors
res@xyDashPatterns = Patterns

res@gsnLeftString = "(m) Areal averaged TP surf temp change"
res@gsnStringFontHeightF = 0.012
res@gsnRightString = "K"

plot(0) = gsn_csm_xy(wks, ispan(1, 12, 1), tdiff, res)
;==================================================================================
    res_lines                   = True                  ; polyline mods desired
    res_lines@gsLineThicknessF  = 5                     ; line thicker
    

    res_text                    = True                  ; text mods desired
    res_text@txFontHeightF      = 0.008                 ; change text size
    res_text@txJust             = "CenterLeft"          ; text justification

    Texts  = (/"SEN_ALL","SEN_ALL6","SEN_IND","SEN_PAC","SEN_ATL",\
           "Pacific&Indian","Pacific&Atlantic","Indian&Atlantic"/)

; ;============================================================
    interval = 0.2
    gx = (/1.5,2.2/)
    text_pos = 2.3
    gy1 = (/0.8,0.8/)
    gy2 = gy1 - interval ;(/15,15/)
    gy3 = gy1
    gy4 = gy2
    gy5 = gy4 - interval
    print(gy1)

    res_lines@gsLineColor = Colors(0)
    res_lines@gsLineDashPattern = Patterns(0)
    res_lines@xyLineThicknessF = Thick(0)
    dum11 = gsn_add_polyline(wks, plot(0), gx, gy1, res_lines)
    dum12 = gsn_add_text(wks, plot(0), Texts(0), text_pos, gy1, res_text)
    res_lines@gsLineColor = Colors(1)
    res_lines@gsLineDashPattern = Patterns(1)
    res_lines@xyLineThicknessF = Thick(1)
    dum21 = gsn_add_polyline(wks, plot(0), gx, gy2, res_lines)
    dum22 = gsn_add_text(wks, plot(0), Texts(1), text_pos, gy2, res_text)


    text_pos = 6.0
    gx = (/5.2,5.9/)    
    res_lines@gsLineColor = Colors(3)
    res_lines@gsLineDashPattern = Patterns(3)
    res_lines@xyLineThicknessF = Thick(3)
    dum31 = gsn_add_polyline(wks, plot(0), gx, gy3, res_lines)
    dum32 = gsn_add_text(wks, plot(0), Texts(3), text_pos, gy3, res_text)
    res_lines@gsLineColor = Colors(2)
    res_lines@gsLineDashPattern = Patterns(2)
    res_lines@xyLineThicknessF = Thick(2)
    dum41 = gsn_add_polyline(wks, plot(0), gx, gy4, res_lines)
    dum42 = gsn_add_text(wks, plot(0), Texts(2), text_pos, gy4, res_text)

    res_lines@gsLineColor = Colors(4)
    res_lines@gsLineDashPattern = Patterns(4)
    res_lines@xyLineThicknessF = Thick(4)
    dum51 = gsn_add_polyline(wks, plot(0), gx, gy5, res_lines)
    dum52 = gsn_add_text(wks, plot(0), Texts(4), text_pos, gy5, res_text)

; frame(wks)
draw(plot(0))

end