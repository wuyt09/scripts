
begin
lonstr = 75
lonend = 105
latstr = 25
latend = 40
yrStrt = 198102
yrLast = 200601


casename = (/"Ctrl_FAMIP","Hist_FAMIP","Hist_IO","Hist_PO","Hist_AO","Hist_PIO","Hist_PAO","Hist_IAO"/)

ncase = dimsizes(casename)

TS = new((/ncase,6,12,96,144/), float)

do i = 0, ncase-1
    diri = "/home/yangsong3/data-model/wuyt/TPbias/TPbias_"+casename(i)+"/a20191206/"
    fili1 = casename(i)+".cam.h0.SHFLX.1979-2005.nc"
    fili2 = casename(i)+".cam.h0.LHFLX.1979-2005.nc"
    fili3 = casename(i)+".cam.h0.FSNS.1979-2005.nc"
    fili4 = casename(i)+".cam.h0.FLNS.1979-2005.nc"
    fili5 = casename(i)+".clm2.h0.FSR.1979-2005.nc"

    f1 = addfile(diri+fili1, "r")
    f2 = addfile(diri+fili2, "r")
    f3 = addfile(diri+fili3, "r")
    f4 = addfile(diri+fili4, "r")
    f5 = addfile(diri+fili5, "r")
    TIME = f1->time
    YYYY = cd_calendar(TIME, -1)
    iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

    x1 = f1->SHFLX(iYYYY,:,:)
    x2 = f2->LHFLX(iYYYY,:,:)
    x3 = f3->FSNS(iYYYY,:,:)
    x4 = f4->FLNS(iYYYY,:,:)
    x5 = f5->FSR(iYYYY,:,:)

    TS(i,0,:,:,:) = clmMonTLL(x1)
    TS(i,1,:,:,:) = clmMonTLL(x2)
    TS(i,2,:,:,:) = clmMonTLL(x3)
    TS(i,3,:,:,:) = clmMonTLL(x4)
    TS(i,4,:,:,:) = clmMonTLL(x5)

    TSX = x1
    TSX = x1 + x2 - x3 + x4 + x5
    TS(i,5,:,:,:) = clmMonTLL(TSX)
end do
printVarSummary(TS)

ts = TS(0:ncase-2,:,:,:,:)
do i = 0, ncase-2
    ts(i,:,:,:,:) = TS(i+1,:,:,:,:) - TS(0,:,:,:,:)
end do
printVarSummary(ts)

;=========================================================
; outline of tp 
;=========================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)
TOPO = linint2_Wrap(topo&lon, topo&lat, topo, True, ts&lon, ts&lat, 0)

ts = mask(ts, TOPO.ge.1500, True)

tdiff = wgt_areaave_Wrap(ts(:,:,:,{latstr:latend},{lonstr:lonend}), 1, 1, 0)
printVarSummary(tdiff)

;=========================================================
; plot 
;=========================================================
Month = ispan(1, 12, 1)
month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/a20191206/"
figname = "12month_series_Ctrl_Basin_HeatFlux"

wks = gsn_open_wks("pdf", figpath+figname)

plot = new(8, graphic)

res = True
res@gsnDraw  = False
res@gsnFrame = False
res@gsnMaximize = True

res@vpWidthF = 0.7
res@vpHeightF = 0.4

res@trXMinF = 1
res@trXMaxF = 12
res@tmXBMode = "Explicit"
res@tmXBValues = Month
res@tmXBLabels = month
res@tmXTOn = False
;res@tiXAxisString = "Month"
res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025
res@gsnStringFontHeightF = 0.035

res@trYMinF = -10
res@trYMaxF = 10
; res@tmYLMode = "Explicit"
; res@tmYLValues = ispan(-, 0, 1)
; res@tmYLLabels = ispan(-5, 0, 1)
res@tmYROn = False
res@tiYAxisString = "W/m~S~2~N~"
res@gsnYRefLine = 0
res@gsnYRefLineDashPattern = 16
res@gsnYRefLineColor = "gray20"
res@gsnYRefLineThicknessF = 5

Colors = (/"red","blue","forestgreen","purple","goldenrod3","black","orange2"/)
;Patterns = (/0,0,0,0,16,16,16/)
res@xyLineThicknessF = 5
res@xyLineColors = Colors
res@xyDashPatterns = 0

leftsti = (/"Hist","IndianOnly","PacificOnly","AtlanticOnly",\
           "Pacific&Indian","Pacific&Atlantic","Indian&Atlantic"/)

;==================================================================================

do i = 0, ncase-2
    res@gsnLeftString = leftsti(i)+" - CTRL"
    plot(i) = gsn_csm_xy(wks, Month, tdiff(i,:,:), res)
end do

;==================================================================================
    res_lines                   = True                  ; polyline mods desired
    res_lines@gsLineThicknessF  = 5                     ; line thicker
    

    res_text                    = True                  ; text mods desired
    res_text@txFontHeightF      = 0.028                 ; change text size
    res_text@txJust             = "CenterLeft"          ; text justification

    Texts  = (/"SH","LH","NetS","NetR","G0","SUM"/)

; ;============================================================
    interval = 2
    gx = (/1.5,2.2/)
    text_pos = 2.4
    gy1 = (/-3,-3/)
    gy2 = gy1-interval ;(/15,15/)
    gy3 = gy2-interval
    gy4 = gy1 ;gy3-interval
    gy5 = gy2 ;gy4-interval
    gy6 = gy3 ;gy5-interval
    print(gy1)

    res_lines@gsLineColor = Colors(0)
    ;res_lines@gsLineDashPattern = Patterns(0)
    dum11 = gsn_add_polyline(wks, plot(0), gx, gy1, res_lines)
    dum12 = gsn_add_text(wks, plot(0), Texts(0), text_pos, gy1, res_text)
    res_lines@gsLineColor = Colors(1)
    ;res_lines@gsLineDashPattern = Patterns(1)
    dum21 = gsn_add_polyline(wks, plot(0), gx, gy2, res_lines)
    dum22 = gsn_add_text(wks, plot(0), Texts(1), text_pos, gy2, res_text)
    res_lines@gsLineColor = Colors(2)
    ;res_lines@gsLineDashPattern = Patterns(2)
    dum31 = gsn_add_polyline(wks, plot(0), gx, gy3, res_lines)
    dum32 = gsn_add_text(wks, plot(0), Texts(2), text_pos, gy3, res_text)

    text_pos = 4.4
    gx = (/3.5,4.2/) 
    res_lines@gsLineColor = Colors(3)
    ;res_lines@gsLineDashPattern = Patterns(3)
    dum41 = gsn_add_polyline(wks, plot(0), gx, gy4, res_lines)
    dum42 = gsn_add_text(wks, plot(0), Texts(3), text_pos, gy4, res_text) 
    res_lines@gsLineColor = Colors(4)
    ;res_lines@gsLineDashPattern = Patterns(4)
    dum51 = gsn_add_polyline(wks, plot(0), gx, gy5, res_lines)
    dum52 = gsn_add_text(wks, plot(0), Texts(4), text_pos, gy5, res_text)
    res_lines@gsLineColor = Colors(5)
    ;res_lines@gsLineDashPattern = Patterns(5)
    dum61 = gsn_add_polyline(wks, plot(0), gx, gy6, res_lines)
    dum62 = gsn_add_text(wks, plot(0), Texts(5), text_pos, gy6, res_text)

    ; dum1 = new(7, graphic)
    ; dum2 = dum1
    ; gy = fspan(-1.0, -1.6, 7)
    ; print(gy)
    ; do i = 0, 6
    ;     res_lines@gsLineColor = Colors(i)
    ;     dum1(i) = gsn_add_polyline(wks, plot(0), gx, (/gy(i),gy(i)/), res_lines)
    ;     dum2(i) = gsn_add_text(wks, plot(0), Texts(i), text_pos, (/gy(i),gy(i)/), res_text) ;
    ; end do


pres = True
pres@gsnMaximize = True
;pres@gsnPanelMainString = "Annual mean: All Flux"
;pres@pmLabelBarOrthogonalPosF = -0.02
;pres@lbLabelFontHeightF = 0.015
gsn_panel(wks, plot, (/2,4/), pres)


end
