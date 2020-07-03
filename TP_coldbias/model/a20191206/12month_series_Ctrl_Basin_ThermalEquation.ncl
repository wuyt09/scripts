
begin
lonstr = 75
lonend = 105
latstr = 25
latend = 40
yrStrt = 198102
yrLast = 200601

casename = (/"Ctrl_FAMIP","Hist_FAMIP","Hist_IO","Hist_PO","Hist_AO","Hist_PIO","Hist_PAO","Hist_IAO"/)

ncase = dimsizes(casename)
LEV = 500

TS = new((/ncase,4,12,96,144/), float)
diri = "/home/yangsong3/data-model/wuyt/TPbias/TEMP/"

do i = 0, ncase-1
    fili1 = casename(i)+".dTdt.1981-2005.mon.mean.nc"
    fili2 = casename(i)+".VgradT.1981-2005.mon.mean.nc"
    fili3 = casename(i)+".Vertical.1981-2005.mon.mean.nc"
    fili4 = casename(i)+".Q1.1981-2005.mon.mean.nc"

    f1 = addfile(diri+fili1, "r")
    f2 = addfile(diri+fili2, "r")
    f3 = addfile(diri+fili3, "r")
    f4 = addfile(diri+fili4, "r")

    x1 = f1->P1(:,{LEV},:,:)
    x2 = f2->P2(:,{LEV},:,:)
    x3 = f3->P3(:,{LEV},:,:)
    x4 = f4->Q1(:,{LEV},:,:)

    TS(i,0,:,:,:) = clmMonTLL(x1)
    TS(i,1,:,:,:) = clmMonTLL(x2)
    TS(i,2,:,:,:) = clmMonTLL(x3)
    TS(i,3,:,:,:) = clmMonTLL(x4)
end do
printVarSummary(TS)

TS(:,0:2,:,:,:) = TS(:,0:2,:,:,:)*86400
TS@units = "K/day"

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
figname = "12month_series_Ctrl_Basin_ThermalEquation"

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

res@trYMinF = -0.4
res@trYMaxF = 0.4
; res@tmYLMode = "Explicit"
; res@tmYLValues = ispan(-, 0, 1)
; res@tmYLLabels = ispan(-5, 0, 1)
res@tmYROn = False
res@tiYAxisString = "K/day"
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

    Texts  = (/"dTdt","VgradT","Vertical","Q1"/)

; ;============================================================
    interval = 0.1
    gx = (/1.5,2.2/)
    text_pos = 2.4
    gy1 = (/0.35,0.35/)
    gy2 = gy1-interval ;(/15,15/)
    gy3 = gy1
    gy4 = gy2 ;gy3-interval
    print(gy1)

    res_lines@gsLineColor = Colors(0)
    ;res_lines@gsLineDashPattern = Patterns(0)
    dum11 = gsn_add_polyline(wks, plot(0), gx, gy1, res_lines)
    dum12 = gsn_add_text(wks, plot(0), Texts(0), text_pos, gy1, res_text)
    res_lines@gsLineColor = Colors(1)
    ;res_lines@gsLineDashPattern = Patterns(1)
    dum21 = gsn_add_polyline(wks, plot(0), gx, gy2, res_lines)
    dum22 = gsn_add_text(wks, plot(0), Texts(1), text_pos, gy2, res_text)

    text_pos = 5.4
    gx = (/4.5,5.2/)
    res_lines@gsLineColor = Colors(2)
    ;res_lines@gsLineDashPattern = Patterns(2)
    dum31 = gsn_add_polyline(wks, plot(0), gx, gy3, res_lines)
    dum32 = gsn_add_text(wks, plot(0), Texts(2), text_pos, gy3, res_text)
    res_lines@gsLineColor = Colors(3)
    ;res_lines@gsLineDashPattern = Patterns(3)
    dum41 = gsn_add_polyline(wks, plot(0), gx, gy4, res_lines)
    dum42 = gsn_add_text(wks, plot(0), Texts(3), text_pos, gy4, res_text) 


pres = True
pres@gsnMaximize = True
pres@gsnPanelMainString = "500hPa TP area avg"

gsn_panel(wks, plot, (/2,4/), pres)


end
