
begin
lonstr = 75 ; 72 ;75
lonend = 105 ;136 ;105
latstr = 25 ;18 ;25
latend = 40 ;54 ; 40
yrStrt = 198102
yrLast = 200601

casename = (/"Ctrl_FAMIP","Hist_FAMIP","Hist_FAMIP_CMIP6","Hist_IO","Hist_PO","Hist_AO","Hist_PIO","Hist_PAO","Hist_IAO"/)

ncase = dimsizes(casename)

;U = new((/ncase,324,96,144/), float)
TS = new((/ncase,12,96,144/), float)

do i = 0, ncase-1
    diri = "/home/yangsong3/data-model/wuyt/TPbias/TPbias_"+casename(i)+"/a20191206/"
    fili = casename(i)+".cam.h0.TREFHT.1979-2005.nc"
    f = addfile(diri+fili, "r")
    TIME = f->time
    YYYY = cd_calendar(TIME, -1)
    iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
    TSX = f->TREFHT ;(iYYYY,:,:)
    TS(i,:,:,:) = clmMonTLL(TSX)
delete([/diri,fili,f,TSX/])
end do
printVarSummary(TS)

;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)
TOPO = linint2_Wrap(topo&lon, topo&lat, topo, True, TS&lon, TS&lat, 0)
TS = mask(TS, TOPO.ge.1500, True)

ts = TS(0:ncase-2,:,:,:)
do i = 0, ncase-2
    ts(i,:,:,:) = TS(i+1,:,:,:) - TS(0,:,:,:)
end do
printVarSummary(ts)

tdiff = wgt_areaave_Wrap(ts(:,:,{25:40},{75:105}), 1, 1, 0)
printVarSummary(tdiff)
t_ave = (tdiff(0,0) + tdiff(0,1) + tdiff(0,11))/3
print(t_ave)
;===================================================================
; Plot
;===================================================================

Month = ispan(1, 12, 1)
month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/a20200910/"
figname = "12month_series_Ctrl_Basin_TAS"
wks = gsn_open_wks("pdf", figpath+figname)

plot = new(1, graphic)

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

res@trYMinF = -2
res@trYMaxF = 1
; res@tmYLMode = "Explicit"
; res@tmYLValues = ispan(-, 0, 1)
; res@tmYLLabels = ispan(-5, 0, 1)
res@tmYROn = False
res@tiYAxisString = "Temperature (~S~o~N~C)"
res@gsnYRefLine = 0
res@gsnYRefLineDashPattern = 16
res@gsnYRefLineColor = "gray20"
res@gsnYRefLineThicknessF = 5

Colors = (/"black","red","blue","forestgreen","purple","chocolate1","goldenrod3","lightblue"/)
Patterns = (/0,0,0,0,0,16,16,16/)
Thick = (/10,10,4,10,4,10,10,4/)
res@xyLineThicknesses = Thick
res@xyLineColors = Colors
res@xyDashPatterns = Patterns

res@gsnLeftString = "Diff with Ctrl run"
res@gsnStringFontHeightF = 0.025

plot(0) = gsn_csm_xy(wks, Month, tdiff, res)


;==================================================================================
    res_lines                   = True                  ; polyline mods desired
    res_lines@gsLineThicknessF  = 5                     ; line thicker
    

    res_text                    = True                  ; text mods desired
    res_text@txFontHeightF      = 0.018                 ; change text size
    res_text@txJust             = "CenterLeft"          ; text justification

    Texts  = (/"All Basin","All Basin CMIP6","IndianOnly","PacificOnly","AtlanticOnly",\
           "Pacific&Indian","Pacific&Atlantic","Indian&Atlantic"/)

; ;============================================================
    interval = 0.2
    gx = (/1.5,2.2/)
    text_pos = 2.3
    gy1 = (/0.8,0.8/)
    gy2 = gy1-interval ;(/15,15/)
    gy3 = gy2-interval
    gy4 = gy3-interval
    gy5 = gy1 ;gy4-interval
    gy6 = gy2 ;gy5-interval
    gy7 = gy3 ;gy6-interval
    gy8 = gy4
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
    res_lines@gsLineColor = Colors(2)
    res_lines@gsLineDashPattern = Patterns(2)
    res_lines@xyLineThicknessF = Thick(2)
    dum31 = gsn_add_polyline(wks, plot(0), gx, gy3, res_lines)
    dum32 = gsn_add_text(wks, plot(0), Texts(2), text_pos, gy3, res_text)
    res_lines@gsLineColor = Colors(3)
    res_lines@gsLineDashPattern = Patterns(3)
    res_lines@xyLineThicknessF = Thick(3)
    dum41 = gsn_add_polyline(wks, plot(0), gx, gy4, res_lines)
    dum42 = gsn_add_text(wks, plot(0), Texts(3), text_pos, gy4, res_text)

    text_pos = 5.8
    gx = (/5,5.7/)    
    res_lines@gsLineColor = Colors(4)
    res_lines@gsLineDashPattern = Patterns(4)
    res_lines@xyLineThicknessF = Thick(4)
    dum51 = gsn_add_polyline(wks, plot(0), gx, gy5, res_lines)
    dum52 = gsn_add_text(wks, plot(0), Texts(4), text_pos, gy5, res_text)
    res_lines@gsLineColor = Colors(5)
    res_lines@gsLineDashPattern = Patterns(5)
    res_lines@xyLineThicknessF = Thick(5)
    dum61 = gsn_add_polyline(wks, plot(0), gx, gy6, res_lines)
    dum62 = gsn_add_text(wks, plot(0), Texts(5), text_pos, gy6, res_text)
    res_lines@gsLineColor = Colors(6)
    res_lines@gsLineDashPattern = Patterns(6)
    res_lines@xyLineThicknessF = Thick(6)
    dum71 = gsn_add_polyline(wks, plot(0), gx, gy7, res_lines)
    dum72 = gsn_add_text(wks, plot(0), Texts(6), text_pos, gy7, res_text)
    res_lines@gsLineColor = Colors(7)
    res_lines@gsLineDashPattern = Patterns(7)
    res_lines@xyLineThicknessF = Thick(7)
    dum81 = gsn_add_polyline(wks, plot(0), gx, gy8, res_lines)
    dum82 = gsn_add_text(wks, plot(0), Texts(7), text_pos, gy8, res_text)

    ; dum1 = new(7, graphic)
    ; dum2 = dum1
    ; gy = fspan(-1.0, -1.6, 7)
    ; print(gy)
    ; do i = 0, 6
    ;     res_lines@gsLineColor = Colors(i)
    ;     dum1(i) = gsn_add_polyline(wks, plot(0), gx, (/gy(i),gy(i)/), res_lines)
    ;     dum2(i) = gsn_add_text(wks, plot(0), Texts(i), text_pos, (/gy(i),gy(i)/), res_text) ;
    ; end do


draw(plot(0))
;draw(plot(1))
frame(wks)

end