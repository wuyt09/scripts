
begin

yrStrt = 1979
yrLast = 2018
mnStrt = 1
mnLast = 12

year = ispan(yrStrt, yrLast, 1)
month = ispan(mnStrt, mnLast, 1)
penta = ispan(1, 73, 1)
nyear = dimsizes(year)
nmonth = dimsizes(month)
npenta = dimsizes(penta)

nlat = 181
nlon = 360
nlev = 19


diri = "/home/yangsong3/wuyt/sysu/scripts/monsoon_onset_pentad/ISM/"
fili1 = "MTG_onset_pentad_1979_2018.txt"
fili2 = "LBQ_onset_pentad_1979_2017.txt"
fili3 = "IMD_onset_pentad_1970_2020.txt"
fili4 = "TTGI_onset_pentad_1979_2018.txt"
fili5 = "WFI_onset_pentad_1979_2018.txt"
fili6 = "HOWI_onset_pentad_1979_2018.txt"

x1 = asciiread(diri+fili1, (/nyear,2/), "integer")
x2 = asciiread(diri+fili2, (/nyear-1,2/), "integer")
x3 = asciiread(diri+fili3, (/51,3/), "integer")
x4 = asciiread(diri+fili4, (/nyear,2/), "integer")
x5 = asciiread(diri+fili5, (/nyear,2/), "integer")
x6 = asciiread(diri+fili6, (/nyear,2/), "integer")
; print(x3(:,2))

x = new((/6,nyear/), integer)
x(0,:) = x1(:,1)
x(1,0:nyear-2) = x2(:,1)
x(2,:) = x3(9:48,2)
x(3,:) = x4(:,1)
x(4,:) = x5(:,1)
x(5,:) = x6(:,1)
printVarSummary(x)

corr = new(6, float)
do i = 0, 5
    corr(i) = escorc(x(3,:), x(i,:))
end do
print(sprintf("%4.2f", corr))
;=============================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/monsoon_onset_pentad/ISM/"
figname = "3_multi_onset_time_corr"

Colors = (/"black","red","blue","goldenrod2","green4","purple"/)
Patterns = (/0,0,0,0,0,0/)
Thick = (/20,15,12,10,7,5/)

wks = gsn_open_wks("pdf", figpath+figname)
plot = new(1, graphic)

res = True
res@gsnDraw  = False
res@gsnFrame = False
res@gsnMaximize = True

res@vpWidthF = 0.7
res@vpHeightF = 0.4

res@trXMinF = yrStrt
res@trXMaxF = yrLast
; res@tmXBMode = "Explicit"
; res@tmXBValues = Month
; res@tmXBLabels = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
; res@tmXTOn = False
; res@tiXAxisString = "Month"

res@trYMinF = 22
res@trYMaxF = 36
res@tmYLMode = "Explicit"
res@tmYLValues = ispan(22, 36, 1)
res@tmYLLabels = "P"+res@tmYLValues
; res@tmYROn = False
; res@tiYAxisString = "Temperature (~S~o~N~C)"


res@xyLineThicknesses = Thick
res@xyLineColors = Colors
res@xyDashPatterns = Patterns

plot(0) = gsn_csm_xy(wks, year, x, res)

;==================================================================================
    res_lines                   = True                  ; polyline mods desired
    res_lines@gsLineThicknessF  = 5                     ; line thicker
    

    res_text                    = True                  ; text mods desired
    res_text@txFontHeightF      = 0.018                 ; change text size
    res_text@txJust             = "CenterLeft"          ; text justification

    Texts  = (/"MTG","LBQ","IMD","TTGI","WFI","HOWI"/)

;============================================================
    gx = (/1983,1985/)
    text_pos = 1985.5

    interval = 1
    gy1 = (/25,25/)
    gy2 = gy1 - interval
    gy3 = gy2 - interval
    ; gy4 = gy3 - interval

    gx2 = (/1991,1993/)
    text_pos2 = 1993.5


    res_lines@gsLineColor = Colors(0)
    res_lines@gsLineDashPattern = Patterns(0)
    dum11 = gsn_add_polyline(wks, plot(0), gx, gy1, res_lines)
    dum12 = gsn_add_text(wks, plot(0), Texts(0), text_pos, gy1, res_text)
    res_lines@gsLineColor = Colors(1)
    res_lines@gsLineDashPattern = Patterns(1)
    dum21 = gsn_add_polyline(wks, plot(0), gx, gy2, res_lines)
    dum22 = gsn_add_text(wks, plot(0), Texts(1), text_pos, gy2, res_text)
    res_lines@gsLineColor = Colors(2)
    res_lines@gsLineDashPattern = Patterns(2)
    dum31 = gsn_add_polyline(wks, plot(0), gx, gy3, res_lines)
    dum32 = gsn_add_text(wks, plot(0), Texts(2), text_pos, gy3, res_text)
    res_lines@gsLineColor = Colors(3)
    res_lines@gsLineDashPattern = Patterns(3)
    dum41 = gsn_add_polyline(wks, plot(0), gx2, gy1, res_lines)
    dum42 = gsn_add_text(wks, plot(0), Texts(3), text_pos2, gy1, res_text)
    res_lines@gsLineColor = Colors(4)
    res_lines@gsLineDashPattern = Patterns(4)
    dum51 = gsn_add_polyline(wks, plot(0), gx2, gy2, res_lines)
    dum52 = gsn_add_text(wks, plot(0), Texts(4), text_pos2, gy2, res_text)
    res_lines@gsLineColor = Colors(5)
    res_lines@gsLineDashPattern = Patterns(5)
    dum61 = gsn_add_polyline(wks, plot(0), gx2, gy3, res_lines)
    dum62 = gsn_add_text(wks, plot(0), Texts(5), text_pos2, gy3, res_text)

    draw(plot(0))
end