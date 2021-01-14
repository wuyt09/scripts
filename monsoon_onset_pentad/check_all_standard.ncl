
begin

bob_lat1 = 5
bob_lat2 = 15
bob_lon1 = 90
bob_lon2 = 100

scs_lat1 = 10
scs_lat2 = 20
scs_lon1 = 110
scs_lon2 = 120

ism_lat1 = 10
ism_lat2 = 20
ism_lon1 = 60
ism_lon2 = 85


yrStrt = 1979
yrLast = 2015
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

diri = "/home/yangsong3/wuyt/sysu/scripts/monsoon_onset_pentad/"
fili1 = "wind_direction_10m_onset_BOB_SCS_ISM.txt"
fili2 = "MTG_onset_BOB_SCS_ISM.txt"
fili3 = "prec_ERA5_onset_BOB_SCS_ISM.txt"

x1 = asciiread(diri+fili1, (/nyear,4/), "integer")
x2 = asciiread(diri+fili2, (/nyear,4/), "integer")
x3 = asciiread(diri+fili3, (/nyear,4/), "integer")


fili4 = "LBQ_onset_pentad_1979_2017_wind.txt"
fili5 = "LBQ_onset_pentad_1979_2017_MTG.txt"
fili6 = "LBQ_onset_pentad_1979_2017_prec.txt"

x4 = asciiread(diri+fili4, (/nyear,4/), "integer")
x5 = asciiread(diri+fili5, (/nyear,4/), "integer")
x6 = asciiread(diri+fili6, (/nyear,4/), "integer")


monsoons = (/"BOB","SCS","ISM"/)
x = new((/6,nyear,3/), integer)
x(0,:,:) = x1(:,1:3)
x(1,:,:) = x2(:,1:3)
x(2,:,:) = x3(:,1:3)
x(3,:,:) = x4(:,1:3)
x(4,:,:) = x5(:,1:3)
x(5,:,:) = x6(:,1:3)

ll = 2
yy = x(:,:,ll)
;=====================================================================
; plot
;=====================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/monsoon_onset_pentad/"
figname = "check_all_standard"

wks = gsn_open_wks("pdf", figpath+figname)
plot = new(1, graphic)
Colors = (/"black","red","blue","black","red","blue"/)
Patterns = (/2,2,2,0,0,0/)
Thick = (/10,10,10,10,10,10/)

res = True
res@gsnFrame = False
res@gsnDraw  = False
res@gsnMaximize = True

res@vpWidthF = 0.7
res@vpHeightF = 0.45

res@xyLineThicknesses = Thick
res@xyLineColors = Colors
res@xyDashPatterns = Patterns

res@trXMinF = yrStrt
res@trXMaxF = yrLast
res@trYMinF = 18
res@trYMaxF = 36
res@tmYLPrecision = 2
res@tiYAxisString = "Pentad"

res@gsnYRefLine = 0
res@gsnYRefLineColor = "grey40"
res@gsnYRefLineThicknessF = 5
res@gsnYRefLineDashPattern = 2

res@gsnRightString = ""
res@gsnStringFontHeightF = 0.025

res@gsnLeftString = monsoons(ll)
plot(0) = gsn_csm_xy(wks, year, yy, res)

;==================================================================================
    res_lines                   = True                  ; polyline mods desired
    res_lines@gsLineThicknessF  = 5                     ; line thicker
    

    res_text                    = True                  ; text mods desired
    res_text@txFontHeightF      = 0.018                 ; change text size
    res_text@txJust             = "CenterLeft"          ; text justification

    Texts  = (/"wind","MTG","prec","wind-L","MTG-L","prec-L"/)

;============================================================
    gx = (/1983,1985/)
    text_pos = 1985.5

    interval = 1
    gy1 = (/35,35/)
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

draw(plot)
end