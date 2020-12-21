
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
yrLast = 2017
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

x = new((/3,nyear,3/), integer)
; x(0,:,:) = x1(:,1:3)
x(1,:,:) = x2(:,1:3)
x(2,:,:) = x3(:,1:3)

x_max = dim_max_n_Wrap(x, 0)
; print("BOB "+x_max(:,0))
; print("SCS "+x_max(:,1))
; print("ISM "+x_max(:,2))
x_max(20,0) = 19
x_max(21,0) = 20

write_table("multi_standard_onset_BOB_SCS_ISM.txt", "w", [/"YEAR","BOB","SCS","ISM"/], "%4s %4s %4s %4s")
write_table("multi_standard_onset_BOB_SCS_ISM.txt", "a", [/year,x_max(:,0),x_max(:,1),x_max(:,2)/], "%0.4i %0.2i %0.2i %0.2i")


yy = new((/6,nyear/), integer)
yy(0,:) = x_max(:,0)
yy(1,:) = x_max(:,1)
yy(2,:) = x_max(:,2)
yy(3,:) = (/27,26,26,23,27,22,22,26,27,25,24,25,24,26,27,25,26,23,28,27,20,22,24,26,27,26,26,23,25,25,21,28,24,23,27,26,27,28,27/)
yy(4,:) = (/29,28,30,30,29,28,30,27,32,29,32,29,32,28,29,30,28,34,28,29,28,27,28,28,28,34,29,30,28,29,29,29,27,27,27,31,33,32,29/)
yy(5,:) = (/33,31,30,30,33,30,30,31,31,29,31,27,32,33,30,30,33,31,33,32,28,30,32,32,32,27,32,29,31,30,31,29,30,32,30,32,32,31,31/)

x_max!0 = "year"
x_max!1 = "component"
;=====================================================================
; plot
;=====================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/monsoon_onset_pentad/"
figname = "3_onset_pentad_for_ASM"

wks = gsn_open_wks("pdf", figpath+figname)
plot = new(1, graphic)
Colors = (/"black","red","blue","black","red","blue"/)
Patterns = (/0,0,0,2,2,2/)
Thick = (/10,10,10,5,5,5/)

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

plot(0) = gsn_csm_xy(wks, year, yy, res)
draw(plot)

end