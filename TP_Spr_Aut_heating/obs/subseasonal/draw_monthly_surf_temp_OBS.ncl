
begin

; yrStrt = 1979
; yrLast = 2015
; year = ispan(yrStrt, yrLast, 1)
Month = ispan(1, 12, 1)
diri = "/home/yangsong3/wuyt/data/Station/TEMP/scripts/"
fili = "3-write_txt2nc.nc"

yrbegin = 1951
yrStrt = 1951
yrLast = 2012
strt_index = yrStrt - yrbegin
last_index = yrLast - yrbegin

f = addfile(diri+fili, "r")
TEMP = f->TAS ;(:,strt_index:last_index,:)
LAT  = f->LAT ;(:,strt_index:last_index,:)
LON  = f->LON ;(:,strt_index:last_index,:)
HGT  = f->HGT ;(:,strt_index:last_index,:)
staID = f->Station
year  = f->Year ;(strt_index:last_index)
month = f->Month
printVarSummary(TEMP)
printVarSummary(LON)

nstations = dimsizes(staID)
nyear = dimsizes(year)
nmonth = dimsizes(month)

;==========================================================
; task1: remove some stations for too much fillvalue
;==========================================================
; Num = new(297, integer)

; do i = 0, nstations-1
;     Num(i) = num(.not.ismissing(TEMP(i,:,:)))
; end do
; ;print(Num)

; index = ind(Num .ge. (nyear-5)*nmonth)
; ;print(staID(index))

hgt = dim_avg_n_Wrap(HGT, (/1,2/))
lat = dim_avg_n_Wrap(LAT, (/1,2/))
lon = dim_avg_n_Wrap(LON, (/1,2/))

index = ind(hgt.ge.1500)
lat0 = LAT(index,:,:) ; cause lat/lon do not change with time
lon0 = LON(index,:,:)
latx = dim_avg_n_Wrap(lat0, (/1,2/))
lonx = dim_avg_n_Wrap(lon0, (/1,2/))
printVarSummary(latx)
print(lonx)

temp = TEMP(index,:,:)
printVarSummary(temp)

ttt_12mn = dim_avg_n_Wrap(temp, (/0,1/))
printVarSummary(ttt_12mn)

;=========================================================
file_path = "./"
file_name = "All_stations_1500m_12mon_clm.txt"
system("rm -rf "+file_path+file_name)

; write_table(file_path+file_name, "w", (["Temp"]), "%10s")
; write_table(file_path+file_name, "a", ttt_12mn, "%10.2f")
asciiwrite(file_path+file_name, ttt_12mn)
;===========================================================

dt = ttt_12mn
do i = 1, 11
    dt(i) = ttt_12mn(i) - ttt_12mn(i-1)
end do
dt(0) = ttt_12mn(0) - ttt_12mn(11)
;dt = abs(dt)
print(dt)
; (0)     -1.048993
; (1)     2.787453
; (2)     4.528405
; (3)     4.620275
; (4)     3.792753
; (5)     2.655948
; (6)     1.42623
; (7)     -0.6836548
; (8)     -3.064156
; (9)     -4.846241
; (10)    -5.854042
; (11)    -4.313979

rc1 = regline(Month(2:4), ttt_12mn(2:4))
rc2 = regline(Month(8:10), ttt_12mn(8:10))
print("rc1 "+rc1+"  rc2 "+rc2) ;rc1 3.69997  rc2 -4.61396

;================================================================
; plot
;================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_Spr_Aut_heating/obs/"
figname = "draw_monthly_surf_temp_OBS"

lonstr = 60 ;20
lonend = 110 ;220
latstr = 20
latend = 50

wks = gsn_open_wks("pdf", figpath+figname)
plot = new(12, graphic)

; res = True
; res@gsnFrame = False
; res@gsnDraw  = False

; res@mpOutlineOn = True
; res@mpMinLonF = lonstr
; res@mpMaxLonF = lonend
; res@mpMinLatF = latstr
; res@mpMaxLatF = latend

; res@gsnStringFontHeightF = 0.035
; res@tmXBLabelFontHeightF = 0.025
; res@tmYLLabelFontHeightF = 0.025

; res@cnFillOn = True
; res@cnLinesOn = False
; res@cnLineLabelsOn = False
; res@cnInfoLabelOn = False

; res@lbLabelBarOn = False
; res@cnLevelSelectionMode = "ManualLevels"
; res@cnLevelSpacingF = 2.
; res@cnMaxLevelValF = 300.
; res@cnMinLevelValF = 250.
; res@cnFillPalette = "CBR_coldhot"

; do i = 0, 11
;     res@gsnLeftString = "Month = "+(i+1)
;     plot(i) = gsn_csm_contour_map(wks, ttt(i,:,:), res)
; end do

; pres = True
; pres@gsnMaximize = True
; pres@gsnPanelLabelBar = True
; pres@gsnPanelMainString = "surf temp"

; gsn_panel(wks, plot, (/4,3/), pres)

;====================================================================
;delete(res)

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
res@tmXBLabels = Month
res@tmXTOn = False
res@tiXAxisString = "Month"

; res@trYMinF = -16
; res@trYMaxF = 20
; res@tmYLMode = "Explicit"
; res@tmYLValues = ispan(-15, 20, 5)
; res@tmYLLabels = ispan(-15, 20, 5)
res@tmYROn = False
res@tiYAxisString = "Temperature (~S~o~N~C)"

;Colors = (/"deeppink2", "blue", "purple","green4","orange2", "black" /);
Colors = (/"black","blue","forestgreen"/)
Patterns = (/0,0,0/)
res@xyLineThicknessF = 5
res@xyLineColors = Colors
res@xyDashPatterns = Patterns
res@gsnLeftString = "Station data"


plot(0) = gsn_csm_xy(wks, Month, ttt_12mn, res)
draw(plot(0))
frame(wks)



end