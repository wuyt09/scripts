
begin

yrStrt = 1979
yrLast = 2015
year = ispan(yrStrt, yrLast, 1)
Month = ispan(1, 12, 1)
LEV = 200

;=================================================================
filitp = "~/wuyt/data/topo.sfc.1x1.nc"
ftp = addfile(filitp, "r")
topo = ftp->topo(::-1,:)
printVarSummary(topo)
;=================================================================

diri = "/home/yangsong3/data-observation/ERA-interim/monthly/"
fili = "air.mon.mean.nc"

f = addfile(diri+fili, "r")
TIME = f->time
YYYY = cd_calendar(TIME, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

T = short2flt(f->t(iYYYY,{LEV},::-1,:))
;T = T - 273.15
printVarSummary(T)

t = clmMonTLL(T)
;print(t(11,{25:40},{75:135}))
;============================================================
; lon = T&longitude
; lat = T&latitude
; ;print(lat)

; dlon = (lon(2)-lon(1))*0.0174533
; dtdx = t
; dtdx = 0
; do i = 0, dimsizes(lat)-1
;     dx = 6378388.*cos(0.0174533*lat(i))*dlon  
;     dtdx(:,i,:) = center_finite_diff_n(t(:,i,:), dx, True, 0, 1)
; end do
; printVarSummary(dtdx)
; dt = wgt_areaave_Wrap(dtdx(:,{25:40},{75:135}), 1, 1, 0)
; ;dt = dim_sum_n_Wrap(dtdx(:,{25:40},{75:135}), (/1,2/))
; printVarSummary(dt)
; ;print(dtdx(11,{25:40},{75:135}))

; dk = dt
; do i = 1, 11
;     dk(i) = dt(i) - dt(i-1)
; end do
; dk(0) = dt(0) - dt(11)

; rc1 = regline(Month(2:5), dt(2:5))
; rc2 = regline(Month(7:10), dt(7:10))
; print("rc1 "+rc1+" rc2 "+rc2)
;============================================================
; t2 = wgt_areaave_Wrap(t(:,{25:40},{105:135}), 1, 1, 0)
; t3 = wgt_areaave_Wrap(t(:,{25:40},{120:150}), 1, 1, 0)
; t = mask(t, topo.ge.1500, True)
; t1 = wgt_areaave_Wrap(t(:,{25:40},{75:105}), 1, 1, 0)
; deltat = t1
deltat = wgt_areaave_Wrap(t(:,{25:40},{75:105}), 1, 1, 0)
printVarSummary(deltat)

dk = deltat
do i = 1, 11
    dk(i) = deltat(i) - deltat(i-1)
end do
dk(0) = deltat(0) - deltat(11)

rc1 = regline(Month(2:5), deltat(2:5))
rc2 = regline(Month(7:10), deltat(7:10))
print("rc1 "+rc1+" rc2 "+rc2)
;============================================================
; pp = new((/2,12/), float)
; pp(0,:) = (/dt/)
; pp(1,:) = (/deltat/) 

;================================================================
; plot
;================================================================
lonstr = 60 ;20
lonend = 110 ;220
latstr = 20
latend = 50

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_Spr_Aut_heating/obs/subseasonal/"
figname = "draw_monthly_T200"

wks = gsn_open_wks("pdf", figpath+figname)
plot = new(12, graphic)

res = True
res@gsnFrame = False
res@gsnDraw  = False

res@mpOutlineOn = True
res@mpMinLonF = lonstr
res@mpMaxLonF = lonend
res@mpMinLatF = latstr
res@mpMaxLatF = latend

res@gsnStringFontHeightF = 0.035
res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025

res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False

res@lbLabelBarOn = False
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 1.
res@cnMaxLevelValF = 230.
res@cnMinLevelValF = 210.
res@cnFillPalette = "CBR_coldhot"

do i = 0, 11
    res@gsnLeftString = "Month = "+(i+1)
    plot(i) = gsn_csm_contour_map(wks, t(i,:,:), res)
end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = "T 200hPa"

gsn_panel(wks, plot, (/4,3/), pres)

delete(res)
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
;res@tiYAxisString = "Temperature (~S~o~N~C)"

res@gsnYRefLine = 0
res@gsnYRefLineThicknessF = 5
res@gsnYRefLineColor = "grey53"
res@gsnYRefLineDashPattern = 2

;Colors = (/"deeppink2", "blue", "purple","green4","orange2", "black" /);
Colors = (/"black","blue","forestgreen"/)
Patterns = (/0,0,0/)
res@xyLineThicknessF = 5
res@xyLineColors = Colors
res@xyDashPatterns = Patterns
res@gsnLeftString = "Temp 500hPa: east - west, "+sprintf("%4.2f", rc1)+" / "+sprintf("%4.2f", rc2)
res@gsnStringFontHeightF = 0.025


plot(0) = gsn_csm_xy(wks, Month, deltat, res)
draw(plot(0))
frame(wks)


end