
begin

yrStrt = 1979
yrLast = 2015
year = ispan(yrStrt, yrLast, 1)
Month = ispan(1, 12, 1)

diri = "/home/yangsong3/wuyt/data/ERA-Interim/monthly/"
fili = "q1.mon.mean.1979.2016.ERAI.nc" ; ERAI 

f = addfile(diri+fili, "r")
TIME = f->time
YYYY = cd_calendar(TIME, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

T = short2flt(f->q1(iYYYY,:,::-1,:))
T = T*1004.07 ;*86400
printVarSummary(T)

filisp = "/home/yangsong3/data-observation/ERA-interim/monthly/Surface_pressure.mon.mean.nc"
fsp = addfile(filisp, "r")
sp0 = short2flt(fsp->sp(iYYYY,::-1,:))
sp = linint2_Wrap(sp0&longitude, sp0&latitude, sp0, True, T&lon, T&lat, 0)
printVarSummary(sp)

filitp = "~/wuyt/data/topo.sfc.1x1.nc"
ftp = addfile(filitp, "r")
topo0 = ftp->topo(::-1,:)
topo = linint2_Wrap(topo0&lon, topo0&lat, topo0, True, T&lon, T&lat, 0)
printVarSummary(topo)

;===========integration========================================
lev = T&lev
plev = lev*100.
psfc = sp
ptop = min(plev)

dP = dpres_plevel_Wrap(plev, psfc, ptop, 0)
dP!1 = "lev"
dP&lev = lev
printVarSummary(dP)

; dp = conform_dims(dimsizes(T), dP, (//))
Tdp = T*dP
copy_VarCoords(T, Tdp)

; levsfc = 400
; levtop = 200

;TTT = dim_sum_n_Wrap(Tdp(:,{levsfc:levtop},:,:), 1)/dim_sum_n_Wrap(dP(:,{levsfc:levtop},:,:),1)
TTT = dim_sum_n_Wrap(Tdp, 1) ;/dim_sum_n_Wrap(dP,1)
copy_VarCoords(T(:,0,:,:), TTT)
TTT = TTT/9.8
printVarSummary(TTT)


tttx = clmMonTLL(TTT)
printVarSummary(tttx)
ttt = mask(tttx, topo.ge.1500, True)
copy_VarCoords(tttx, ttt)
print(ttt(4,{30},{115}))

ttt_12mn = wgt_areaave_Wrap(ttt(:,{25:40},{75:105}), 1, 1, 0)

dt = ttt_12mn
do i = 1, 11
    dt(i) = ttt_12mn(i) - ttt_12mn(i-1)
end do
dt(0) = ttt_12mn(0) - ttt_12mn(11)
;dt = abs(dt)
print(dt)

rc1 = regline(Month(2:5), ttt_12mn(2:5))
rc2 = regline(Month(7:10), ttt_12mn(7:10))
print("rc1 "+rc1+"  rc2 "+rc2) ;rc1 1.48554  rc2 -1.20139

;================================================================
; plot
;================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_Spr_Aut_heating/obs/subseasonal/"
figname = "draw_monthly_Q1_ERAI"

lonstr = 60 ;20
lonend = 110 ;220
latstr = 20
latend = 50

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
res@cnLevelSpacingF = 10.
res@cnMaxLevelValF = 250.
res@cnMinLevelValF = -250.
res@cnFillPalette = "CBR_coldhot"

do i = 0, 11
    res@gsnLeftString = "Month = "+(i+1)
    plot(i) = gsn_csm_contour_map(wks, tttx(i,:,:), res)
end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
;pres@gsnPanelMainString = "tropo ("+levsfc+"hPa-"+levtop+"hPa) temp"

gsn_panel(wks, plot, (/4,3/), pres)

;====================================================================
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
res@tiYAxisString = "(W/m~S~2~N~)"

;Colors = (/"deeppink2", "blue", "purple","green4","orange2", "black" /);
Colors = (/"black","blue","forestgreen"/)
Patterns = (/0,0,0/)
res@xyLineThicknessF = 5
res@xyLineColors = Colors
res@xyDashPatterns = Patterns
res@gsnLeftString = "Q1"
res@gsnStringFontHeightF = 0.025


plot(0) = gsn_csm_xy(wks, Month, ttt_12mn, res)
draw(plot(0))
frame(wks)



end