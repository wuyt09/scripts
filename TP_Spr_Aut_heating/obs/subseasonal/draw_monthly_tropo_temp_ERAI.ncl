
begin

yrStrt = 1979
yrLast = 2015
year = ispan(yrStrt, yrLast, 1)
Month = ispan(1, 12, 1)
month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)

diri = "/home/yangsong3/data-observation/ERA-interim/monthly/"
fili = "air.mon.mean.nc"

f = addfile(diri+fili, "r")
TIME = f->time
YYYY = cd_calendar(TIME, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

T = short2flt(f->t(iYYYY,:,::-1,:))
printVarSummary(T)

filisp = "Surface_pressure.mon.mean.nc"
fsp = addfile(diri+filisp, "r")
sp = short2flt(fsp->sp(iYYYY,:,:))
printVarSummary(sp)

filitp = "~/wuyt/data/topo.sfc.1x1.nc"
ftp = addfile(filitp, "r")
topo = ftp->topo(::-1,:)
printVarSummary(topo)

;===========integration========================================
lev = T&level
plev = lev*100.
psfc = sp
ptop = min(plev)

dP = dpres_plevel_Wrap(plev, psfc, ptop, 0)
dP!1 = "lev"
dP&lev = lev
printVarSummary(dP)

; dp = conform_dims(dimsizes(T), dP, 1)
Tdp = T*dP
copy_VarCoords(T, Tdp)

levsfc = 500
levtop = 200

TTT = dim_sum_n_Wrap(Tdp(:,{levsfc:levtop},:,:), 1)/dim_sum_n_Wrap(dP(:,{levsfc:levtop},:,:),1)
copy_VarCoords(T(:,0,:,:), TTT)
printVarSummary(TTT)


ttty = clmMonTLL(TTT)
ttt_zonal = dim_avg_n_Wrap(ttty, 1)
ttt_zonaly = conform_dims(dimsizes(ttty), ttt_zonal, (/0,2/))
tttx = ttty - ttt_zonaly
copy_VarCoords(ttty, tttx)
printVarSummary(tttx)
ttt = mask(tttx, topo.ge.1000, True)
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
figname = "draw_monthly_tropo_temp_ERAI"

lonstr = 60 ;20
lonend = 150 ;220
latstr = 10
latend = 60

wks = gsn_open_wks("pdf", figpath+figname)
plot = new(8, graphic)
plot_topo = plot

 res_topo            = True
 res_topo@cnFillOn   = False
 res_topo@cnLinesOn  = True
 res_topo@gsnFrame   = False
 res_topo@gsnDraw    = False
 res_topo@cnLineLabelsOn         = False
 res_topo@cnLineThicknessF       = 3
 res_topo@cnLineColor            = "grey45"
 res_topo@cnInfoLabelOn          = False
 res_topo@cnLevelSelectionMode   ="ExplicitLevels"
 res_topo@cnLevels               = (/1500/)
 res_topo@gsnLeftString = ""
 res_topo@gsnRightString = ""

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
res@cnLevelSpacingF = 2.
res@cnMaxLevelValF = 30.
res@cnMinLevelValF = -30.
res@cnFillPalette = "CBR_coldhot"
; res@cnLevelSpacingF = 1.
; res@cnMaxLevelValF = 255.
; res@cnMinLevelValF = 240.
; res@cnFillPalette = "sunshine_9lev"

do i = 0, 3
    res@gsnLeftString = month(i+2)
    plot(i) = gsn_csm_contour_map(wks, tttx(i+2,:,:), res)
    plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(i), plot_topo(i))
    res@gsnLeftString = month(i+7)
    plot(i+4) = gsn_csm_contour_map(wks, tttx(i+7,:,:), res)
    plot_topo(i+4) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(i+4), plot_topo(i+4))
end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = "Tropo ("+levsfc+"hPa-"+levtop+"hPa) temp"
pres@pmLabelBarOrthogonalPosF = -0.02

gsn_panel(wks, plot, (/2,4/), pres)

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
res@tiYAxisString = "Temperature (~S~o~N~C)"

;Colors = (/"deeppink2", "blue", "purple","green4","orange2", "black" /);
Colors = (/"black","blue","forestgreen"/)
Patterns = (/0,0,0/)
res@xyLineThicknessF = 5
res@xyLineColors = Colors
res@xyDashPatterns = Patterns
res@gsnLeftString = "ERAI tropo temp"
res@gsnStringFontHeightF = 0.025

plot(0) = gsn_csm_xy(wks, Month, ttt_12mn, res)
draw(plot(0))
frame(wks)



end