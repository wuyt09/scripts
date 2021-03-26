
begin

lonstr = 0
lonend = 130
latstr = 15
latend = 90


yrStrt = 1979
yrLast = 2019
year = ispan(yrStrt, yrLast, 1)
nyear = dimsizes(year)
season = "JJA"

LEV = (/200,500,850/)

latstr_east = 35
latend_east = 55
lonstr_east = 84 ;80
lonend_east = 115 ;120

latstr_west = 40
latend_west = 55
lonstr_west = 10
lonend_west = 50

latstr_south = 15
latend_south = 25
lonstr_south = 60
lonend_south = 110

latstr_north = 45
latend_north = 60
lonstr_north = 70
lonend_north = 110

latstr_tp = 25 ;28
latend_tp = 40 ;38
lonstr_tp = 75 ;80
lonend_tp = 105

region = (/"West","East","South","North","TP"/)
;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)
printVarSummary(topo)

;===================================================================

diri = "/home/yangsong3/data-observation/ERA5-monthly/pressure/"
fili = "hgt_mon_1x1_1979_2019.nc"

f = addfile(diri+fili, "r")

Time = f->time
YYYY = cd_calendar(Time, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

Hgt = short2flt(f->z(iYYYY,{LEV},:,:))
Hgt = Hgt/9.8
Hgt@units = "gpm"

hgt = month_to_season(Hgt, season)
printVarSummary(hgt)

hgt_zonal = dim_avg_n_Wrap(hgt, 3)
printVarSummary(hgt_zonal)

x = conform_dims(dimsizes(hgt), hgt_zonal, (/0,1,2/))
hgt_new = hgt - x
copy_VarCoords(hgt, hgt_new)
hgt_new&time = year
printVarSummary(hgt_new)

rc = regCoef_n(year, hgt_new, 0, 0)
copy_VarCoords(hgt_new(0,:,:,:), rc)
; hgt_new = hgt

tval = onedtond(rc@tval, dimsizes(rc))
df   = onedtond(rc@nptxy, dimsizes(rc))
b    = tval
b    = 0.5
prob = betainc(df/(df+tval^2), df/2.0, b)
copy_VarCoords(rc, prob)

rc@long_name = "regression coefficient"
prob@long_name = "probability"

trend = rc*10
copy_VarCoords(rc, trend)
; trend@units = "K/decade"
printVarSummary(trend)

sig = 0.05
;===============================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_surrounding/"
figname = "trend_"+season+"_H500_ERA5"
wks = gsn_open_wks("pdf", figpath+figname)

plot = new(3, graphic)
plot_topo = plot
plot_prob = plot
dum1 = plot
dum2 = plot
dum3 = plot
dum4 = plot

 res_topo            = True
 res_topo@cnFillOn   = False
 res_topo@cnLinesOn  = True
 res_topo@gsnFrame   = False
 res_topo@gsnDraw    = False
 res_topo@cnLineLabelsOn         = False
 res_topo@cnLineThicknessF       = 5
 res_topo@cnLineColor            = "forestgreen"
 res_topo@cnInfoLabelOn          = False
 res_topo@cnLevelSelectionMode   ="ExplicitLevels"
 res_topo@cnLevels               = (/2000/)
 res_topo@gsnLeftString = ""
 res_topo@gsnRightString = ""

;----------------------------------------------
; draw prob
;----------------------------------------------
 res3  = True
 res3@gsnDraw   = False
 res3@gsnFrame  = False

 ;res3@cnFillOn  = False 
 res3@cnLinesOn = False
 res3@cnLineLabelsOn = False
 res3@cnInfoLabelOn  = False
 ;res3@lbLabelBarOn   = False

 res3@gsnRightString = ""

 res3@cnLevelSelectionMode = "ManualLevels"
 res3@cnMinLevelValF = 0.00
 res3@cnMaxLevelValF = 0.10
 res3@cnLevelSpacingF = 0.01

 ; res3@cnFillPattern = 17 ;(/17/)
 res3@cnFillColor  =  "black"
 res3@gsnLeftString = ""
 res3@cnFillDotSizeF = 0.003
 res3@cnFillScaleF  = 1.2
 ;----------------------------------------------
;-----------------------------------------------------
res = True
res@gsnDraw = False
res@gsnFrame = False
res@gsnMaximize = True

res@mpOutlineOn = True
res@mpFillOn = True
res@mpGeophysicalLineThicknessF = 3
res@mpGeophysicalLineColor = "grey25" ;"red" 
; res@mpFillDrawOrder = "PostDraw"
res@mpCenterLonF = 90
res@mpMinLatF = latstr
res@mpMaxLatF = latend
; res@mpMinLonF = lonstr
; res@mpMaxLonF = lonend
; res@mpCenterLonF = (lonstr+lonend)/2
res@gsnPolar = "NH"

res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False
; res@cnFillMode = "CellFill"

res@gsnStringFontHeightF = 0.03
; res@tmXBLabelFontHeightF = 0.03
; res@tmYLLabelFontHeightF = 0.03
res@gsnPolarLabelFontHeightF= 0.02

res@cnFillPalette = "BlueWhiteOrangeRed"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -10
res@cnMaxLevelValF = 10
res@cnLevelSpacingF = 2

res@lbLabelBarOn = False

res@pmLabelBarHeightF = 0.1
res@pmLabelBarWidthF = 0.7
res@pmLabelBarOrthogonalPosF = 0.1


; res@gsnLeftString = "Clim"
; plot(0) = gsn_csm_contour_map(wks, hgt_clm, res)


;-----------------------------------
gres = True
gres@gsLineColor = "purple"
gres@gsLineThicknessF = 4
gres@gsLineDashPattern = 16

 glon1 = (/lonend_north,lonend_north,lonstr_north,lonstr_north,lonend_north/)
 glat1 = (/latstr_north,latend_north,latend_north,latstr_north,latstr_north/) ; north
 glon2 = (/lonend_east,lonend_east,lonstr_east,lonstr_east,lonend_east/)
 glat2 = (/latstr_east,latend_east,latend_east,latstr_east,latstr_east/) ; east
 glon3 = (/lonend_south,lonend_south,lonstr_south,lonstr_south,lonend_south/)
 glat3 = (/latstr_south,latend_south,latend_south,latstr_south,latstr_south/) ; south
 glon4 = (/lonend_west,lonend_west,lonstr_west,lonstr_west,lonend_west/)
 glat4 = (/latstr_west,latend_west,latend_west,latstr_west,latstr_west/) ; west
;-----------------------------------
do i = 0, 2
    res@gsnLeftString = LEV(i)+"hPa "+season+" trend"
    plot(i) = gsn_csm_contour_map_polar(wks, trend(i,:,:), res)
    plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(i), plot_topo(i))

    plot_prob(i) = gsn_csm_contour(wks,prob(i,:,:),res3)
    plot_prob(i) = ShadeLtContour(plot_prob(i),sig,17)
    overlay(plot(i), plot_prob(i))

    ; dum1 = gsn_add_polyline(wks, plot(0), glon1, glat1, gres)
    dum2(i) = gsn_add_polyline(wks, plot(i), glon2, glat2, gres)
    dum3(i) = gsn_add_polyline(wks, plot(i), glon3, glat3, gres)
    dum4(i) = gsn_add_polyline(wks, plot(i), glon4, glat4, gres)
end do


pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
pres@lbLabelFontHeightF = 0.012
; pres@gsnPanelMainString = "H500"

; gsn_panel(wks, plot, (/nx/3+1,3/), pres)
gsn_panel(wks, plot, (/1,3/), pres)
end