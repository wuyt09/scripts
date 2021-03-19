
begin

lonstr = 0
lonend = 360 ;130
latstr = 0 ;10
latend = 70

yrStrt = 1979
yrLast = 2019
year = ispan(yrStrt, yrLast, 1)
nyear = dimsizes(year)
season = "JJA"

LEV = 500

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

hgt_zonal = dim_avg_n_Wrap(hgt, 1)
printVarSummary(hgt_zonal)

x = conform_dims(dimsizes(hgt), hgt_zonal, (/0,2/))
hgt_new = hgt - x
copy_VarCoords(hgt, hgt_new)
printVarSummary(hgt_new)

hgt_clm = dim_avg_n_Wrap(hgt_new, 0)
printVarSummary(hgt_clm)

hgt_F10 = dim_avg_n_Wrap(hgt_new(0:9,:,:), 0)
hgt_L10 = dim_avg_n_Wrap(hgt_new(30:39,:,:), 0)
printVarSummary(hgt_L10)
; print(hgt_new(:,{30},{120}))
; exit
;===============================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_surrounding/"
figname = season+"_H500_Clim_F10_L10"
wks = gsn_open_wks("pdf", figpath+figname)

plot = new(3, graphic)
plot_plot = new(3, graphic)

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

;-----------------------------------------------------
res = True
res@gsnDraw = False
res@gsnFrame = False
res@gsnMaximize = True

res@mpOutlineOn = True
res@mpFillOn = True
res@mpGeophysicalLineThicknessF = 5
res@mpGeophysicalLineColor = "grey25" ;"red" 
; res@mpFillDrawOrder = "PostDraw"
; res@mpCenterLonF = 180
res@mpMinLatF = latstr
res@mpMaxLatF = latend
res@mpMinLonF = lonstr
res@mpMaxLonF = lonend
res@mpCenterLonF = (lonstr+lonend)/2


res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False
; res@cnFillMode = "CellFill"

res@gsnStringFontHeightF = 0.025
res@tmXBLabelFontHeightF = 0.02
res@tmYLLabelFontHeightF = 0.02

res@cnFillPalette = "BlueWhiteOrangeRed"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -300
res@cnMaxLevelValF = 300
res@cnLevelSpacingF = 50

; res@lbLabelBarOn = False

res@pmLabelBarHeightF = 0.1
res@pmLabelBarWidthF = 0.7
res@pmLabelBarOrthogonalPosF = 0.1


res@gsnLeftString = "Clim"
plot(0) = gsn_csm_contour_map(wks, hgt_clm, res)

res@gsnLeftString = "F10"
plot(1) = gsn_csm_contour_map(wks, hgt_F10, res)

res@gsnLeftString = "L10"
plot(2) = gsn_csm_contour_map(wks, hgt_L10, res)


pres = True
pres@gsnMaximize = True
gsn_panel(wks, plot, (/3,1/), pres)
end