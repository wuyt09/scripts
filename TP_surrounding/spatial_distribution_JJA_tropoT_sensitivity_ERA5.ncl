
begin

lonstr = 0
lonend = 360
latstr = 10 ;10
latend = 60

yrStrt = 1979
yrLast = 2018
year = ispan(yrStrt, yrLast, 1)
nyear = dimsizes(year)
season = "JJA"

; LEVs = 500
; LEVe = 200
; LEVs = 850
; LEVe = 200
; LEVs = 850
; LEVe = 300
; LEVs = 1000
; LEVe = 300
LEVs = (/500,850,850,1000/)
LEVe = (/200,200,300,300/)

latstr_east = 25
latend_east = 45
lonstr_east = 105
lonend_east = 120

latstr_west = 30
latend_west = 55
lonstr_west = 30
lonend_west = 70

latstr_south = 10
latend_south = 25
lonstr_south = 75
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
fili = "temp.monthly.197901-201812.nc"
f = addfile(diri+fili, "r")

Time = f->time
YYYY = cd_calendar(Time, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

t = short2flt(f->t(iYYYY,::-1,::-1,:))
printVarSummary(t)

lat = t&latitude
lon = t&longitude
nlat = dimsizes(lat)
nlon = dimsizes(lon)

filips = "/home/yangsong3/data-observation/ERA5-monthly/surface/sur_pressure_mon_1x1_1979_2019.nc"
fps = addfile(filips, "r")
ps = short2flt(fps->sp(iYYYY,::-1,:))
printVarSummary(ps)


plev = tofloat(t&level)
plev = plev*100
plev!0 = "lev"
plev&lev = t&level
plev@units = "Pa"
; print(plev)

dp = dpres_plevel_Wrap(plev, ps, min(plev), 0) ;LEVs(i)*100
; dP = dpres_plevel_Wrap(plev, max(plev), min(plev), 0)
; dp = conform_dims(dimsizes(t), dP, 1)
; copy_VarCoords(t, dp)
printVarSummary(dp)
; print(dp(90,:,{35},{90}))
; exit

tdp = t*dp
copy_VarCoords(t, tdp)
printVarSummary(tdp)

t_tropo = new((/dimsizes(iYYYY),4,nlat,nlon/), float)

do i = 0, 3
    t_tropo(:,i,:,:) = dim_sum_n_Wrap(tdp(:,{LEVs(i):LEVe(i)},:,:), 1)/dim_sum_n_Wrap(dp(:,{LEVs(i):LEVe(i)},:,:), 1)
    copy_VarCoords(t(:,0,:,:), t_tropo(:,i,:,:))
end do


t2m_clm = month_to_season(t_tropo, season)
printVarSummary(t2m_clm)

rc = regCoef_n(year, t2m_clm, 0, 0)
copy_VarCoords(t2m_clm(0,:,:,:), rc)

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
trend@units = "K/decade"
printVarSummary(trend)

sig = 0.05
;===============================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_surrounding/"
figname = "spatial_distribution_"+season+"_tropoT_trend_sensitivity_ERA5"
wks = gsn_open_wks("pdf", figpath+figname)

plot = new(4, graphic)
plot_topo = plot
plot_prob = new(4, graphic)
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
 ;----------------------------------------------

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

res@cnFillPalette = "CBR_coldhot"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -0.5
res@cnMaxLevelValF = 0.5
res@cnLevelSpacingF = 0.1
res@lbLabelBarOn = False


;-----------------------------------
gres = True
gres@gsLineColor = "blue"
gres@gsLineThicknessF = 5.0
gres@gsLineDashPattern = 16
 glon1 = (/110,110,70,70,100/)
 glat1 = (/45,60,60,45,45/) ; north
 glon2 = (/120,120,105,105,120/)
 glat2 = (/25,45,45,25,25/) ; east
 glon3 = (/110,110,75,75,110/)
 glat3 = (/10,25,25,10,10/) ; south
 glon4 = (/70,70,30,30,70/)
 glat4 = (/30,55,55,30,30/) ; west
;-----------------------------------


do i = 0, 3
    res@gsnLeftString = season+" trend ("+LEVs(i)+"-"+LEVe(i)+"hPa)"
    plot(i) = gsn_csm_contour_map(wks, trend(i,:,:), res)
    plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(i), plot_topo(i))

    plot_prob(i) = gsn_csm_contour(wks,prob(i,:,:),res3)
    plot_prob(i) = ShadeLtContour(plot_prob(i),sig,17)
    overlay(plot(i), plot_prob(i))

    ; dum1(i) = gsn_add_polyline(wks, plot(i), glon1, glat1, gres)
    ; dum2(i) = gsn_add_polyline(wks, plot(i), glon2, glat2, gres)
    ; dum3(i) = gsn_add_polyline(wks, plot(i), glon3, glat3, gres)
    ; dum4(i) = gsn_add_polyline(wks, plot(i), glon4, glat4, gres)
end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
; pres@pmLabelBarOrthogonalPosF = 0.12
pres@lbLabelFontHeightF = 0.015

gsn_panel(wks, plot, (/2,2/), pres)

end