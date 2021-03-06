begin

lonstr = 0
lonend = 130
latstr = 10
latend = 60

yrStrt = 1979
yrLast = 2019
year = ispan(yrStrt, yrLast, 1)
nyear = dimsizes(year)

;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)
printVarSummary(topo)


;===================================================================
diri = "/home/yangsong3/data-observation/ERA5-monthly/surface/"
fili = "t2m_mon_1x1_1979_2019.nc"
; fili = "skt_mon_1x1_1979_2019.nc"
f = addfile(diri+fili, "r")

Time = f->time
YYYY = cd_calendar(Time, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

t2m = short2flt(f->t2m(iYYYY,::-1,:))
printVarSummary(t2m)


t2m_clm = month_to_annual(t2m, 1)
printVarSummary(t2m_clm)

rc = regCoef_n(year, t2m_clm, 0, 0)
copy_VarCoords(t2m_clm(0,:,:), rc)
printVarSummary(rc)

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
;===============================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_surrounding/"
figname = "trend_ann_surfT_ERA5"
wks = gsn_open_wks("pdf", figpath+figname)

plot = new(1, graphic)
plot_topo = plot
plot_prob = new(1, graphic)

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
res@cnMinLevelValF = -1
res@cnMaxLevelValF = 1
res@cnLevelSpacingF = 0.2
; res@lbLabelBarOn = False

res@gsnLeftString = "Annual trend ("+yrStrt+"-"+yrLast+")"
res@pmLabelBarOrthogonalPosF = 0.12
res@lbLabelFontHeightF = 0.015
; res@tiMainString = "surface temperature"

plot(0) = gsn_csm_contour_map(wks, trend, res)
plot_topo(0) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(0), plot_topo(0))

plot_prob(0) = gsn_csm_contour(wks,prob,res3)
plot_prob(0) = ShadeLtContour(plot_prob(0),sig,17)
overlay(plot(0), plot_prob(0))


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


dum1 = gsn_add_polyline(wks, plot(0), glon1, glat1, gres)
dum2 = gsn_add_polyline(wks, plot(0), glon2, glat2, gres)
dum3 = gsn_add_polyline(wks, plot(0), glon3, glat3, gres)
dum4 = gsn_add_polyline(wks, plot(0), glon4, glat4, gres)
draw(plot(0))

end