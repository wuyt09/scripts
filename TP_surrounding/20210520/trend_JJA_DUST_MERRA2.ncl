
begin

lonstr = -20
lonend = 160
latstr = 0
latend = 70

yrStrt = 1980
yrLast = 2019
year = ispan(yrStrt, yrLast, 1)
nyear = dimsizes(year)
month = ispan(1, 12, 1)


season = "JJA"

;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)
printVarSummary(topo)
;===================================================================

diri = "/home/yangsong3/data-observation/MERRA2/monthly/aer_AOT/"
x = new((/nyear*12,361,576/), float)


do i = 0, nyear-1
    do j = 0, 11
        nyr = year(i)
        nmn = month(j)
        fili = systemfunc("ls "+diri+"MERRA2_*.tavgM_2d_aer_Nx."+nyr+sprinti("%0.2i", nmn)+".nc4")

        f = addfile(fili, "r")

        x_med = f->DUEXTTAU(0,:,:)
        x(j+i*12,:,:) = x_med
    end do
end do
printVarSummary(x)
x!0 = "time"
x&time = ispan(1, nyear*12, 1)

t2m_clm = month_to_season(x, season)
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
trend@units = "1/decade"
printVarSummary(trend)

sig = 0.05

;===============================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_surrounding/20210520/"
figname = "trend_"+season+"_DUST_MERRA2"
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

; res@cnFillPalette = "CBR_coldhot"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -0.03
res@cnMaxLevelValF = 0.03
res@cnLevelSpacingF = 0.005
; res@lbLabelBarOn = False

res@gsnLeftString = "Dust "+season+" trend ("+yrStrt+"-"+yrLast+")"
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
 glon1 = (/120,120,85,85,120/)
 glat1 = (/45,60,60,45,45/) ; north
 glon2 = (/140,140,105,105,140/)
 glat2 = (/30,45,45,30,30/) ; east
 ; glon3 = (/100,100,50,50,100/)
 ; glat3 = (/10,25,25,10,10/) ; south
 glon3 = (/100,100,50,50,100/)
 glat3 = (/10,25,25,10,10/) ; south
 glon4 = (/60,60,-10,-10,60/)
 glat4 = (/30,55,55,30,30/) ; west
;-----------------------------------


dum1 = gsn_add_polyline(wks, plot(0), glon1, glat1, gres)
dum2 = gsn_add_polyline(wks, plot(0), glon2, glat2, gres)
dum3 = gsn_add_polyline(wks, plot(0), glon3, glat3, gres)
dum4 = gsn_add_polyline(wks, plot(0), glon4, glat4, gres)
draw(plot(0))
end