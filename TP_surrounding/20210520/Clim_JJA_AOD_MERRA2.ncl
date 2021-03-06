
begin

lonstr = -20
lonend = 160
latstr = 0
latend = 70

yrStrt = 1980
yrLast = 2019
year = ispan(yrStrt, yrLast, 1)

season = "JJA"

;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)
printVarSummary(topo)
;===================================================================
diri = "/home/yangsong3/data-observation/MERRA2/monthly/aod/"
fili = "aod.198001-202103.mon.nc4"

f = addfile(diri+fili, "r")

Time = f->time
YYYY = cd_calendar(Time, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
x = f->AODANA(iYYYY,:,:)
printVarSummary(x)

t2m_clmx = month_to_season(x, season)
t2m_clm = dim_avg_n_Wrap(t2m_clmx, 0)
printVarSummary(t2m_clm)

;===============================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_surrounding/20210520/"
figname = "Clim_"+season+"_AOD_MERRA2"
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

res@cnFillPalette = "WhiteBlue"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = 0.1
res@cnMaxLevelValF = 1
res@cnLevelSpacingF = 0.1
; res@lbLabelBarOn = False

res@gsnLeftString = "AOD "+season+" Clim ("+yrStrt+"-"+yrLast+")"
res@pmLabelBarOrthogonalPosF = 0.12
res@lbLabelFontHeightF = 0.015
; res@tiMainString = "surface temperature"

plot(0) = gsn_csm_contour_map(wks, t2m_clm, res)
plot_topo(0) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(0), plot_topo(0))

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
 glon3 = (/90,90,70,70,90/)
 glat3 = (/10,28,28,10,10/) ; south
 glon4 = (/60,60,-10,-10,60/)
 glat4 = (/30,55,55,30,30/) ; west
;-----------------------------------


dum1 = gsn_add_polyline(wks, plot(0), glon1, glat1, gres)
dum2 = gsn_add_polyline(wks, plot(0), glon2, glat2, gres)
dum3 = gsn_add_polyline(wks, plot(0), glon3, glat3, gres)
dum4 = gsn_add_polyline(wks, plot(0), glon4, glat4, gres)
draw(plot(0))
end