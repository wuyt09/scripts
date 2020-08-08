
begin
yrStrt = 1979
yrLast = 2015
year = ispan(yrStrt, yrLast, 1)
nyear = dimsizes(year)

monsoon = (/"BOB","SCS","ISM"/)
month = "May"

lonstr = 60
lonend = 110
latstr = 20
latend = 45
;===============================================================
; read onset pentad
;===============================================================
X = asciiread("./onset_pentad_ERAI.txt", (/nyear,4/), "float")
printVarSummary(X)
x = X(:,1:3)

print(x)

;---------------------------------------------------------
;topo data
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)
;---------------------------------------------------------

;=========================================================
; Q1
;=========================================================

diri = "/home/yangsong3/wuyt/data/ERA-Interim/monthly/"
fili = "q1.mon.mean.1979.2016.ERAI.nc"
f = addfile(diri+fili, "r")

TIME = f->time
YYYY = cd_calendar(TIME, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
q1 = f->q1(iYYYY,:,::-1,:)
printVarSummary(q1)

lat = q1&lat
lon = q1&lon
lev = q1&lev
nlat = dimsizes(lat)
nlon = dimsizes(lon)
nlev = dimsizes(lev)

fsp = addfile("/home/yangsong3/data-observation/ERA-interim/monthly/Surface_pressure.mon.mean.nc", "r")
PS = short2flt(fsp->sp(iYYYY,::-1,:))
ps = linint2_Wrap(PS&longitude, PS&latitude, PS, True, lon, lat, 0)
printVarSummary(ps)

;-----------------------------------------
; vertical integration
;-----------------------------------------

plev = lev
plev = plev*100
ptop = min(plev)
dp = dpres_plevel_Wrap(plev, ps, ptop, 0)
q1dp = q1*dp
copy_VarCoords(q1, q1dp)
Q1 = dim_sum_n_Wrap(q1dp, 1)
printVarSummary(Q1)

Q1_mon = Q1(4::12,:,:)
printVarSummary(Q1_mon)


;============================================================
; Correlation
;============================================================
cc = new((/3,nlat,nlon/), float)
prob = cc
do i = 0, 2
    cc(i,:,:) = escorc_n(x(:,i), Q1_mon, 0, 0)
    prob(i,:,:) = rtest(cc(i,:,:), nyear, 0)
end do
copy_VarCoords(Q1_mon(0,:,:), cc(0,:,:))
copy_VarCoords(cc, prob)
printVarSummary(cc)
printVarSummary(prob)

;============================================================
; figure
;============================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20200704/"
figname = "Corr_onset_pentad_with_TP_Q1_"+month+"_ERAI"

wks = gsn_open_wks("pdf", figpath+figname)
plot = new(3, graphic)
plot_topo = plot
plot_prob = plot
;----------------------------------------------------
 res_topo            = True
 res_topo@cnFillOn   = False
 res_topo@cnLinesOn  = True
 res_topo@gsnFrame   = False
 res_topo@gsnDraw    = False
 res_topo@cnLineLabelsOn         = False
 res_topo@cnLineThicknessF       = 4
 res_topo@cnLineDashPattern      = 16
 res_topo@cnLineColor            = "blue" ; gray25
 res_topo@cnInfoLabelOn          = False
 res_topo@cnLevelSelectionMode   ="ExplicitLevels"
 res_topo@cnLevels               = (/1500/)
 res_topo@gsnLeftString = ""
 res_topo@gsnRightString = ""
;----------------------------------------------------

res = True
res@gsnDraw  = False
res@gsnFrame = False
res@gsnMaximize = True

res@mpOutlineOn = True
res@mpFillOn = True
res@mpGeophysicalLineThicknessF = 2
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

res@gsnStringFontHeightF = 0.035
res@tmXBLabelFontHeightF = 0.02
res@tmYLLabelFontHeightF = 0.02

res@lbLabelBarOn = False
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 0.1
res@cnMaxLevelValF = 0.6
res@cnMinLevelValF = -0.6
res@cnFillPalette = "CBR_coldhot"
;res@gsnRightString = "K"

;------------------------
; draw prob
;------------------------
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

;res3@cnFillPattern = 17 ;(/17/)
res3@cnFillColor  =  "black" ;"green4"
res3@gsnLeftString = ""
res3@cnFillDotSizeF = 0.003
res3@cnFillScaleF  = 1 ;1.2

sig = 0.95
res@gsnRightString = sig*100+"%"

do i = 0, 2
    res@gsnLeftString = monsoon(i)
    plot(i) = gsn_csm_contour_map(wks, cc(i,:,:), res)
    plot_prob(i) = gsn_csm_contour(wks, prob(i,:,:), res3)
    plot_prob(i) = ShadeLtContour(plot_prob(i),(1-sig),6)
    overlay(plot(i), plot_prob(i))
    plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(i), plot_topo(i))
end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = "Correlation coefficient of "+month+" mean: onset time & TP Q1"
gsn_panel(wks, plot, (/3,1/), pres)

end