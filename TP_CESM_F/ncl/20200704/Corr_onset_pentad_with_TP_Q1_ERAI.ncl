
begin
yrStrt = 1979
yrLast = 2015
year = ispan(yrStrt, yrLast, 1)
nyear = dimsizes(year)

monsoon = (/"BOB","SCS","ISM"/)

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

Q1_mon = Q1(2::12,:,:)
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
copy_VarCoords(Q1_mon, cc(0,:,:))
copy_VarCoords(cc, prob)
printVarSummary(cc)
printVarSummary(prob)

;============================================================
; figure
;============================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20200704/"
figname = "Corr_onset_pentad_with_TP_Q1_March_ERAI"

wks = gsn_open_wks("pef", figpath+figname)
plot = new(3, graphic)

plot_prob = plot

res = True
res@gsnDraw  = False
res@gsnFrame = False
res@gsnMaximize = True

res@mpOutlineOn = True
res@mpFillOn = True
res@mpGeophysicalLineThicknessF = 2
res@mpGeophysicalLineColor = "grey25" ;"red"
res@mpFillDrawOrder = "PostDraw"
res@mpCenterLonF = 180

res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False

res@gsnStringFontHeightF = 0.035
res@tmXBLabelFontHeightF = 0.02
res@tmYLLabelFontHeightF = 0.02

; res@lbLabelBarOn = False
; res@cnLevelSelectionMode = "ManualLevels"
; res@cnLevelSpacingF = 0.1
; res@cnMaxLevelValF = 0.8
; res@cnMinLevelValF = -0.8
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
    res@gsnLeftString = models(ip_tbias(i))
    res@gsnCenterString = sprintf("%4.2f", t_bias(ip_tbias(i)))
    plot(i) = gsn_csm_contour_map(wks, cc(ip_tbias(i),:,:), res)
    plot_prob(i) = gsn_csm_contour(wks, prob(ip_tbias(i),:,:), res3)
    plot_prob(i) = ShadeLtContour(plot_prob(i),(1-sig),6)
    overlay(plot(i), plot_prob(i))
end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = "Correlation coefficient of "+season+" mean: model Ts & model SST "
gsn_panel(wks, plot, (/7, 4/), pres)

end