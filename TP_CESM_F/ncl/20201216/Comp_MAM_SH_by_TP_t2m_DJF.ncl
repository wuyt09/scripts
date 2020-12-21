
begin

yrStrt = 1979
yrLast = 2018
mnStrt = 1
mnLast = 12

year = ispan(yrStrt, yrLast, 1)
month = ispan(mnStrt, mnLast, 1)
penta = ispan(1, 73, 1)
nyear = dimsizes(year)
nmonth = dimsizes(month)
npenta = dimsizes(penta)
season = (/"MAM"/)

fx = "time_series_TP_t2m_DJF_MAM.txt"
x = asciiread(fx, (/nyear,3/), "float")
inx = x(:,1)

hig_ind = ind(inx.ge.1.)
low_ind = ind(inx.le.-1.)
print(year(hig_ind))
print(year(low_ind))


diri = "/home/yangsong3/data-observation/ERA5-monthly/surface/"
fili = "sensible_heat_flux_mon_1x1_1979_2020.nc"
f = addfile(diri+fili, "r")

Time = f->time
YYYY = cd_calendar(Time, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
var = short2flt(f->sshf(iYYYY,0,:,:))
var = var/86400.
var@units = "W/m~s~2~N~"
var = var * (-1.)
printVarSummary(var)

var_sea = month_to_season(var, season)
printVarSummary(var_sea)

var_hig = var_sea(hig_ind,:,:)
var_low = var_sea(low_ind,:,:)

;=============ttest====================================

var_hig_ave = dim_avg_n_Wrap(var_hig, 0)
var_low_ave = dim_avg_n_Wrap(var_low, 0)

var_hig_var = dim_variance_n_Wrap(var_hig, 0)
var_low_var = dim_variance_n_Wrap(var_low, 0)
printVarSummary(var_hig_ave)

prob = var_hig_ave
prob = ttest(var_hig_ave, var_hig_var, dimsizes(hig_ind), var_low_ave, var_low_var, dimsizes(low_ind), False, 0)
printVarSummary(prob)

diff_var = var_hig_ave
diff_var = var_hig_ave - var_low_ave
printVarSummary(diff_var)

;===================================================================
; plot
;===================================================================
lonstr = 60
lonend = 110
latstr = 15
latend = 45

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20201216/"
figname = "Comp_MAM_SH_by_TP_t2m_DJF"
wks = gsn_open_wks("pdf", figpath+figname)

plot = new(3, graphic)


res = True
res@gsnDraw = False
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
; res@cnFillMode = "CellFill"

res@gsnStringFontHeightF = 0.025
res@tmXBLabelFontHeightF = 0.02
res@tmYLLabelFontHeightF = 0.02

res@cnFillPalette = "cmp_flux"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -100
res@cnMaxLevelValF = 100
res@cnLevelSpacingF = 10
; res@lbLabelBarOn = False

res@gsnLeftString = "(a) High DJF TP t2m"
plot(0) = gsn_csm_contour_map(wks, var_hig_ave, res)

res@gsnLeftString = "(b) Low DJF TP t2m"
plot(1) = gsn_csm_contour_map(wks, var_low_ave, res)

res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -10
res@cnMaxLevelValF = 10
res@cnLevelSpacingF = 1
res@lbLabelBarOn = True
res@gsnLeftString = "(a) - (b)"
plot(2) = gsn_csm_contour_map(wks, diff_var, res)

pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
pres@pmLabelBarOrthogonalPosF = -0.015
gsn_panel(wks, plot, (/3,1/), pres)


end