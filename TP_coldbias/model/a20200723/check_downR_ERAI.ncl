
begin
yrStrt = 1979
yrLast = 2016
season = "JJA"
var = (/"downR","upR","netR"/)

diri = "/home/yangsong3/data-observation/wuyt/ERAI/monthly/1x1/"
fili1 = "net_thermal_sur_mon_7916_1.nc"
fili2 = "surface_thermal_downwards_monthly_1979_2018.nc"

f1 = addfile(diri+fili1, "r")
f2 = addfile(diri+fili2, "r")

TIME = f1->time
YYYY = cd_calendar(TIME, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
netR = f1->str(iYYYY,:,:)
downR = f2->strd(iYYYY,:,:)
upR = downR + netR
copy_VarCoords(downR, upR)
printVarSummary(upR)

x = new((/456,3,181,360/), float)
x(:,0,:,:) = downR
x(:,1,:,:) = upR
x(:,2,:,:) = netR
copy_VarCoords(netR, x(:,0,:,:))
x!1 = "Var"
x&Var = var
printVarSummary(x)

x_Sea = month_to_season(x, season)
x_sea = dim_avg_n_Wrap(x_Sea, 0)
printVarSummary(x_sea)


;=====================================================================
; figure
;=====================================================================

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/a20200723/"
figname = "check_downR_ERAI"
wks = gsn_open_wks("pdf", figpath+figname)

plot = new(3, graphic)

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
res@mpMaxLatF = 30
res@mpMinLatF = -30

res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False

res@gsnStringFontHeightF = 0.035
res@tmXBLabelFontHeightF = 0.02
res@tmYLLabelFontHeightF = 0.02

res@pmLabelBarOrthogonalPosF = 0.3
;res@lbLabelBarOn = False
; res@cnFillPalette = "CBR_coldhot"

res@gsnRightString = season

res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 10
res@cnMaxLevelValF = 450
res@cnMinLevelValF = 350
res@gsnLeftString = var(0)
plot(0) = gsn_csm_contour_map(wks, x_sea(0,:,:), res)

res@cnLevelSpacingF = 10
res@cnMaxLevelValF = 400
res@cnMinLevelValF = 300
res@gsnLeftString = var(1)
plot(1) = gsn_csm_contour_map(wks, x_sea(1,:,:), res)

res@cnLevelSpacingF = 5
res@cnMaxLevelValF = -30
res@cnMinLevelValF = -80
res@gsnLeftString = var(2)
plot(2) = gsn_csm_contour_map(wks, x_sea(2,:,:), res)

pres = True
pres@gsnMaximize = True
; pres@gsnPanelLabelBar = True
gsn_panel(wks, plot, (/3,1/), pres)




end