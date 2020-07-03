
begin

lonstr = 75
lonend = 105
latstr = 25
latend = 40
yrStrt = 198102
yrLast = 200601


casename = (/"Ctrl_FAMIP","Hist_FAMIP"/) ;,"Hist_IO","Hist_PO","Hist_AO","Hist_PIO","Hist_PAO","Hist_IAO"

ncase = dimsizes(casename)
season = (/"MAM","JJA","SON","DJF"/)
season_new = (/"Ann","MAM","JJA","SON","DJF"/)
year = ispan(1981, 2005, 1)
varname = (/"netR","downR","Ts"/)

var = new((/ncase,5,25,3,96,144/), float)
do i = 0, ncase-1
    diri = "/home/yangsong3/data-model/wuyt/TPbias/TPbias_"+casename(i)+"/a20191206/"
    fili1 = casename(i)+".cam.h0.FLNS.1979-2005.nc"
    fili2 = casename(i)+".cam.h0.FLDS.1979-2005.nc"
    fili3 = casename(i)+".cam.h0.TREFHT.1979-2005.nc"

    f1 = addfile(diri+fili1, "r")
    f2 = addfile(diri+fili2, "r")
    f3 = addfile(diri+fili3, "r")
    TIME = f1->time
    YYYY = cd_calendar(TIME, -1)
    iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)\

    x1 = f1->FLNS(iYYYY,:,:)
    x2 = f2->FLDS(iYYYY,:,:)
    x3 = f3->TREFHT(iYYYY,:,:)
    printVarSummary(x1)
    printVarSummary(x2)
    printVarSummary(x3)

    x = new((/300,3,96,144/), float)
    x(:,0,:,:) = x1
    x(:,1,:,:) = x2
    x(:,2,:,:) = x3
    x!1 = "var"
    x&var = varname
    copy_VarCoords(x1, x(:,0,:,:))

    y = new((/5,25,3,96,144/), float)
    y(0,:,:,:,:) = month_to_annual(x, 1)
    y(1:4,:,:,:,:) = month_to_seasonN(x, season)

    var(i,:,:,:,:,:) = y
    copy_VarCoords(x1(0,:,:), var(0,0,0,0,:,:))
end do

var!0 = "case"
var!1 = "season"
var!2 = "year"
var!3 = "var"
var&case = casename
var&season = season_new
var&year = year
var&var = varname

; var = dim_standardize_n_Wrap(var, 0, 2)
printVarSummary(var)

;=========================================================
; outline of tp 
;=========================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)
TOPO = linint2_Wrap(topo&lon, topo&lat, topo, True, var&lon, var&lat, 0)
var = mask(var, TOPO.ge.1500, True)
var_area = wgt_areaave_Wrap(var(:,:,:,:,{latstr:latend},{lonstr:lonend}), 1, 1, 0)
printVarSummary(var_area)


pp = new((/ncase+1,5,25,3/), float)
pp(0:1,:,:,:) = var_area
pp(2,:,:,:) = var_area(1,:,:,:) - var_area(0,:,:,:)
copy_VarCoords(var_area(0,:,:,:), pp(0,:,:,:))

cc_12 = new((/ncase+1,5/), float)
cc_13 = new((/ncase+1,5/), float)
do i = 0, ncase
    do j = 0, 4
        cc_12(i,j) = escorc(pp(i,j,:,0), pp(i,j,:,2))        
        cc_13(i,j) = escorc(pp(i,j,:,1), pp(i,j,:,2))
    end do
end do
copy_VarCoords(pp(:,:,0,0), cc_12)
copy_VarCoords(pp(:,:,0,0), cc_13)
write_table("./corr_netR_Tas.txt", "w", [/casename(0), casename(1), "Diff"/], "%10s %10s  %10s")
write_table("./corr_netR_Tas.txt", "a", [/cc_12(0,:),cc_12(1,:),cc_12(2,:)/], "%10.2f %10.2f %10.2f")

write_table("./corr_downR_Tas.txt", "w", [/casename(0), casename(1), "Diff"/], "%10s %10s %10s")
write_table("./corr_downR_Tas.txt", "a", [/cc_13(0,:),cc_13(1,:),cc_13(2,:)/], "%10.2f %10.2f %10.2f")
exit
;=========================================================
; plot 
;=========================================================
do j = 0, 2
title = (/"Ctrl","Hist","Diff"/)
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/a20200605/"
figname = "Ann_4season_netR_downR_Ts_"+title(j)

Colors = (/"green4","purple","red"/) ;,"","goldenrod3","black","orange2"
leftsti = season_new
Texts  = varname


wks = gsn_open_wks("pdf", figpath+figname)

plot = new(5, graphic)

res = True
res@gsnDraw  = False
res@gsnFrame = False
res@gsnMaximize = True

res@vpWidthF = 0.7
res@vpHeightF = 0.4

res@trXMinF = min(year)
res@trXMaxF = max(year)
; res@tmXBMode = "Explicit"
; res@tmXBValues = Month
; res@tmXBLabels = month
res@tmXTOn = False
; res@tiXAxisString = "Month"
res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025
res@gsnStringFontHeightF = 0.035

; res@trYMinF = -2
; res@trYMaxF = 2
; res@tmYLMode = "Explicit"
; res@tmYLValues = ispan(-, 0, 1)
; res@tmYLLabels = ispan(-5, 0, 1)
res@tmYROn = False

; res@gsnYRefLine = 0
; res@gsnYRefLineDashPattern = 16
; res@gsnYRefLineColor = "gray20"
; res@gsnYRefLineThicknessF = 5

res@tmXMajorGrid                = True
res@tmYMajorGrid                = True
res@tmXMajorGridThicknessF      = 1.0
res@tmYMajorGridThicknessF      = 1.0
res@tmXMajorGridLineDashPattern = 2
res@tmYMajorGridLineDashPattern = 2

res@xyLineThicknessF = 5
res@xyDashPatterns = 0
res@xyLineColor = Colors(0)
res@tiYAxisString = "netR, W/m~S~2~C~"
res@tiYAxisFontColor = res@xyLineColor
res@tmYLLabelFontColor = res@xyLineColor

res2 = True
res2@xyLineColor = Colors(1)
res2@xyLineThicknessF = 5
res2@tiYAxisString = "downR, W/m~S~2~C~"
res2@tiYAxisFontColor = res2@xyLineColor
res2@tmYRLabelFontColor = res2@xyLineColor


res3 = res2
res3@xyLineColor = Colors(2)
res3@tiYAxisString = "Ts, K"
res3@tiYAxisFontColor = res3@xyLineColor
res3@tmYRLabelFontColor = res3@xyLineColor

;==================================================================================

do i = 0, 4
    res@gsnLeftString = leftsti(i)
    plot(i) = gsn_csm_xy3(wks, year, pp(j,i,:,0), pp(j,i,:,1), pp(j,i,:,2), res, res2, res3)
end do

; ;==================================================================================
;     res_lines                   = True                  ; polyline mods desired
;     res_lines@gsLineThicknessF  = 5                     ; line thicker
    

;     res_text                    = True                  ; text mods desired
;     res_text@txFontHeightF      = 0.028                 ; change text size
;     res_text@txJust             = "CenterLeft"          ; text justification



; ; ;============================================================
;     interval = 0.4
;     gx = (/1982,1984/)
;     text_pos = 1985
;     gy1 = (/1.8,1.8/)
;     gy2 = gy1-interval ;(/15,15/)
;     gy3 = gy2-interval
;     gy4 = gy1 ;gy3-interval
;     gy5 = gy2 ;gy4-interval
;     gy6 = gy3 ;gy5-interval
;     print(gy1)

;     res_lines@gsLineColor = Colors(0)
;     ;res_lines@gsLineDashPattern = Patterns(0)
;     dum11 = gsn_add_polyline(wks, plot(0), gx, gy1, res_lines)
;     dum12 = gsn_add_text(wks, plot(0), Texts(0), text_pos, gy1, res_text)
;     res_lines@gsLineColor = Colors(1)
;     ;res_lines@gsLineDashPattern = Patterns(1)
;     dum21 = gsn_add_polyline(wks, plot(0), gx, gy2, res_lines)
;     dum22 = gsn_add_text(wks, plot(0), Texts(1), text_pos, gy2, res_text)
;     res_lines@gsLineColor = Colors(2)
;     ;res_lines@gsLineDashPattern = Patterns(2)
;     dum31 = gsn_add_polyline(wks, plot(0), gx, gy3, res_lines)
;     dum32 = gsn_add_text(wks, plot(0), Texts(2), text_pos, gy3, res_text)


pres = True
pres@gsnMaximize = True
pres@gsnPanelMainString = title(j)
;pres@pmLabelBarOrthogonalPosF = -0.02
;pres@lbLabelFontHeightF = 0.015
gsn_panel(wks, plot, (/3,2/), pres)

end do
end