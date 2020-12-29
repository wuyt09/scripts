
begin


lonstr = 75 ; 72 ;75
lonend = 105 ;136 ;105
latstr = 25 ;18 ;25
latend = 40 ;54 ; 40


;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)
printVarSummary(topo)

;===================================================================
; OBS grid data
;===================================================================
yrStrt = 1979
yrLast = 2018
year = ispan(yrStrt, yrLast, 1)
nyear = dimsizes(year)
month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
nm = dimsizes(month)

;===================================================================
diri = "/home/yangsong3/data-observation/ERA5-monthly/surface/"
fili = "t2m_mon_1x1_1979_2019.nc"
f1 = addfile(diri+fili, "r")

Time = f1->time
YYYY = cd_calendar(Time, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
time = f1->time(iYYYY)
ts_obs = short2flt(f1->t2m(iYYYY,::-1,:))
printVarSummary(ts_obs)
nlon = dimsizes(ts_obs&longitude)
nlat = dimsizes(ts_obs&latitude)

ts_obs = dtrend_msg_n(time, ts_obs, False, False, 0)
; TP_obs = mask(ts_obs, topo.ge.1500, True)
; copy_VarCoords(ts_obs, TP_obs)
TP_obs_sea = ts_obs(4::12,{latstr:latend},{lonstr:lonend})
Left = reshape(TP_obs_sea, (/nyear,dimsizes(TP_obs_sea&longitude)*dimsizes(TP_obs_sea&latitude)/))

Left!0 = "time"
Left!1 = "col"
printVarSummary(Left)
;===================================================================

diri2 = "/home/yangsong3/data-observation/ERA5-monthly/surface/"
filip = "tp_mon_1x1_1979_2019.nc"

fp = addfile(diri2+filip, "r")
Time = fp->time
YYYY = cd_calendar(Time, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
prec = short2flt(fp->tp(iYYYY,::-1,:))
prec = prec*1000
prec@units = "mm/day"

printVarSummary(prec)

prec_sea = prec(5::12,:,:)
prec_sea_EA = prec_sea(:,{18:45},{100:130})
Right = reshape(prec_sea_EA, (/nyear,dimsizes(prec_sea_EA&latitude)*dimsizes(prec_sea_EA&longitude)/))

Right!0 = "time"
Right!1 = "col"
printVarSummary(Right)
;===================================================================

nsvd = 1
homlft = new((/nsvd,dimsizes(Left(0,:))/), float)
hetlft = new((/nsvd,dimsizes(Left(0,:))/), float)
homrgt = new((/nsvd,dimsizes(Right(0,:))/), float)
hetrgt = new((/nsvd,dimsizes(Right(0,:))/), float)
svd = svdstd(Left(col|:,time|:), Right(col|:,time|:), nsvd, homlft, hetlft, homrgt, hetrgt)
print("svdstd: percent variance= "+svd)
printVarSummary(svd)

svLeft = new((/nsvd,dimsizes(Left(0,:))/), float)
svRight = new((/nsvd,dimsizes(Right(0,:))/), float)
pcVar = svdstd_sv(Left(col|:,time|:), Right(col|:,time|:), nsvd, svLeft, svRight)
printVarSummary(pcVar)
print("svdstd_sv: percent variance= "+pcVar)
; print("svdstd_sv: singular values= "+pcVar@sv)
printVarSummary(svLeft)
printVarSummary(svRight)

xL = reshape(svLeft, (/dimsizes(TP_obs_sea&latitude),dimsizes(TP_obs_sea&longitude)/))
xR = reshape(svRight, (/dimsizes(prec_sea_EA&latitude),dimsizes(prec_sea_EA&longitude)/))
copy_VarCoords(TP_obs_sea(0,:,:), xL)
copy_VarCoords(prec_sea_EA(0,:,:), xR)
printVarSummary(xL)

lines = new((/2,nyear/), float)
lines(0,:) = svd@ak
lines(1,:) = svd@bk
corr = escorc(lines(0,:), lines(1,:))
;======================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20201216/"
figname = "SVD_May_TP_t2m_June_EA_prec"
wks = gsn_open_wks("pdf", figpath+figname)

plot = new(2, graphic)
plot_topo = plot

 res_topo            = True
 res_topo@cnFillOn   = False
 res_topo@cnLinesOn  = True
 res_topo@gsnFrame   = False
 res_topo@gsnDraw    = False
 res_topo@cnLineLabelsOn         = False
 res_topo@cnLineThicknessF       = 3
 res_topo@cnLineColor            = "grey40"
 res_topo@cnInfoLabelOn          = False
 res_topo@cnLevelSelectionMode   ="ExplicitLevels"
 res_topo@cnLevels               = (/2000/)
 res_topo@gsnLeftString = ""
 res_topo@gsnRightString = ""



res = True
res@gsnDraw = False
res@gsnFrame = False
; res@gsnMaximize = True
res@gsnAddCyclic = False

res@mpOutlineOn = True
res@mpFillOn = True
res@mpGeophysicalLineThicknessF = 2
res@mpGeophysicalLineColor = "grey25" ;"red" 
; res@mpFillDrawOrder = "PostDraw"
; res@mpCenterLonF = 180

res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False
; res@cnFillMode = "CellFill"

res@gsnStringFontHeightF = 0.025
res@tmXBLabelFontHeightF = 0.02
res@tmYLLabelFontHeightF = 0.02

res@mpMinLatF = latstr
res@mpMaxLatF = latend
res@mpMinLonF = lonstr
res@mpMaxLonF = lonend
res@cnFillPalette = "CBR_coldhot"
; res@cnLevelSelectionMode = "ManualLevels"
; res@cnLevelSpacingF = 0.02
; res@cnMaxLevelValF = 0.1
; res@cnMinLevelValF = -0.1
res@gsnLeftString = "May TP t2m"
res@gsnRightString = sprintf("%4.2f", svd)+"%"
plot(0) = gsn_csm_contour_map(wks, xL, res)
plot_topo(0) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(0), plot_topo(0))

res@mpMinLatF = 22
res@mpMaxLatF = 40
res@mpMinLonF = 106
res@mpMaxLonF = 122
res@cnFillPalette = "CBR_drywet"
; res@cnLevelSelectionMode = "ManualLevels"
; res@cnLevelSpacingF = 0.02
; res@cnMaxLevelValF = 0.1
; res@cnMinLevelValF = -0.1
res@gsnLeftString = "June EA prec"
plot(1) = gsn_csm_contour_map(wks, xR, res)


pres = True
pres@gsnMaximize = True
gsn_panel(wks, plot, (/1,2/), pres)

delete(res)
Colors = (/"black","red","blue","black","red","blue"/)
Patterns = (/0,0,0,2,2,2/)

res = True
res@gsnFrame = False
res@gsnDraw  = False
res@gsnMaximize = True

res@vpWidthF = 0.7
res@vpHeightF = 0.45

res@xyLineThicknessF = 10
res@xyLineColors = Colors
res@xyDashPatterns = Patterns

res@trXMinF = yrStrt
res@trXMaxF = yrLast

res@gsnYRefLine = (/0/)
res@gsnYRefLineColor = "grey40"
res@gsnYRefLineThicknessF = 5
res@gsnYRefLineDashPattern = 2

res@gsnRightString = "Corr: "+sprintf("%4.2f", corr)
res@gsnStringFontHeightF = 0.025

plot(0) = gsn_csm_xy(wks, year, lines, res)

;==================================================================================
    res_lines                   = True                  ; polyline mods desired
    res_lines@gsLineThicknessF  = 5                     ; line thicker
    

    res_text                    = True                  ; text mods desired
    res_text@txFontHeightF      = 0.018                 ; change text size
    res_text@txJust             = "CenterLeft"          ; text justification

    Texts  = (/"TP t2m","EA prec"/)

;============================================================
    gx = (/1980,1982/)
    text_pos = 1983
    gy1 = (/28,28/)
    interval = 5
    gy2 = gy1 - interval
    gy3 = gy2 - interval

    res_lines@gsLineColor = Colors(0)
    res_lines@gsLineDashPattern = 0
    dum11 = gsn_add_polyline(wks, plot(0), gx, gy1, res_lines)
    dum12 = gsn_add_text(wks, plot(0), Texts(0), text_pos, gy1, res_text)
    res_lines@gsLineColor = Colors(1)
    ;res_lines@gsLineDashPattern = 2
    dum21 = gsn_add_polyline(wks, plot(0), gx, gy2, res_lines)
    dum22 = gsn_add_text(wks, plot(0), Texts(1), text_pos, gy2, res_text)

    draw(plot(0))
end