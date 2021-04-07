
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
ts_obs = new((/nyear*nm,163,283/), float,-99.)

diriobs = "/home/yangsong3/data-observation/CN05.1/"

filit = "CN05.1_Tm_1961_2018_month_025x025.nc"
filip = "CN05.1_Pre_1961_2018_month_025x025.nc"
ft = addfile(diriobs+filit, "r")
fp = addfile(diriobs+filip, "r")

Time = ft->time
YYYY = ut_calendar(Time, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
t_time = ft->time(iYYYY)
ts_obs = ft->tm(iYYYY,:,:)
pr_obs = fp->pre(iYYYY,:,:)

ts_obs = where(ismissing(ts_obs), 0, ts_obs)
pr_obs = where(ismissing(pr_obs), 0, pr_obs)
ts_obs = dtrend_msg_n(t_time, ts_obs, False, False, 0)
pr_obs = dtrend_msg_n(t_time, pr_obs, False, False, 0)

printVarSummary(ts_obs)
; exit()
; lat = fspan(14.75, 55.25, 163)
; lon = fspan(69.75, 140.25, 283)
; print(lat)
; print(lon)

; ts_obs!0 = "time"
; ts_obs!1 = "lat"
; ts_obs!2 = "lon"
; ts_obs&time = ispan(1, nyear*nm, 1)
; ts_obs&lat = lat
; ts_obs&lon = lon
; ts_obs&lat@units = "degrees_north"
; ts_obs&lon@units = "degrees_east"
; printVarSummary(ts_obs)

; TOPO = topo({14.75:55.25},{69.25:140.25})
; TOPOx = linint2_Wrap(TOPO&lon, TOPO&lat, TOPO, False, fspan(69.75, 140.25, 283), fspan(14.75, 55.25, 163), 0)
; ts_obs = mask(ts_obs, TOPOx.ge.1500, True)
; copy_VarCoords(ts_obs, TP_obs)
TP_obs_sea = ts_obs(4::12,{latstr:latend},{lonstr:lonend})
Left = reshape(TP_obs_sea, (/nyear,dimsizes(TP_obs_sea&lon)*dimsizes(TP_obs_sea&lat)/))

Left!0 = "time"
Left!1 = "col"
printVarSummary(Left)
;===================================================================

prec_sea = pr_obs(5::12,:,:)
prec_sea_EA = prec_sea(:,{18:45},{100:130})
; prec_sea_EA = prec_sea(:,{22:35},{106:120})
Right = reshape(prec_sea_EA, (/nyear,dimsizes(prec_sea_EA&lat)*dimsizes(prec_sea_EA&lon)/))

Right!0 = "time"
Right!1 = "col"
printVarSummary(Right)  
;===================================================================
; check
; p2d = Right ;(col|:,time|:)
; p2d_dims = dimsizes(p2d) ; I assume this is time x p1d
; ntim = p2d_dims(0)
; np1d = p2d_dims(1)
; do i=0,np1d-1
;   pmin = min(p2d(:,i))
;   pmax = max(p2d(:,i))
;   print(pmin+" "+pmax)
;   if(ismissing(pmin) .or. pmin.eq.pmax) then
;     print("p2d contains all constant values at np1d = " + i)
;   end if
; end do
; exit
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

xL = reshape(svLeft, (/dimsizes(TP_obs_sea&lat),dimsizes(TP_obs_sea&lon)/))
xR = reshape(svRight, (/dimsizes(prec_sea_EA&lat),dimsizes(prec_sea_EA&lon)/))
copy_VarCoords(TP_obs_sea(0,:,:), xL)
copy_VarCoords(prec_sea_EA(0,:,:), xR)
printVarSummary(xL)


lines = new((/2,nyear/), float)
lines(0,:) = svd@ak
lines(1,:) = svd@bk
corr = escorc(lines(0,:), lines(1,:))
;======================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_surf_temp_EA_prec/20210316/"
figname = "SVD_May_TP_t2m_June_EA_prec_CN051"
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
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 0.01
res@cnMaxLevelValF = 0.05
res@cnMinLevelValF = -0.05
plot(0) = gsn_csm_contour_map(wks, xL, res)
plot_topo(0) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(0), plot_topo(0))

res@mpMinLatF = 22
res@mpMaxLatF = 40
res@mpMinLonF = 106
res@mpMaxLonF = 122
res@cnFillPalette = "CBR_drywet"
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 0.05
res@cnMaxLevelValF = 0.2
res@cnMinLevelValF = -0.2
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