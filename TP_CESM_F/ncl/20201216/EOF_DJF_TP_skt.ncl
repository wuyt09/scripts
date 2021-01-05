
begin

lonstr = 75
lonend = 105
latstr = 25 ;-10
latend = 40

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
season = "DJF"

;===================================================================
diri = "/home/yangsong3/data-observation/ERA5-monthly/surface/"
fili = "skt_mon_1x1_1979_2019.nc"
f1 = addfile(diri+fili, "r")

Time = f1->time
YYYY = cd_calendar(Time, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
time = f1->time(iYYYY)
ts_obs = short2flt(f1->skt(iYYYY,::-1,:))
printVarSummary(ts_obs)
nlon = dimsizes(ts_obs&longitude)
nlat = dimsizes(ts_obs&latitude)

; ts_obs = mask(ts_obs, topo.ge.1500, True)
ts_obs = dtrend_msg_n(time, ts_obs, False, False, 0)
TS_obs_sea = month_to_season(ts_obs, season)
TP_obs_sea = TS_obs_sea(:,{latstr:latend},{lonstr:lonend})
printVarSummary(TP_obs_sea)

nlat = dimsizes(TP_obs_sea&latitude)
nlon = dimsizes(TP_obs_sea&longitude)
neval = 1
optEOF = True
t_eof = new((/neval,nlat,nlon/), float)
t_eof_ts = new((/neval,nyear/), float)

x = eofunc_n_Wrap(TP_obs_sea, neval, optEOF, 0)
t_eof = x
pcvar = x@pcvar
t_eof_ts = eofunc_ts_n_Wrap(TP_obs_sea, x, False, 0)
printVarSummary(t_eof)
printVarSummary(t_eof_ts)

std_val = dim_stddev_n_Wrap(t_eof_ts, 1)
print(std_val)

inx_hig = ind(t_eof_ts(0,:) .ge. std_val)
inx_low = ind(t_eof_ts(0,:) .le. -1.*std_val)
print(year(inx_hig))
print(year(inx_low))
;=================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20201216/"
figname = "EOF_"+season+"_TP_skt"

wks = gsn_open_wks("pdf", figpath+figname)
nplot = 2
plot = new(nplot, graphic)
plot_topo = new(1, graphic)

;----------------------------------------------------
 res_topo            = True
 res_topo@cnFillOn   = False
 res_topo@cnLinesOn  = True
 res_topo@gsnFrame   = False
 res_topo@gsnDraw    = False
 res_topo@cnLineLabelsOn         = False
 res_topo@cnLineThicknessF       = 4
 res_topo@cnLineDashPattern      = 16
 res_topo@cnLineColor            = "gray25" ; 
 res_topo@cnInfoLabelOn          = False
 res_topo@cnLevelSelectionMode   ="ExplicitLevels"
 res_topo@cnLevels               = (/1500/)
 res_topo@gsnLeftString = ""
 res_topo@gsnRightString = ""
;----------------------------------------------------

res = True
res@gsnDraw  = False
res@gsnFrame = False
res@gsnAddCyclic = False

res@mpOutlineOn = True
res@mpGeophysicalLineThicknessF = 2
res@mpGeophysicalLineColor = "grey25" ;"red" 
res@mpFillDrawOrder = "PreDraw"
res@mpFillOn = True

res@mpMinLatF = latstr
res@mpMaxLatF = latend
res@mpMinLonF = lonstr
res@mpMaxLonF = lonend

res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@lbLabelBarOn = True
res@cnInfoLabelOn = False

res@gsnStringFontHeightF = 0.035
res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025

res@cnFillPalette = "CBR_coldhot"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -0.1
res@cnMaxLevelValF = 0.1
res@cnLevelSpacingF = 0.02
res@pmLabelBarOrthogonalPosF = 0.12

;*******************************************
; time series (principal component) plot
;*******************************************
  ; t_eof_ts@long_name = "Amplitude"

  rts           = True
  rts@gsnDraw   = False       ; don't draw yet
  rts@gsnFrame  = False       ; don't advance frame yet
 ;rts@gsnScale  = True        ; force text scaling               

; these four resources allow the user to stretch the plot size, and
; decide exactly where on the page to draw it.

  rts@vpHeightF = 0.4        ; Changes the aspect ratio
  rts@vpWidthF  = 0.8
  ; rts@vpXF      = 0.10        ; change start locations
  ; rts@vpYF      = 0.75        ; the plot

  rts@gsnYRefLine           = 0.              ; reference line   
  rts@gsnXYBarChart         = True
  rts@gsnAboveYRefLineColor = "red"           ; above ref line fill red
  rts@gsnBelowYRefLineColor = "blue"          ; below ref line fill blue

  rts@gsnStringFontHeightF = 0.035
  rts@tmXBLabelFontHeightF = 0.025
  rts@tmYLLabelFontHeightF = 0.025
  rts@tiYAxisString = ""

  rts@trXMinF = yrStrt-1
  rts@trXMaxF = yrLast+1

  ; rts@gsnYRefLine = (/-std_val(0),0,std_val(0)/)
  ; rts@gsnYRefLineColor = "grey40"
  ; rts@gsnYRefLineThicknessF = 3
  ; rts@gsnYRefLineDashPattern = 2



    res@gsnLeftString = "EOF1"
    res@gsnCenterString = season
    res@gsnRightString = sprintf("%4.2f", pcvar)+"%"
    plot(0) = gsn_csm_contour_map(wks, t_eof(0,:,:), res)
    ; plot(i) = gsn_csm_contour_map(wks, t_clm(j,:,:), res)
    plot_topo(0) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(0), plot_topo(0))

    rts@gsnLeftString = "PC1"
    rts@gsnCenterString = season
    rts@gsnRightString = sprintf("%4.2f", pcvar)+"%"
    plot(1) = gsn_csm_xy(wks, year, t_eof_ts(0,:), rts)


pres = True
pres@gsnMaximize = True
; pres@gsnPanelLabelBar = True
pres@pmLabelBarOrthogonalPosF = -0.02

gsn_panel(wks, plot, (/1,2/), pres)

end