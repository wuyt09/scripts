
begin

lonstr = 60
lonend = 120
latstr = 20 ;-10
latend = 50

yrStrt = 1979
yrLast = 2015
year = ispan(yrStrt, yrLast, 1)
nyear = dimsizes(year)
;season = (/"MAM","SON"/)
season = (/"MAM","JJA","SON"/)
Month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
nmonth = dimsizes(Month)

;================================================================
filitp = "~/wuyt/data/topo.sfc.1x1.nc"
ftp = addfile(filitp, "r")
topo = ftp->topo(::-1,:)
printVarSummary(topo)
;================================================================

diri = "/home/yangsong3/data-observation/ERA-interim/monthly/"
fili = "air.mon.mean.nc"

f = addfile(diri+fili, "r")
TIME = f->time
YYYY = cd_calendar(TIME, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

T = short2flt(f->t(iYYYY,:,::-1,:))
printVarSummary(T)
printMinMax(T, True)

filip = "Surface_pressure.mon.mean.nc"
fp = addfile(diri+filip, "r")
PS = short2flt(fp->sp(iYYYY,::-1,:))

; print(T&latitude)
nlat = dimsizes(T&latitude)
nlon = dimsizes(T&longitude)
;===========integration========================================
levs = 600
levt = 300

lev = T&level
plev = lev*100.
psfc = PS
ptop = min(plev)

dP = dpres_plevel_Wrap(plev, psfc, ptop, 0)
dP!1 = "lev"
dP&lev = lev
printVarSummary(dP)

Tdp = T*dP
copy_VarCoords(T, Tdp)
TTT = dim_sum_n_Wrap(Tdp(:,{levs:levt},:,:), 1)/dim_sum_n_Wrap(dP(:,{levs:levt},:,:),1)
copy_VarCoords(T(:,0,:,:), TTT)
printVarSummary(TTT)
printMinMax(TTT, True)


; t = clmMonTLL(TTT)
t = reshape(TTT, (/nyear,nmonth,nlat,nlon/))
copy_VarCoords(T(0,0,:,:), t(0,0,:,:))
t!0 = "year"
t&year = year
t!1 = "month"
t&month = Month
printVarSummary(t)

neval = 1
optEOF = True
t_eof = new((/neval,nmonth,nlat,nlon/), float)
t_eof_ts = new((/neval,nmonth,nyear/), float)
pcvar = new((/nmonth/),float)
do i = 0, nmonth-1
    x = eofunc_n_Wrap(t(:,i,:,:), neval, optEOF, 0)
    t_eof(:,i,:,:) = x
    pcvar(i) = x@pcvar
    t_eof_ts(:,i,:) = eofunc_ts_n_Wrap(t(:,i,:,:), x, False, 0)
end do

printVarSummary(t_eof)
printVarSummary(t_eof_ts)

;=================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20200704/"
figname = "TTT_Mar_Apr_May_EOF_ERAI"

wks = gsn_open_wks("pdf", figpath+figname)
nplot = 6
plot = new(nplot, graphic)
plot_topo = new(3, graphic)

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
res@cnMinLevelValF = -0.01
res@cnMaxLevelValF = 0.01
res@cnLevelSpacingF = 0.002

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

  rts@trXMinF = yrStrt-1
  rts@trXMaxF = yrLast+1


do i = 0, 2
    j = i+2
    res@gsnLeftString = Month(j)
    res@gsnCenterString = "EOF1"
    res@gsnRightString = sprintf("%4.2f", pcvar(j))+"%"
    plot(i) = gsn_csm_contour_map(wks, t_eof(0,j,:,:), res)
    plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(i), plot_topo(i))

    rts@gsnLeftString = Month(j)
    rts@gsnCenterString = "PC1"
    rts@gsnRightString = sprintf("%4.2f", pcvar(j))+"%"
    plot(i+3) = gsn_csm_xy(wks, year, t_eof_ts(0,j,:), rts)
end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
pres@pmLabelBarOrthogonalPosF = -0.02

gsn_panel(wks, plot, (/2,3/), pres)




; ; panel plot only resources
;   rtsP                     = True             ; modify the panel plot
;   rtsP@gsnMaximize         = True             ; large format
      
;   resP@gsnPanelMainString  = "SCRIPPS: Pierce"
;   do n=0,neof-1
;      rts@gsnLeftString  = "EOF "+(n+1)
;      rts@gsnRightString = sprintf("%5.1f", eof@pcvar(n)) +"%"
;      plot(n) = gsn_csm_xy (wks,x&time,eof_ts(n,:),rts)
;   end do
;   gsn_panel(wks,plot,(/neof,1/),rtsP)        ; draw all 'neof' as one plot

end