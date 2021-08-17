begin

lonstr = 68
lonend = 107
latstr = 8
latend = 42

yrStrt = 1979
yrLast = 2017
year = ispan(yrStrt, yrLast, 1)
nyear = dimsizes(year)

diri = "/home/yangsong3/data-observation/GPCP/"
fili = "precip.mon.mean.nc"
f = addfile(diri+fili, "r")

TIME = f->time
YYYY = cd_calendar(TIME, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
precip = f->precip(iYYYY,{latstr:latend},{lonstr:lonend})

prec_med = (precip(6::12,:,:) + precip(7::12,:,:))/2
copy_VarCoords(precip(6::12,:,:), prec_med)
prec = dim_standardize_n_Wrap(prec_med, 1, 0)
printVarSummary(prec)



neof   = 1        ; number of EOFs
optEOF = True       
optEOF@jopt = 1   ; This is the default; most commonly used; no need to specify.
;;optEOF@jopt = 1   ; **only** if the correlation EOF is desired
optETS = False
optETS@jopt = 1

eof = eofunc_Wrap(prec(lat|:,lon|:,time|:), neof, optEOF)
pcvar = eof@pcvar
eof_ts = eofunc_ts_Wrap(prec(lat|:,lon|:,time|:), eof, optETS)
printVarSummary(eof)
printVarSummary(eof_ts)


;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)

 res_topo            = True
 res_topo@cnFillOn   = False
 res_topo@cnLinesOn  = True
 res_topo@gsnFrame   = False
 res_topo@gsnDraw    = False
 res_topo@cnLineLabelsOn         = False
 res_topo@cnLineThicknessF       = 3
 res_topo@cnLineColor            = (/"grey40"/)
 res_topo@cnInfoLabelOn          = False
 res_topo@cnLevelSelectionMode   ="ExplicitLevels"
 res_topo@cnLevels               = (/2000/)
 res_topo@gsnLeftString = ""
 res_topo@gsnRightString = ""
;===================================================================
; Plot
;===================================================================
filipath = "/home/yangsong3/wuyt/sysu/figures/model_precip/"
filiname = "EOF_JA_precip_GPCP_"+yrStrt+"_"+yrLast
wks = gsn_open_wks("pdf", filipath+filiname)
plot = new(2, graphic)
plot_topo = plot

res = True
res@gsnDraw  = False
res@gsnFrame = False
res@gsnMaximize = True
res@gsnAddCyclic = False

res@mpMinLatF = 10 ;latstr
res@mpMaxLatF = 40 ;latend
res@mpMinLonF = 70 ;lonstr
res@mpMaxLonF = 105 ;lonend

res@mpOutlineOn = True
res@mpFillOn = True
res@mpGeophysicalLineThicknessF = 2
res@mpGeophysicalLineColor = "grey25" ;"red"
; res@mpFillDrawOrder = "PostDraw"
res@mpCenterLonF = 180

res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False

res@gsnStringFontHeightF = 0.025
res@tmXBLabelFontHeightF = 0.018
res@tmYLLabelFontHeightF = 0.018

res@lbLabelBarOn = False
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 0.02
res@cnMaxLevelValF = 0.2
res@cnMinLevelValF = -0.2
; res@cnLevelSpacingF = 1.
; res@cnMaxLevelValF = 10
; res@cnMinLevelValF = -10
res@cnFillPalette = "MPL_BrBG"

res@gsnRightString = "EOF1"

res@gsnLeftString = "GPCP ("+yrStrt+"-"+yrLast+")"
plot(0) = gsn_csm_contour_map(wks, eof(0,:,:), res)
plot_topo(0) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(0), plot_topo(0))



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

  rts@gsnStringFontHeightF = 0.03
  rts@tmXBLabelFontHeightF = 0.025
  rts@tmYLLabelFontHeightF = 0.025

  rts@trXMinF = yrStrt-1
  rts@trXMaxF = yrLast+1

  rts@gsnRightString = "PC1"

  rts@gsnLeftString = "GPCP ("+yrStrt+"-"+yrLast+")"
  rts@gsnCenterString = sprintf("%4.2f", pcvar)+"%"
  plot(1) = gsn_csm_xy(wks, year, eof_ts(0,:), rts)

pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
pres@lbLabelFontHeightF = 0.012
gsn_panel(wks, plot, (/2,1/), pres)

end