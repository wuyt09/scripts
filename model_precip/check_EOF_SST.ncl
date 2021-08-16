begin

lonstr = 120
lonend = 330
latstr = -20
latend = 20

yrStrt = 1959
yrLast = 2000
year = ispan(yrStrt, yrLast, 1)


diri = "/home/yangsong3/data-observation/SST/NOAA_ERSST/"
fili = "sst.mnmean.v5.nc"

f = addfile(diri+fili, "r")
TIME = f->time
YYYY = cd_calendar(TIME, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
time = f->time(iYYYY)

sst = f->sst(iYYYY,{latstr:latend},{lonstr:lonend})
sst = dtrend_msg_n(time, sst, False, False, 0)

neof   = 1        ; number of EOFs
optEOF = True       
optEOF@jopt = 0   ; This is the default; most commonly used; no need to specify.
;;optEOF@jopt = 1   ; **only** if the correlation EOF is desired
optETS = False

eof = eofunc_Wrap(sst(lat|:,lon|:,time|:), neof, optEOF)
pcvar = eof@pcvar
eof_ts = eofunc_ts_Wrap(sst(lat|:,lon|:,time|:), eof, optETS)
printVarSummary(eof)
printVarSummary(eof_ts)

;===================================================================
; Plot
;===================================================================
filipath = "/home/yangsong3/wuyt/sysu/figures/model_precip/"
filiname = "check_EOF_SST"
wks = gsn_open_wks("pdf", filipath+filiname)
plot = new(2, graphic)
plot_topo = plot

res = True
res@gsnDraw  = False
res@gsnFrame = False
res@gsnMaximize = True
res@gsnAddCyclic = False

res@mpMinLatF = latstr
res@mpMaxLatF = latend
res@mpMinLonF = lonstr
res@mpMaxLonF = lonend

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
res@cnLevelSpacingF = 0.01
res@cnMaxLevelValF = 0.1
res@cnMinLevelValF = -0.1
; res@cnLevelSpacingF = 0.2
; res@cnMaxLevelValF = 1
; res@cnMinLevelValF = -1
res@cnFillPalette = "BlueWhiteOrangeRed"

res@gsnRightString = "EOF1"

; res@gsnLeftString = "APHRO ("+yrStrt+"-"+yrLast+")"
plot(0) = gsn_csm_contour_map(wks, eof(0,:,:), res)
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

  ; rts@trXMinF = yrStrt-1
  ; rts@trXMaxF = yrLast+1

  rts@gsnRightString = "PC1"

  ; rts@gsnLeftString = "APHRO ("+yrStrt+"-"+yrLast+")"
  rts@gsnCenterString = sprintf("%4.2f", pcvar)+"%"
  plot(1) = gsn_csm_xy(wks, time, eof_ts(0,:), rts)
  ; plot(1) = gsn_csm_y(wks, eof_ts(0,:), rts)

pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
pres@lbLabelFontHeightF = 0.012
gsn_panel(wks, plot, (/2,1/), pres)


end