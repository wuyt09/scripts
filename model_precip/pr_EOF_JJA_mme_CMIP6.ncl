
begin

lonstr = 60
lonend = 110
latstr = 5
latend = 45

season = "JJA"

diriobs = "/home/yangsong3/data-observation/GPCP/"
filiobs = "precip.mon.mean.nc"
fobs = addfile(diriobs+filiobs, "r")

plat = fobs->lat
plon = fobs->lon



yrStrt = 1861
yrLast = 2014
year6 = ispan(yrStrt, yrLast, 1)
nyr6 = dimsizes(year6)

dirigcm = "/home/yangsong3/wuyt/sysu/scripts/model_precip/CMIP6/"

; models := (/"CAMS-CSM1-0","CanESM5","CESM2","CESM2-WACCM","CNRM-CM6-1","CNRM-ESM2-1",\
;            "E3SM-1-0","EC-Earth3","EC-Earth3-Veg","FGOALS-f3-L","FGOALS-g3","GISS-E2-1-G-CC","GISS-E2-1-G",\
;            "HadGEM3-GC31-LL","IPSL-CM6A-LR","MPI-ESM1-2-HR","NESM3","SAM0-UNICON","UKESM1-0-LL"/)
models = (/"CAMS-CSM1-0","CanESM5"/)

nm = dimsizes(models)
pr_gcm = new((/nm,nyr6,72,144/), float)
pr_gcm!0 = "model"
pr_gcm&model = models

do i = 0, nm-1
    filigcm = systemfunc("ls "+dirigcm+"pr_Amon_"+models(i)+"_historical_r1i1p1f*_g*_*.nc")
    fgcm = addfile(filigcm, "r")
    TIME := fgcm->time
    YYYY := cd_calendar(TIME, -1)/100
    iYYYY := ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
    PR_GCMs = fgcm->pr(iYYYY,:,:)
    PR_GCMs = PR_GCMs * 86400
    PR_GCMx = linint2_Wrap(PR_GCMs&lon, PR_GCMs&lat, PR_GCMs, True, plon, plat, 0)
    pr_gcm(i,:,:,:) = month_to_season(PR_GCMx, season)
    ; pr_gcm_ave(i,:,:) = dim_avg_n_Wrap(pr_gcm, 0)
    ; copy_VarCoords(pr_gcm(0,:,:), pr_gcm_ave(i,:,:))

    delete([/PR_GCMs,PR_GCMx/])
    print("end of "+models(i))
    print("============================================")
end do
printVarSummary(pr_gcm)

pr_gcm_mme = dim_avg_n_Wrap(pr_gcm, 0)
printVarSummary(pr_gcm_mme)

TP_IN = pr_gcm_mme(:,{latstr:latend},{lonstr:lonend})

neval = 1

x = eofcor_Wrap(TP_IN(lat|:,lon|:,time|:), neval)
printVarSummary(x)

optEOF = True
x_eof = eofunc_Wrap(TP_IN(lat|:,lon|:,time|:), neval, optEOF)
pcvar = x_eof@pcvar
x_ts = eofunc_ts_Wrap(TP_IN(lat|:,lon|:,time|:), x_eof, False)
printVarSummary(x_eof)
printVarSummary(x_ts)

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
 res_topo@cnLineColor            = "grey40"
 res_topo@cnInfoLabelOn          = False
 res_topo@cnLevelSelectionMode   ="ExplicitLevels"
 res_topo@cnLevels               = (/2000/)
 res_topo@gsnLeftString = ""
 res_topo@gsnRightString = ""

;===================================================================
; Plot
;===================================================================
filipath = "/home/yangsong3/wuyt/sysu/figures/model_precip/"
filiname = "pr_EOF_JJA_mme_CMIP6"
wks = gsn_open_wks("pdf", filipath+filiname)
plot = new(2, graphic)
plot_topo = plot
plotv = plot
; plot_prob = plot

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
; res@mpCenterLonF = 180

res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False

res@gsnStringFontHeightF = 0.02
res@tmXBLabelFontHeightF = 0.018
res@tmYLLabelFontHeightF = 0.018

; res@lbLabelBarOn = False
; res@pmLabelBarOrthogonalPosF = 0.2
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 0.02
res@cnMaxLevelValF = 0.2
res@cnMinLevelValF = -0.2
res@cnFillPalette = "MPL_BrBG"
; res@gsnRightString = "mm/day"

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

; gres = True
; gres@gsLineColor = "red"
; gres@gsLineThicknessF = 3.0
; glat1 = (/40,40,28,28,40/)
; glon1 = (/85,105,105,85,85/)
; glat2 = (/28,28,18,18,28/)
; glon2 = (/72,85,85,72,72/)

res@gsnLeftString = "EOF1"
plot(0) = gsn_csm_contour_map(wks, x_eof(0,:,:), res)
plot_topo(0) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(0), plot_topo(0))
; dum1 = gsn_add_polyline(wks,plot(0),glon1,glat1,gres)
; dum2 = gsn_add_polyline(wks,plot(0),glon2,glat2,gres)

rts@gsnLeftString = "PC1"
rts@gsnRightString = sprintf("%4.2f", pcvar)+"%"
plot(1) = gsn_csm_xy(wks, year6, x_ts(0,:), rts)

pres = True
pres@gsnMaximize = True
; pres@gsnPanelLabelBar = True
pres@lbLabelFontHeightF = 0.012
gsn_panel(wks, plot, (/2,1/), pres)

; draw(plot(0))


end