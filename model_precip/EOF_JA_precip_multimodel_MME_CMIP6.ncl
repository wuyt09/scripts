begin

lonstr = 68
lonend = 107
latstr = 7
latend = 43

yrStrt = 1951
yrLast = 2014
year = ispan(yrStrt, yrLast, 1)
nyear = dimsizes(year)

;=======================================================
diriobs = "/home/yangsong3/data-observation/GPCP/"
filiobs = "precip.mon.mean.nc"
fobs = addfile(diriobs+filiobs, "r")
lat = fobs->lat({latstr:latend})
lon = fobs->lon({lonstr:lonend})
nlat = dimsizes(lat)
nlon = dimsizes(lon)
;=======================================================================
dirigcm = "/home/yangsong3/wuyt/sysu/scripts/model_precip/CMIP6/"

; models := (/"CAMS-CSM1-0","CanESM5","CESM2","CESM2-WACCM","CNRM-CM6-1","CNRM-ESM2-1",\
;            "E3SM-1-0","EC-Earth3","EC-Earth3-Veg","FGOALS-f3-L","FGOALS-g3","GISS-E2-1-G-CC","GISS-E2-1-G",\
;            "HadGEM3-GC31-LL","IPSL-CM6A-LR","MPI-ESM1-2-HR","NESM3","SAM0-UNICON","UKESM1-0-LL"/)
models = (/"BCC-CSM2-MR","BCC-ESM1","CAMS-CSM1-0","CanESM5","CESM2-FV2","CESM2","CESM2-WACCM-FV2",\
          "CESM2-WACCM","CIESM","E3SM-1-0","E3SM-1-1-ECA","E3SM-1-1","EC-Earth3","EC-Earth3-Veg",\
          "FGOALS-f3-L","FGOALS-g3","GISS-E2-1-G-CC","GISS-E2-1-G","GISS-E2-1-H","IPSL-CM6A-LR",\
          "KACE-1-0-G","MCM-UA-1-0","MIROC6","MPI-ESM1-2-HR","MRI-ESM2-0","NESM3","SAM0-UNICON"/)

nm := dimsizes(models)
pcvar = new(nm, float)
prec = new((/nm,nyear,nlat,nlon/), float)
prec!0 = "model"
prec&model = models
prec!1 = "year"
prec&year = year
prec!2 = "lat"
prec&lat = lat
prec!3 = "lon"
prec&lon = lon

eof = new((/nm,nlat,nlon/), float)
eof!0 = "model"
eof&model = models
eof!1 = "lat"
eof&lat = lat
eof!2 = "lon"
eof&lon = lon

eof_ts = new((/nm,nyear/), float)
eof_ts!0 = "model"
eof_ts&model = models
eof_ts!1 = "year"
eof_ts&year = year

neof   = 1        ; number of EOFs
optEOF = True       
optEOF@jopt = 0   ; This is the default; most commonly used; no need to specify.
;;optEOF@jopt = 1   ; **only** if the correlation EOF is desired
optETS = False

do i = 0, nm-1
    filigcm = systemfunc("ls "+dirigcm+"pr_Amon_"+models(i)+"_historical_r1i1p1f*_*.nc")
    fgcm = addfile(filigcm, "r")
    TIME := fgcm->time
    YYYY := cd_calendar(TIME, -1)/100
    iYYYY := ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
    PR_GCMs = fgcm->pr(iYYYY,{latstr:latend},{lonstr:lonend})
    PR_GCMs!1 = "lat"
    PR_GCMs!2 = "lon"
    PR_GCMs = PR_GCMs * 86400
    PR_GCMx = linint2_Wrap(PR_GCMs&lon, PR_GCMs&lat, PR_GCMs, True, lon, lat, 0)
    pr_jul = PR_GCMx(6::12,:,:)
    pr_aug = PR_GCMx(7::12,:,:)
    prec_med1 = (pr_aug + pr_jul)/2
    prec_med = dim_standardize_n_Wrap(prec_med1, 1, 0)
    copy_VarCoords(pr_aug, prec_med)

    prec(i,:,:,:) = (/prec_med/)
    delete([/PR_GCMs,PR_GCMx,pr_jul,pr_aug,prec_med1/]) ;

    eof_med = eofunc(prec_med(lat|:,lon|:,time|:), neof, optEOF)
    pcvar(i) = eof_med@pcvar
    eof_ts_med = eofunc_ts(prec_med(lat|:,lon|:,time|:), eof_med, optETS)
    eof(i,:,:) = (/eof_med(0,:,:)/)
    eof_ts(i,:) = (/eof_ts_med(0,:)/)
    delete([/prec_med,eof_med,eof_ts_med/])
    print("end of "+models(i))
end do
printVarSummary(prec)
printVarSummary(eof_ts)

prec_mme = dim_avg_n_Wrap(prec, 0)

eof_mme_med = eofunc(prec_mme(lat|:,lon|:,year|:), neof, optEOF)
eof_mme = eof_mme_med(0,:,:)
eof_ts_mme_med = eofunc_ts(prec_mme(lat|:,lon|:,year|:), eof_mme_med, optETS)
eof_ts_mme = eof_ts_mme_med(0,:)
copy_VarCoords(prec_mme(0,:,:), eof_mme)
copy_VarCoords(prec_mme(:,0,0), eof_ts_mme)
printVarSummary(eof_mme)


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
 res_topo@cnLineColor            = (/"grey40","black"/)
 res_topo@cnInfoLabelOn          = False
 res_topo@cnLevelSelectionMode   ="ExplicitLevels"
 res_topo@cnLevels               = (/2000,3000/)
 res_topo@gsnLeftString = ""
 res_topo@gsnRightString = ""
;===================================================================
; Plot
;===================================================================
filipath = "/home/yangsong3/wuyt/sysu/figures/model_precip/"
filiname = "EOF_JA_precip_multimodel_MME_CMIP6_"+yrStrt+"_"+yrLast
wks = gsn_open_wks("pdf", filipath+filiname)
plot = new(nm+1, graphic)
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

res@gsnStringFontHeightF = 0.03
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
do i = 0, nm-1
    res@gsnLeftString = models(i)
    plot(i) = gsn_csm_contour_map(wks, eof(i,:,:), res)
    plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(i), plot_topo(i))
end do

res@gsnLeftString = "MME"
plot(nm) = gsn_csm_contour_map(wks, eof_mme, res)
plot_topo(nm) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(nm), plot_topo(nm))


pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
pres@lbLabelFontHeightF = 0.012
gsn_panel(wks, plot, (/4,7/), pres)


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

  rts@gsnRightString = "PC1"

  do i = 0, nm-1
    rts@gsnLeftString = models(i)
    rts@gsnCenterString = sprintf("%4.2f", pcvar(i))+"%"
    plot(i) = gsn_csm_xy(wks, year, eof_ts(i,:), rts)
  end do
  
  rts@gsnLeftString = "MME"
  plot(nm) = gsn_csm_xy(wks, year, eof_ts_mme, rts)
  gsn_panel(wks, plot, (/4,7/), pres)

end