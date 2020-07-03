
begin

casename = (/"CTRL","NOQ1","NOTP"/)
ncase = dimsizes(casename)
pentad_prec = new((/ncase,73,96,144/), float)

do i = 0, ncase-1
  diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(i)+"_F2000_F19G16_CAM4/pro/daily/"
  cname = "F2000_"+casename(i)+".PRECC.daily.06-30.nc"
  lname = "F2000_"+casename(i)+".PRECL.daily.06-30.nc"

  fc = addfile(diri+cname, "r")
  fl = addfile(diri+lname, "r")

  cp = fc->PRECC
  lp = fl->PRECL

  PREC = cp
  PREC = cp + lp
  PREC = PREC*1000.*86400
  PREC@units = "mm/day"

  prec = dim_avg_n_Wrap(PREC, 0) ; 365,96,144
  printVarSummary(prec)
  pentad_p = prec(0:72,:,:)

  do p = 0,72  ; a year has 73 penta
    pentad_p(p,:,:) = dim_sum_n_Wrap(prec(5*p:5*p+4,:,:), 0)/5.0
  end do
pentad_prec(i,:,:,:) = pentad_p
delete([/diri,cname,lname,fc,fl,cp,lp,PREC,prec,pentad_p/])
end do
printVarSummary(pentad_prec)

prec_tlat = dim_avg_n_Wrap(pentad_prec(:,:,:,{110:120}), 3)
prec_tlat&time = ispan(1, 73, 1)
printVarSummary(prec_tlat)
prec_tlat!0 = "case"
prec_tlat&case = casename
prec_tlatx = prec_tlat(case|:,lat|:,time|:)


;-------------------------------------------------------------------
Month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20200513/"
figname = "precip_evolution_110-120E"

title = (/"CTRL","NOHEAT","NOTP"/)

yaxis = ispan(1, 73, 2)

wks = gsn_open_wks("pdf", figpath+figname)

nplot = ncase
plot = new(nplot, graphic)


res          = True                ; plot mods desired
res@gsnDraw  = False
res@gsnFrame = False

res@cnFillOn         = True                ; color on
res@cnFillPalette    = "CBR_wet"    ; set color map
res@lbLabelBarOn = False
res@cnInfoLabelOn = False

res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = 5
res@cnMaxLevelValF = 20
res@cnLevelSpacingF  = 2                   ; contour spacing

res@tmXBMode = "Manual"
res@tmXBTickEndF = 73
res@tmXBTickStartF = 1
res@tmXBTickSpacingF = 6
res@tmXBPrecision = 2
; res@tmXBMode = "Explicit"
; res@tmXBLabels = yaxis
; res@tmXBValues = yaxis
; res@tmXBLabelStride = 3
; res@tmXBMinorOn = True
; res@tmXBPrecision = 2
; res@tmYLMode = "Manual"
; res@tmYLTickEndF = 73
; res@tmYLTickStartF = 1
; res@tmYLTickSpacingF = 6
; ; res@tmYLMinorValues = 2
; res@tmYLMinLabelSpacingF = 2

; res@tmYLMode = "Explicit"
; res@tmYLValues = (/10,20,30,40,50/)
; res@tmYLLabels = (/"10N","20N","30N","40N","50N"/)
; res@tmYLMinorOn = False

; res@tiXAxisString    = "Pentad"              ; y axis title

res@gsnStringFontHeightF = 0.03
; res@lbLabelStride    = 2                   ; every other label
; res@cnLinesOn        = False               ; turn off contour lines
do i = 0, ncase-1
    res@gsnLeftString = title(i)
    ; plot(i) = gsn_csm_time_lat(wks, prec_tlat(i,:,{10:50}), res)
    plot(i) = gsn_csm_lat_time(wks, prec_tlatx(i,{10:50},:), res)
end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True

gsn_panel(wks, plot, (/1,3/), pres)

end