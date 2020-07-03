
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

prec_tlat = dim_avg_n_Wrap(pentad_prec(:,:,{10:20},:), 2)
prec_tlat&time = ispan(1, 73, 1)
printVarSummary(prec_tlat)

; x = new((/3,73,144/), float)
; x(0:1,:,:) = prec_tlat
; x(2,:,:) = prec_tlat(1,:,:) - prec_tlat(0,:,:)
; copy_VarCoords(prec_tlat(0,:,:), x(0,:,:))

;-------------------------------------------------------------------
Month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
title = (/"CTRL","NOHEAT","NOTP"/)

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20200513/"
figname = "precip_evolution_10-20N"

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
res@cnLineLabelsOn = False

res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = 5
res@cnMaxLevelValF = 20
res@cnLevelSpacingF  = 2                   ; contour spacing

res@trYMinF = 1
res@trYMaxF = 73
res@tmYLMode = "Explicit"
res@tmYLLabels = yaxis
res@tmYLValues = yaxis
res@tmYLLabelStride = 3
res@tmYLMinorOn = True
; res@tmYLMode = "Manual"
; res@tmYLTickEndF = 73
; res@tmYLTickStartF = 1
; res@tmYLTickSpacingF = 6
; ; res@tmYLMinorValues = 2
; res@tmYLMinLabelSpacingF = 2

res@tmXBMode = "Explicit"
res@tmXBValues = (/40,60,80,100,120,140,160,180/)
res@tmXBLabels = (/"40E","60E","80E","100E","120E","140E","160E","180E"/)
res@tmXBMinorOn = False

res@tmYLPrecision = 2
res@tiYAxisString    = "Pentad"              ; y axis title
res@gsnStringFontHeightF = 0.03

; res@lbLabelStride    = 2                   ; every other label
; res@cnLinesOn        = False               ; turn off contour lines
do i = 0, ncase-1
    res@gsnLeftString = title(i)
    plot(i) = gsn_csm_hov(wks, prec_tlat(i,:,{40:180}), res)
end do

; res@cnFillPalette = "CBR_drywet"
; res@cnLevelSelectionMode = "ManualLevels"
; res@cnMinLevelValF = -5
; res@cnMaxLevelValF = 5
; res@cnLevelSpacingF  = 1  
; res@gsnLeftString = title(2)
; plot(2) = gsn_csm_hov(wks, x(2,:,{40:180}), res)


pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True

gsn_panel(wks, plot, (/1,3/), pres)

end