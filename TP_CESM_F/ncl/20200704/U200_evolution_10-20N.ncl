
begin

casename = (/"CTRL","NOTPQ1","NOTPQ1_MAR","NOTPQ1_APR","NOTPQ1_MAY"/)
ncase = dimsizes(casename)
pentad_prec = new((/ncase,73,96,144/), float)
Z3 = new((/25,365,96,144/), float)
LEV = 200

yy = ispan(6, 30, 1)

do i = 0, ncase-1
  diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(i)+"_F2000_F19G16_CAM4/pro/daily/"
  
  do j = 0, dimsizes(yy)-1
    yr = sprinti("%0.2i", yy(j))
    cname = "F2000_"+casename(i)+".plevel.U."+yr+".daily.nc"
    fc = addfile(diri+cname, "r")
    Z3(j,:,:,:) = fc->U(:,{LEV},:,:)
  end do
  print("end of reading")
  prec = dim_avg_n_Wrap(Z3, 0) ; 365,96,144
  printVarSummary(prec)
  pentad_p = prec(0:72,:,:)

  do p = 0,72  ; a year has 73 penta
    pentad_p(p,:,:) = dim_sum_n_Wrap(prec(5*p:5*p+4,:,:), 0)/5.0
  end do
pentad_prec(i,:,:,:) = pentad_p
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
title = casename ;(/"CTRL","NOHEAT"/) ;,"NOTP"

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20200704/"
figname = "U200_evolution_10-20N"

yaxis = ispan(1, 73, 2)

wks = gsn_open_wks("pdf", figpath+figname)

nplot = ncase
plot = new(nplot, graphic)
plots = plot


res          = True                ; plot mods desired
res@gsnDraw  = False
res@gsnFrame = False

res@cnFillOn         = True                ; color on
res@cnFillPalette    = "BlueWhiteOrangeRed"    ; set color map
res@lbLabelBarOn = False
res@cnInfoLabelOn = False
res@cnLinesOn  = False
res@cnLineLabelsOn = False

res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -20
res@cnMaxLevelValF = 20
res@cnLevelSpacingF  = 4                  ; contour spacing

; res@trYMinF = 1
; res@trYMaxF = 73
; res@tmYLMode = "Explicit"
; res@tmYLLabels = yaxis
; res@tmYLValues = yaxis
; res@tmYLLabelStride = 3
; res@tmYLMinorOn = True

res@tmXBMinorOn = False
res@tmXBMode = "Explicit"
res@tmXBValues =(/50,60,70,80,90,100,110,120,130,140/)  ;(/60,80,100,120,140/)
res@tmXBLabels = (/"50E","60E","70E","80E","90E","100E","110E","120E","130E","140E"/)  ;(/"60E","80E","100E","120E","140E"/)
res@tmXBLabelFontHeightF = 0.02

res@tmYLPrecision = 2
res@tiYAxisString    = "Pentad"              ; y axis title
res@gsnStringFontHeightF = 0.03

; res@lbLabelStride    = 2                   ; every other label
; res@cnLinesOn        = False               ; turn off contour lines
do i = 0, ncase-1
    res@gsnLeftString = title(i)
    plot(i) = gsn_csm_hov(wks, prec_tlat(i,13:65,{50:140}), res)
end do

; res@cnFillPalette = "BlueWhiteOrangeRed"
; res@cnLevelSelectionMode = "ManualLevels"
; res@cnMinLevelValF = -20
; res@cnMaxLevelValF = 20
; res@cnLevelSpacingF  = 4  
; res@gsnLeftString = title(2)
; plot(2) = gsn_csm_hov(wks, x(2,18:35,{90:160}), res)


  delete(res@gsnLeftString)
  delete(res@gsnRightString)
  res@cnLinesOn = True
  res@cnFillOn  = False
  res@cnLineLabelsOn = False
  res@cnLevelSelectionMode = "ExplicitLevels"
  res@cnLevels = (/0/)
  res@cnMonoLineThickness = True
  res@cnLineThicknessF = 5
  res@gsnLeftString = ""
  res@gsnRightString = ""
  ;res@gsnContourZeroLineThicknessF = 5
  do i = 0, ncase-1
    plots(i) = gsn_csm_hov(wks, prec_tlat(i,13:65,{50:140}), res)
    overlay(plot(i), plots(i))
  end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = "U200 (10-20N)"

gsn_panel(wks, plot, (/2,3/), pres)

end