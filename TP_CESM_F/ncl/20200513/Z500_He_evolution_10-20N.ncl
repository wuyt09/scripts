
begin

casename = (/"CTRL","NOQ1"/) ;,"NOTP"
ncase = dimsizes(casename)
pentad_prec = new((/ncase,73,96,144/), float)
Z3 = new((/25,365,96,144/), float)
LEV = 500

yy = ispan(6, 30, 1)

do i = 0, ncase-1
  diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(i)+"_F2000_F19G16_CAM4/pro/daily/"
  
  do j = 0, dimsizes(yy)-1
    yr = sprinti("%0.2i", yy(j))
    cname = "F2000_"+casename(i)+".plevel.Z3."+yr+".daily.nc"
    fc = addfile(diri+cname, "r")
    Z = fc->Z3(:,{LEV},:,:)
    Z_zonal = dim_avg_n_Wrap(Z(:,{0:20},:), 1)
    ZX = conform_dims(dimsizes(Z), Z_zonal, (/0,2/))
    Z3(j,:,:,:) = Z - ZX
    copy_VarCoords(Z, Z3(j,:,:,:))
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
title = (/"CTRL","NOHEAT"/) ;,"NOTP"

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20200513/"
figname = "Z500_He_evolution_10-20N"

yaxis = ispan(1, 73, 2)

wks = gsn_open_wks("pdf", figpath+figname)

nplot = ncase
plot = new(nplot, graphic)
plots = plot


res          = True                ; plot mods desired
res@gsnDraw  = False
res@gsnFrame = False

res@cnFillOn         = True                ; color on
res@cnFillPalette    = "WhiteYellowOrangeRed"    ; set color map
res@lbLabelBarOn = False
res@cnInfoLabelOn = False
res@cnLinesOn  = False
res@cnLineLabelsOn = False

; res@cnLevelSelectionMode = "ManualLevels"
; res@cnMinLevelValF = 5870
; res@cnMaxLevelValF = 5900
; res@cnLevelSpacingF  = 3                  ; contour spacing

; res@trYMinF = 1
; res@trYMaxF = 73
; res@tmYLMode = "Explicit"
; res@tmYLLabels = yaxis
; res@tmYLValues = yaxis
; res@tmYLLabelStride = 3
; res@tmYLMinorOn = True

res@tmXBMinorOn = False
res@tmXBMode = "Explicit"
res@tmXBValues =(/90,100,110,120,130,140,150,160/)  ;(/60,80,100,120,140/)
res@tmXBLabels = (/"90E","100E","110E","120E","130E","140E","150E","160E"/)  ;(/"60E","80E","100E","120E","140E"/)
res@tmXBLabelFontHeightF = 0.02

res@tmYLPrecision = 2
res@tiYAxisString    = "Pentad"              ; y axis title
res@gsnStringFontHeightF = 0.03

; res@lbLabelStride    = 2                   ; every other label
; res@cnLinesOn        = False               ; turn off contour lines
do i = 0, ncase-1
    res@gsnLeftString = title(i)
    ; plot(i) = gsn_csm_hov(wks, prec_tlat(i,18:35,{90:160}), res)
    ; plot(i) = gsn_csm_hov(wks, prec_tlat(i,19:54,{90:160}), res)
    plot(i) = gsn_csm_hov(wks, prec_tlat(i,18:35,{90:160}), res)
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
  res@cnMonoLineThickness = True
  res@cnLineThicknessF = 5
  res@gsnLeftString = ""
  res@gsnRightString = ""
  ;res@gsnContourZeroLineThicknessF = 5

  res@cnLevels = (/0/)
  plots(0) = gsn_csm_hov(wks, prec_tlat(0,18:35,{50:140}), res)
  overlay(plot(0), plots(0))
  ; res@cnLevels = (/5885/)
  plots(1) = gsn_csm_hov(wks, prec_tlat(1,18:35,{50:140}), res)
  overlay(plot(1), plots(1))



pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = "Z500 (10-20N): WPSH"

gsn_panel(wks, plot, (/1,3/), pres)

end