
begin

lonstr = 50
lonend = 150
latstr = 10 ;-10
latend = 50
lonW = 90
lonE = 120
casename = (/"CTRL","NOTPQ1","NOTPQ1_MAR","NOTPQ1_APR","NOTPQ1_MAY"/)
ncase = dimsizes(casename)
pentad_u = new((/ncase,73,96,144/), float)
U = new((/25,365,96,144/), float)
LEV = 200

yy = ispan(6, 30, 1)

do i = 0, ncase-1
  diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(i)+"_F2000_F19G16_CAM4/pro/daily/"
  
  do j = 0, dimsizes(yy)-1
    yr = sprinti("%0.2i", yy(j))
    cname = "F2000_"+casename(i)+".plevel.U."+yr+".daily.nc"
    fc = addfile(diri+cname, "r")
    U(j,:,:,:) = fc->U(:,{LEV},:,:)
  end do
  print("end of reading")
  u = dim_avg_n_Wrap(U, 0) ; 365,96,144
  printVarSummary(u)
  lat = tofloat(u&lat)
  dlat = (lat(2)-lat(1))*0.0174533
  dy = 6378388.*dlat
  dudy = center_finite_diff_n(u, dy, True, 0, 1)
  copy_VarCoords(u, dudy)
  ; printVarSummary(dudy)

  pentad_p = dudy(0:72,:,:)

  do p = 0,72  ; a year has 73 penta
    pentad_p(p,:,:) = dim_sum_n_Wrap(dudy(5*p:5*p+4,:,:), 0)/5.0
  end do
pentad_u(i,:,:,:) = pentad_p
end do
pentad_u = pentad_u*1.e+6
pentad_u!0 = "case"
pentad_u&case= casename
pentad_u!1 = "time"
pentad_u&time = ispan(1, 73, 1)
printVarSummary(pentad_u)

pentad_ux = pentad_u(case|:,lat|:,time|:,lon|:)
prec_tlat = dim_avg_n_Wrap(pentad_ux(:,:,:,{lonW:lonE}), 3)
printVarSummary(prec_tlat)


; x = new((/3,73,144/), float)
; x(0:1,:,:) = prec_tlat
; x(2,:,:) = prec_tlat(1,:,:) - prec_tlat(0,:,:)
; copy_VarCoords(prec_tlat(0,:,:), x(0,:,:))

;-------------------------------------------------------------------
Month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
title = casename ;(/"CTRL","NOHEAT"/) ;,"NOTP"

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20200704/"
figname = "U200_jet_stream_axis_evolution_110-120E"

yaxis = ispan(1, 73, 2)
pens = 13
pene = 66 ;48

wks = gsn_open_wks("pdf", figpath+figname)

nplot = ncase
plot = new(nplot, graphic)
plots = plot


res          = True                ; plot mods desired
res@gsnDraw  = False
res@gsnFrame = False
res@vpHeightF = 0.4
res@vpWidthF  = 0.7

res@cnFillOn         = True                ; color on
res@cnFillPalette    = "BlueWhiteOrangeRed"    ; set color map
res@lbLabelBarOn = False
res@cnInfoLabelOn = False
res@cnLinesOn  = False
res@cnLineLabelsOn = False

res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -40
res@cnMaxLevelValF = 40
res@cnLevelSpacingF  = 5                  ; contour spacing


; res@trYMinF = 1
; res@trYMaxF = 73
; res@tmYLMode = "Explicit"
; res@tmYLLabels = yaxis
; res@tmYLValues = yaxis
; res@tmYLLabelStride = 3
; res@tmYLMinorOn = True

res@tmYLMinorOn = False
res@tmYLMode = "Explicit"
res@tmYLValues = ispan(latstr, latend, 5)
res@tmYLLabels = res@tmYLValues+"N" ;(/"10N","15N","20N","25N","30N","35N","40N","45N","50N"/)
res@tmYLLabelFontHeightF = 0.02

res@tmXBPrecision = 2
res@tiXAxisString    = "Pentad"              ; y axis title
res@gsnStringFontHeightF = 0.03

res@gsnRightString = "10^-6 S~S~-1~N~"
; res@lbLabelStride    = 2                   ; every other label
; res@cnLinesOn        = False               ; turn off contour lines
do i = 0, ncase-1
    res@gsnLeftString = title(i)
    plot(i) = gsn_csm_contour(wks, prec_tlat(i,{latstr:latend},pens-1:pene-1), res)
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
    plots(i) = gsn_csm_contour(wks, prec_tlat(i,{latstr:latend},pens-1:pene-1), res)
    overlay(plot(i), plots(i))
  end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = "dU/dy ("+lonW+"-"+lonE+"E)"

gsn_panel(wks, plot, (/2,3/), pres)

end