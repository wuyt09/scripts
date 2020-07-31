
begin

lonstr = 50
lonend = 150
latstr = 10 ;-10
latend = 50
lonW = 60
lonE = 120
casename = (/"CTRL","NOTPQ1","NOTPQ1_MAR","NOTPQ1_APR","NOTPQ1_MAY"/)
ncase = dimsizes(casename)
dudy = new((/ncase,12,96,144/), float)
LEV = 200

yy = ispan(6, 30, 1)

do i = 0, ncase-1
  diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(i)+"_F2000_F19G16_CAM4/pro/"
  uname = "F2000_"+casename(i)+".plevel.U.mon.0601-3012.nc"
  fu = addfile(diri+uname, "r")
  U = fu->U(:,{LEV},:,:)
  u = clmMonTLL(U)
  printVarSummary(u)
  lat = tofloat(u&lat)
  dlat = (lat(2)-lat(1))*0.0174533
  dy = 6378388.*dlat
  dudy(i,:,:,:) = center_finite_diff_n(u, dy, True, 0, 1)
  copy_VarCoords(u, dudy(0,:,:,:))
  ; printVarSummary(dudy)
end do

dudy!0 = "case"
dudy&case= casename
dudy!1 = "month"
dudy&month = ispan(1, 12, 1)
printVarSummary(dudy)

dudy_x = dudy(case|:,lat|:,month|:,lon|:)
dudy_tlat = dim_avg_n_Wrap(dudy_x(:,:,:,{lonW:lonE}), 3)
printVarSummary(dudy_tlat)


; x = new((/3,73,144/), float)
; x(0:1,:,:) = prec_tlat
; x(2,:,:) = prec_tlat(1,:,:) - prec_tlat(0,:,:)
; copy_VarCoords(prec_tlat(0,:,:), x(0,:,:))

;-------------------------------------------------------------------
Month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
title = casename ;(/"CTRL","NOHEAT"/) ;,"NOTP"

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20200704/"
figname = "U200_jet_stream_axis_evolution_110-120E_monthly"


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

; res@cnLevelSelectionMode = "ManualLevels"
; res@cnMinLevelValF = -20
; res@cnMaxLevelValF = 20
; res@cnLevelSpacingF  = 4                  ; contour spacing


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
res@tiXAxisString    = "month"              ; y axis title
res@gsnStringFontHeightF = 0.03

; res@lbLabelStride    = 2                   ; every other label
; res@cnLinesOn        = False               ; turn off contour lines
do i = 0, ncase-1
    res@gsnLeftString = title(i)
    plot(i) = gsn_csm_contour(wks, dudy_tlat(i,{latstr:latend},:), res)
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
    plots(i) = gsn_csm_contour(wks, dudy_tlat(i,{latstr:latend},:), res)
    overlay(plot(i), plots(i))
  end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = "dU/dy ("+lonW+"-"+lonE+"E)"

gsn_panel(wks, plot, (/2,3/), pres)

end