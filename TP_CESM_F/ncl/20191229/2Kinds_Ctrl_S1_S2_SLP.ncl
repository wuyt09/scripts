
begin

lonstr = 0
lonend = 360
latstr = 0
latend = 60

casename = (/"CTRL","NOQ1","NOTP_NOQ1"/) ;,"NOTP"
ncase = dimsizes(casename)
month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
nm = dimsizes(month)
season = (/"MAM","JJA","SON","DJF"/)

ts = new((/ncase,12,96,144/), float)


do i = 0, dimsizes(casename)-1
    diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(i)+"_F2000_F19G16_CAM4/pro/"
    fili = "F2000_"+casename(i)+".PSL.mon.0601-3012.nc"

    f = addfile(diri+fili, "r")
    TS = f->PSL
    ts(i,:,:,:) = clmMonTLL(TS)
    ts = ts/100
end do
printVarSummary(ts)

px = ts(0:1,:,:,:)

px(0,:,:,:) = ts(0,:,:,:) - ts(1,:,:,:)
px(1,:,:,:) = ts(1,:,:,:) - ts(2,:,:,:)

printVarSummary(px)

;----------------------------------------------------------------
tpmaskpath = "/home/yangsong3/wuyt/sysu/scripts/TP_CESM_F/model_need_input_data/"
tpmaskname = "TPmask_1000m_C20190410.nc"
ftp = addfile(tpmaskpath+tpmaskname, "r")

tpmask = ftp->tp
printVarSummary(tpmask)
px = mask(px, tpmask.ge.0.9, False)
;==================================================================
;==================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20191229/"

figname = "2Kind_Ctrl_S1_S2_SLP"

wks = gsn_open_wks("pdf", figpath+figname)
nplot = 12*2
plot = new(nplot, graphic)
plot_topo = plot
plotv = plot
plot_prob = plot


   res_topo            = True
   res_topo@cnFillOn   = False
   res_topo@cnLinesOn  = True
   res_topo@gsnFrame   = False
   res_topo@gsnDraw    = False
   res_topo@cnLineLabelsOn         = False
   res_topo@cnLineThicknessF       = 3
   res_topo@cnLineColor            = "blue"
   res_topo@cnInfoLabelOn          = False
   res_topo@cnLevelSelectionMode   ="ExplicitLevels"
   res_topo@cnLevels               = (/0.9/)
   res_topo@gsnLeftString = ""
   res_topo@gsnRightString = ""


  res = True
  res@gsnDraw  = False
  res@gsnFrame = False

  res@mpOutlineOn = True
  res@mpFillOn = True
  res@mpGeophysicalLineThicknessF = 2
  res@mpGeophysicalLineColor = "grey45" ;"grey25"
  res@mpFillDrawOrder = "PreDraw"

  res@mpMinLatF = latstr
  res@mpMaxLatF = latend
  res@mpMinLonF = lonstr
  res@mpMaxLonF = lonend
  res@mpCenterLonF = (lonstr+lonend)/2

  res@cnFillOn = True
  res@cnLinesOn = False
  res@cnLineLabelsOn = False
  res@cnInfoLabelOn = False

  res@gsnStringFontHeightF = 0.035
  res@tmXBLabelFontHeightF = 0.025
  res@tmYLLabelFontHeightF = 0.025
  res@pmLabelBarOrthogonalPosF = 0.15

  ;------------------------------------------
  
    ;------------------------
    ; draw prob
    ;------------------------
    res3  = True
    res3@gsnDraw   = False
    res3@gsnFrame  = False

    res3@cnFillOn  = False 
    res3@cnLinesOn = False
    res3@cnLineLabelsOn = False
    res3@cnInfoLabelOn  = False
    res3@cnFillColor = "blue"

    res3@gsnRightString = ""

    res3@cnLevelSelectionMode = "ManualLevels"
    res3@cnMinLevelValF = 0.00
    res3@cnMaxLevelValF = 0.10
    res3@cnLevelSpacingF = 0.01

    opt  = True
    opt@gsnShadeFillType = "pattern"
    opt@gsnShadeLow = 17
    opt@gsnShadeFillScaleF = 0.8
    opt@gsnShadeFillDotSizeF = 0.003
  ;------------------------------------------

  ;   res@lbLabelBarOn = False
  ;   res@cnLevelSelectionMode = "ManualLevels"

  ; res@cnLevelSpacingF = 0.5
  ; res@cnMaxLevelValF = 5.
  ; res@cnMinLevelValF = -5.
  ; res@cnFillPalette = "CBR_coldhot"

  res@gsnLeftString = "Thermal"
  do i = 0, dimsizes(month)-1
    res@gsnRightString = month(i)
    plot(2*i) = gsn_csm_contour_map(wks, px(0,i,:,:), res)
    plot_topo(2*i) = gsn_csm_contour(wks, tpmask, res_topo)
    overlay(plot(2*i), plot_topo(2*i))
  end do
  
  res@gsnLeftString = "Dynamical"
  do i = 0, dimsizes(month)-1
    res@gsnRightString = month(i)
    plot(2*i+1) = gsn_csm_contour_map(wks, px(1,i,:,:), res)
    plot_topo(2*i+1) = gsn_csm_contour(wks, tpmask, res_topo)
    overlay(plot(2*i+1), plot_topo(2*i+1))

    ; plot_prob(2*i+1) = gsn_csm_contour(wks, probp(1,i,:,:), res3)
    ; plot_prob(2*i+1) = gsn_contour_shade(plot_prob(2*i+1), siglvl, -999., opt)
    ; ;overlay(plot(2*i+1), plot_prob(2*i+1))
  end do
  
  pres = True
  pres@gsnMaximize = True
  pres@gsnPanelMainString = "SLP"
  pres@gsnPanelLabelBar = True
  ;pres@pmLabelBarOrthogonalPosF = -0.01
  pres@lbLabelFontHeightF = 0.01
  gsn_panel(wks, plot, (/12,2/), pres)


;==================================================================
;==================================================================
ppx = px
monthx = month
do i = 0, 8
    ppx(:,i,:,:) = px(:,i+3,:,:)
end do
do i = 9, 11
    ppx(:,i,:,:) = px(:,i-9,:,:)
end do


p_warm = ppx(:,0,:,:)
p_warm = dim_avg_n_Wrap(ppx(:,0:5,:,:), 1)
p_cold = ppx(:,0,:,:)
p_cold = dim_avg_n_Wrap(ppx(:,6:11,:,:), 1)

res@gsnLeftString = "Thermal"
res@gsnRightString = "Warm season"
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 1
res@cnMaxLevelValF = 10.
res@cnMinLevelValF = -10.
res@cnFillPalette = "CBR_coldhot"
plot(0) = gsn_csm_contour_map(wks, p_warm(0,:,:), res)
plot_topo(0) = gsn_csm_contour(wks, tpmask, res_topo)
overlay(plot(0), plot_topo(0))

res@gsnLeftString = "Dynamical"
plot(1) = gsn_csm_contour_map(wks, p_warm(1,:,:), res)
plot_topo(1) = gsn_csm_contour(wks, tpmask, res_topo)
overlay(plot(1), plot_topo(1))

res@gsnLeftString = "Thermal"
res@gsnRightString = "Cold season"
plot(2) = gsn_csm_contour_map(wks, p_cold(0,:,:), res)
plot_topo(2) = gsn_csm_contour(wks, tpmask, res_topo)
overlay(plot(2), plot_topo(2))

res@gsnLeftString = "Dynamical"
plot(3) = gsn_csm_contour_map(wks, p_cold(1,:,:), res)
plot_topo(3) = gsn_csm_contour(wks, tpmask, res_topo)
overlay(plot(3), plot_topo(3))

pres@gsnPanelMainString = "SLP"
pres@gsnPanelLabelBar = True
;pres@pmLabelBarOrthogonalPosF = -0.01
pres@lbLabelFontHeightF = 0.01
gsn_panel(wks, plot, (/2,2/), pres)


end