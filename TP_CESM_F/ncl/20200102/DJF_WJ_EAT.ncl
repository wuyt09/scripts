
begin

lonstr = 0
lonend = 200
latstr = 0
latend = 60

casename = (/"CTRL","NOQ1","NOTP_NOQ1"/) ;,"NOTP"
ncase = dimsizes(casename)
; month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
; nm = dimsizes(month)
season = (/"DJF"/)

LEV = (/200,500/)
z = new((/ncase,96,144/), float)
u = z

do i = 0, dimsizes(casename)-1
    diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(i)+"_F2000_F19G16_CAM4/pro/"

    filiz = "F2000_"+casename(i)+".plevel.Z3.mon.0601-3012.nc"
    fz = addfile(diri+filiz, "r")
    Z1 = fz->Z3(:,{LEV(1)},:,:)
    z1 = month_to_season(Z1, season)
    z(i,:,:) = dim_avg_n_Wrap(z1, 0)

    filiu = "F2000_"+casename(i)+".plevel.U.mon.0601-3012.nc"
    fu = addfile(diri+filiu, "r")
    U1 = fu->U(:,{LEV(0)},:,:)
    u1 = month_to_season(U1, season)
    u(i,:,:) = dim_avg_n_Wrap(u1, 0)

end do
; z = z/9.8
; z@units = "gpm"
z!0 = "case"
z&case = casename
printVarSummary(z)
u!0 = "case"
u&case = casename

;----------------------------------------------------------------
tpmaskpath = "/home/yangsong3/wuyt/sysu/scripts/TP_CESM_F/model_need_input_data/"
tpmaskname = "TPmask_1000m_C20190410.nc"
ftp = addfile(tpmaskpath+tpmaskname, "r")

tpmask = ftp->tp
printVarSummary(tpmask)
;----------------------------------------------------------------

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20200102/"

figname = "DJF_WJ_EAT"

wks = gsn_open_wks("pdf", figpath+figname)
nplot = ncase
plot = new(nplot, graphic)
plotc = plot
plot_topo = plot

;=====================================================
   res_topo            = True
   res_topo@cnFillOn   = False
   res_topo@cnLinesOn  = True
   res_topo@gsnFrame   = False
   res_topo@gsnDraw    = False
   res_topo@cnLineLabelsOn         = False
   res_topo@cnLineThicknessF       = 5
   res_topo@cnLineDashPattern      = 2
   res_topo@cnLineColor            = "gray40"
   res_topo@cnInfoLabelOn          = False
   res_topo@cnLevelSelectionMode   ="ExplicitLevels"
   res_topo@cnLevels               = (/0.9/)
   res_topo@gsnLeftString = ""
   res_topo@gsnRightString = ""
;=====================================================

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

    res@gsnStringFontHeightF = 0.03
    res@tmXBLabelFontHeightF = 0.02
    res@tmYLLabelFontHeightF = 0.02
    ;res@pmLabelBarOrthogonalPosF = 0.12
    res@lbLabelBarOn = False
    ;res@lbLabelFontHeightF = 0.02
    res@cnFillPalette = "WhiteYellowOrangeRed"

;=====================================================
    cnres = True
    cnres@gsnDraw = False
    cnres@gsnFrame = False

    cnres@cnFillOn = False
    cnres@cnLinesOn = True
    cnres@cnLineLabelsOn = True
    cnres@cnLineLabelFontColor = "red"
    cnres@cnLineColor = "green4"
    cnres@cnLineThicknessF = 5
    cnres@cnInfoLabelOn = False
    ;cnres@cnLineDrawOrder = "PostDraw"
    cnres@cnLabelDrawOrder = "PostDraw"
    ;cnres@cnLevelFlag = "LineAndLabel"
    cnres@cnLineLabelInterval = 1

;=====================================================

; res@cnLevelSelectionMode = "ExplicitLevels"
; res@cnLevels = (/12460,12470,12480,12500,12520,12540,12560/)
; ;res@cnFillColors = ()
res@cnLevelSelectionMode = "ManualLevels"
res@cnMaxLevelValF = 80
res@cnMinLevelValF = 40
res@cnLevelSpacingF = 10


cnres@cnLevelSelectionMode = "ManualLevels"
cnres@cnMaxLevelValF = 5500
cnres@cnMinLevelValF = 5100
cnres@cnLevelSpacingF = 100

leftsti = (/"Ctrl","NT","NTNH"/)
do i = 0, ncase-1
  res@gsnLeftString = leftsti(i)+"_U200_H500"
  plot(i) = gsn_csm_contour_map(wks, u(i,:,:), res)
  plotc(i) = gsn_csm_contour(wks, z(i,:,:), cnres)
  plot_topo(i) = gsn_csm_contour(wks, tpmask, res_topo)
  overlay(plot(i), plotc(i))
  overlay(plot(i), plot_topo(i))
end do

;=====================================================

pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True


gsn_panel(wks, plot, (/ncase,1/), pres)


end