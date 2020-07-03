
begin

lonstr = 0
lonend = 160
latstr = 0
latend = 60

casename = (/"CTRL","NOQ1","NOTP_NOQ1"/) ;,"NOTP"
ncase = dimsizes(casename)
; month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
; nm = dimsizes(month)
season = (/"JJA"/)

LEV = (/200,500/)
z = new((/ncase,2,96,144/), float)

do i = 0, dimsizes(casename)-1
    diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(i)+"_F2000_F19G16_CAM4/pro/"
    filiz = "F2000_"+casename(i)+".plevel.Z3.mon.0601-3012.nc"

    fz = addfile(diri+filiz, "r")
    Z1 = fz->Z3(:,{LEV},:,:)

    z1 = month_to_season(Z1, season)
    z(i,:,:,:) = dim_avg_n_Wrap(z1, 0)
end do
; z = z/9.8
; z@units = "gpm"
z!0 = "case"
z&case = casename
printVarSummary(z)


;----------------------------------------------------------------
tpmaskpath = "/home/yangsong3/wuyt/sysu/scripts/TP_CESM_F/model_need_input_data/"
tpmaskname = "TPmask_1000m_C20190410.nc"
ftp = addfile(tpmaskpath+tpmaskname, "r")

tpmask = ftp->tp
printVarSummary(tpmask)
;----------------------------------------------------------------

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20200102/"

figname = "JJA_SAH_WPSH"

wks = gsn_open_wks("eps", figpath+figname)
nplot = 4
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
    res@pmLabelBarOrthogonalPosF = 0.12
    ;res@lbLabelBarOn = False
    res@lbLabelFontHeightF = 0.02
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

res@cnLevelSelectionMode = "ExplicitLevels"
res@cnLevels = (/12460,12480,12500,12520,12540,12560/)
;res@cnFillColors = ()

cnres@cnLevelSelectionMode = "ExplicitLevels"
cnres@cnLevels = (/12460,12480,12500,12520,12540,12560/) ;(/12460,12470,12480,12500,12520,12540,12560/)


res@gsnLeftString = "(a) H200 (NH & Ctrl)"
plot(0) = gsn_csm_contour_map(wks, z(1,{LEV(0)},:,:), res)
plotc(0) = gsn_csm_contour(wks, z(0,{LEV(0)},:,:), cnres)
plot_topo(0) = gsn_csm_contour(wks, tpmask, res_topo)
overlay(plot(0), plotc(0))
overlay(plot(0), plot_topo(0))

res@gsnLeftString = "(b) H200 (NTNH & NH)"
plot(1) = gsn_csm_contour_map(wks, z(2,{LEV(0)},:,:), res)
plotc(1) = gsn_csm_contour(wks, z(1,{LEV(0)},:,:), cnres)
plot_topo(1) = gsn_csm_contour(wks, tpmask, res_topo)
overlay(plot(1), plotc(1))
overlay(plot(1), plot_topo(1))
;=====================================================

lonstr = 70
lonend = 200
latstr = 0
latend = 60

res@mpMinLatF = latstr
res@mpMaxLatF = latend
res@mpMinLonF = lonstr
res@mpMaxLonF = lonend
res@mpCenterLonF = (lonstr+lonend)/2

res@cnLevelSelectionMode = "ManualLevels"
res@cnMaxLevelValF = 5900
res@cnMinLevelValF = 5885
res@cnLevelSpacingF = 5

cnres@cnLevelSelectionMode = "ManualLevels"
cnres@cnMaxLevelValF = 5900 ;5905
cnres@cnMinLevelValF = 5885
cnres@cnLevelSpacingF = 5

res@gsnLeftString = "(c) H500 (NH & Ctrl)"
plot(2) = gsn_csm_contour_map(wks, z(1,{LEV(1)},:,:), res)
plotc(2) = gsn_csm_contour(wks, z(0,{LEV(1)},:,:), cnres)
plot_topo(2) = gsn_csm_contour(wks, tpmask, res_topo)
overlay(plot(2), plotc(2))
overlay(plot(2), plot_topo(2))

res@gsnLeftString = "(d) H500 (NTNH & NH)"
plot(3) = gsn_csm_contour_map(wks, z(2,{LEV(1)},:,:), res)
plotc(3) = gsn_csm_contour(wks, z(1,{LEV(1)},:,:), cnres)
plot_topo(3) = gsn_csm_contour(wks, tpmask, res_topo)
overlay(plot(3), plotc(3))
overlay(plot(3), plot_topo(3))


;=====================================================
pres = True
pres@gsnMaximize = True

gsn_panel(wks, plot, (/2,2/), pres)


end