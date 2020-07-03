
begin

lonstr = 0
lonend = 160
latstr = 0
latend = 60

casename = (/"CTRL","NOQ1"/) ;,"NOTP_NOQ1","NOTP"
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

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20200513/"

figname = "JJA_SAH_WPSH"

wks = gsn_open_wks("pdf", figpath+figname)
nplot = 2
plot = new(nplot, graphic)
plotc = plot
plotd = plot
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
    res@tmXBLabelFontHeightF = 0.015
    res@tmYLLabelFontHeightF = 0.015
    res@pmLabelBarOrthogonalPosF = 0.1
    res@lbLabelBarOn = True
    res@lbLabelFontHeightF = 0.02
    res@cnFillPalette = "WhiteYellowOrangeRed"

;=====================================================
    cnres = True
    cnres@gsnDraw = False
    cnres@gsnFrame = False

    cnres@cnFillOn = False
    cnres@cnLinesOn = True
    cnres@cnLineLabelsOn = False
    cnres@cnLineLabelFontColor = "red"
    cnres@cnLineThicknessF = 5
    cnres@cnInfoLabelOn = True
    cnres@cnInfoLabelOrthogonalPosF = -0.5
    ;cnres@cnLineDrawOrder = "PostDraw"
    cnres@cnLabelDrawOrder = "PostDraw"
    ;cnres@cnLevelFlag = "LineAndLabel"
    cnres@cnLineLabelInterval = 1

;=====================================================

cnres@cnLevelSelectionMode = "ExplicitLevels"
cnres@cnLevels = (/12460,12480,12500,12520,12540,12560/)
;res@cnFillColors = ()

res@cnLevelSelectionMode = "ExplicitLevels"
res@cnLevels = (/12460,12480,12500,12520,12540,12560/) ;(/12460,12470,12480,12500,12520,12540,12560/)


res@gsnLeftString = "(a) H200 (Ctrl & No_Heat & NoTP)"
plot(0) = gsn_csm_contour_map(wks, z(0,{LEV(0)},:,:), res)
cnres@cnLineColor = "green4"
plotc(0) = gsn_csm_contour(wks, z(1,{LEV(0)},:,:), cnres)
; cnres@cnLineColor = "blue"
; plotd(0) = gsn_csm_contour(wks, z(2,{LEV(0)},:,:), cnres)
plot_topo(0) = gsn_csm_contour(wks, tpmask, res_topo)
overlay(plot(0), plotc(0))
; overlay(plot(0), plotd(0))
overlay(plot(0), plot_topo(0))
;=====================================================

lonstr = 70
lonend = 200
latstr = 0
latend = 50

; res@vpHeightF = 0.5
res@mpMinLatF = latstr
res@mpMaxLatF = latend
res@mpMinLonF = lonstr
res@mpMaxLonF = lonend
res@mpCenterLonF = (lonstr+lonend)/2

res@cnLevelSelectionMode = "ManualLevels"
res@cnMaxLevelValF = 5900
res@cnMinLevelValF = 5885
res@cnLevelSpacingF = 5
res@lbLabelBarOn = True

cnres@cnLevelSelectionMode = "ManualLevels"
cnres@cnMaxLevelValF = 5900 ;5905
cnres@cnMinLevelValF = 5885
cnres@cnLevelSpacingF = 5

res@gsnLeftString = "(b) H500 (Ctrl & No_Heat & NoTP)"
plot(1) = gsn_csm_contour_map(wks, z(0,{LEV(1)},:,:), res)
cnres@cnLineColor = "green4"
plotc(1) = gsn_csm_contour(wks, z(1,{LEV(1)},:,:), cnres)
; cnres@cnLineColor = "blue"
; plotd(1) = gsn_csm_contour(wks, z(2,{LEV(1)},:,:), cnres)
plot_topo(1) = gsn_csm_contour(wks, tpmask, res_topo)
overlay(plot(1), plotc(1))
; overlay(plot(1), plotd(1))
overlay(plot(1), plot_topo(1))
;=====================================================
pres = True
pres@gsnMaximize = True
; pres@gsnPanelLabelBar = False

gsn_panel(wks, plot, (/2,1/), pres)

end