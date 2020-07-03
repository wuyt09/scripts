
begin


lonstr = 40
lonend = 150
latstr = 0
latend = 60

casename = (/"CTRL","NOTP"/) ;"NOQ1","NOTP_NOQ1",
ncase = dimsizes(casename)
month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
nm = dimsizes(month)
season = (/"MAM","JJA","SON","DJF"/)
ns = dimsizes(season)

LEV = 200
t = new((/ncase,4,96,144/), float)

do i = 0, dimsizes(casename)-1
    diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(i)+"_F2000_F19G16_CAM4/pro/"
    filit = "F2000_"+casename(i)+".plevel.T.mon.0601-3012.nc"
    filips = "F2000_"+casename(i)+".PS.mon.0601-3012.nc"

    ft = addfile(diri+filit, "r")
    T = ft->T
    fps = addfile(diri+filips, "r")
    PS = fps->PS

    plev = T&lev*100
    plev!0 = "lev"
    plev&lev = T&lev
    plev@units = "Pa"
    psfc = PS
    ptop = min(plev)
    dp = dpres_plevel_Wrap(plev, psfc, ptop, 0)
    Tdp = T
    Tdp = T*dp
    t1 = dim_sum_n_Wrap(Tdp(:,{400:200},:,:), 1)/dim_sum_n_Wrap(dp(:,{400:200},:,:), 1)
    copy_VarCoords(T(:,0,:,:), t1)

    TSUM = month_to_seasonN(t1, season)
    t(i,:,:,:) = dim_avg_n_Wrap(TSUM, 1)
end do
;t = t - 273.15
tx = t(0,:,:,:)
tx = t(1,:,:,:) - t(0,:,:,:)

;----------------------------------------------------------------
tpmaskpath = "/home/yangsong3/wuyt/sysu/scripts/TP_CESM_F/model_need_input_data/"
tpmaskname = "TPmask_1000m_C20190410.nc"
ftp = addfile(tpmaskpath+tpmaskname, "r")

tpmask = ftp->tp
printVarSummary(tpmask)
;==================================================================
;==================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20200106/"
figname = "TropoT_Ctrl_TPforcing"

wks = gsn_open_wks("pdf", figpath+figname)
nplot = 4*2
plot = new(nplot, graphic)
plot_topo = plot
plotv = plot


 res_topo            = True
 res_topo@cnFillOn   = False
 res_topo@cnLinesOn  = True
 res_topo@gsnFrame   = False
 res_topo@gsnDraw    = False
 res_topo@cnLineLabelsOn         = False
 res_topo@cnLineThicknessF       = 4
 res_topo@cnLineDashPattern      = 16
 res_topo@cnLineColor            = "gray"
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

res@lbLabelBarOn = False
res@cnLevelSelectionMode = "ManualLevels"

do i = 0, ns-1

    res@gsnRightString = season(i)
    res@gsnLeftString = "Ctrl"
    res@cnLevelSpacingF = 5.
    res@cnMaxLevelValF = 260.
    res@cnMinLevelValF = 220.
    res@cnFillPalette = "WhiteYellowOrangeRed"
    plot(i) = gsn_csm_contour_map(wks, t(0,i,:,:), res)
    plot_topo(i) = gsn_csm_contour(wks, tpmask, res_topo)
    overlay(plot(i), plot_topo(i))

    res@gsnLeftString = "NT - Ctrl"
    res@cnLevelSpacingF = 0.4
    res@cnMaxLevelValF = 2.
    res@cnMinLevelValF = -2.
    res@cnFillPalette = "BlueWhiteOrangeRed"
    plot(i+4) = gsn_csm_contour_map(wks, tx(i,:,:), res)
    plot_topo(i+4) = gsn_csm_contour(wks, tpmask, res_topo)
    overlay(plot(i+4), plot_topo(i+4))

end do

pres = True
pres@gsnMaximize = True
;pres@gsnPanelMainString = "UV"+LEV+"_Precip"
pres@gsnPanelLabelBar = True
;pres@pmLabelBarOrthogonalPosF = -0.01
pres@lbLabelFontHeightF = 0.01

gsn_panel(wks, plot(0:3), (/1,4/), pres)
gsn_panel(wks, plot(4:7), (/1,4/), pres)






end