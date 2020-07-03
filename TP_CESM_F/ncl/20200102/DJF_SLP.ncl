
begin

lonstr = 50
lonend = 130
latstr = 20
latend = 70

casename = (/"CTRL","NOQ1","NOTP_NOQ1"/) ;,"NOTP"
ncase = dimsizes(casename)
; month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
; nm = dimsizes(month)
season = (/"DJF"/)

LEV = (/900/)
z = new((/ncase,96,144/), float)
u = z
v = z

do i = 0, dimsizes(casename)-1
    diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(i)+"_F2000_F19G16_CAM4/pro/"

    filiz = "F2000_"+casename(i)+".PSL.mon.0601-3012.nc"
    fz = addfile(diri+filiz, "r")
    Z1 = fz->PSL
    z1 = month_to_season(Z1, season)
    z(i,:,:) = dim_avg_n_Wrap(z1, 0)
end do
z = z/100
z@units = "hPa"
z!0 = "case"
z&case = casename
printVarSummary(z)

zx = z(0:1,:,:)
zx(0,:,:) = z(0,:,:) - z(1,:,:)
zx(1,:,:) = z(1,:,:) - z(2,:,:)
;----------------------------------------------------------------
tpmaskpath = "/home/yangsong3/wuyt/sysu/scripts/TP_CESM_F/model_need_input_data/"
tpmaskname = "TPmask_1000m_C20190410.nc"
ftp = addfile(tpmaskpath+tpmaskname, "r")

tpmask = ftp->tp
printVarSummary(tpmask)
;----------------------------------------------------------------

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20200102/"

figname = "DJF_SLP"

wks = gsn_open_wks("eps", figpath+figname)
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
   res_topo@cnLineThicknessF       = 4
   res_topo@cnLineDashPattern      = 16
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
    res@cnFillPalette = "WhiteBlue"
;=====================================================

; res@cnLevelSelectionMode = "ExplicitLevels"
; res@cnLevels = (/12460,12470,12480,12500,12520,12540,12560/)
; ;res@cnFillColors = ()
res@cnLevelSelectionMode = "ManualLevels"
res@cnMaxLevelValF = 1040
res@cnMinLevelValF = 1020
res@cnLevelSpacingF = 5

; leftsti = (/"Ctrl","NT","NTNH"/)
; do i = 0, ncase-1
;   res@gsnLeftString = leftsti(i)+"_SLP"
;   plot(i) = gsn_csm_contour_map(wks, z(i,:,:), res)
;   plot_topo(i) = gsn_csm_contour(wks, tpmask, res_topo)
;   overlay(plot(i), plot_topo(i))

; end do


;=====================================================

pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
pres@lbLabelFontHeightF = 0.015
;gsn_panel(wks, plot, (/ncase,1/), pres)

;=====================================================
gres = True
gres@gsLineColor = "green4"
gres@gsLineThicknessF = 5.0
 glon = (/120,120,70,70,120/)
 glat = (/40,60,60,40,40/)



;=====================================================
res@cnLevelSelectionMode = "ManualLevels"
res@cnMaxLevelValF = 6
res@cnMinLevelValF = -6
res@cnLevelSpacingF = 1
res@cnFillPalette = "BlueWhiteOrangeRed"
res@gsnSpreadColors = True
res@gsnSpreadColorEnd   = 0
res@gsnSpreadColorStart = 253 ;(255-2)
Texts = (/"(a) Thermal","(b) Dynamical"/)
do i = 0, 1
  res@gsnLeftString = Texts(i)
  plot(i) = gsn_csm_contour_map(wks, zx(i,:,:), res)
  plot_topo(i) = gsn_csm_contour(wks, tpmask, res_topo)
  overlay(plot(i), plot_topo(i))
end do

dum1 = gsn_add_polyline(wks, plot(0), glon, glat, gres)
dum2 = gsn_add_polyline(wks, plot(1), glon, glat, gres)

gsn_panel(wks, plot(0:1), (/1,2/), pres)
end