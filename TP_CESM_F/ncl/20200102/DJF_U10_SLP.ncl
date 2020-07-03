
begin

lonstr = 70
lonend = 150
latstr = 20
latend = 80

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

    filiu = "F2000_"+casename(i)+".plevel.U.mon.0601-3012.nc"
    fu = addfile(diri+filiu, "r")
    U1 = fu->U(:,{LEV},:,:)
    u1 = month_to_season(U1, season)
    u(i,:,:) = dim_avg_n_Wrap(u1, 0)

    filiv = "F2000_"+casename(i)+".plevel.V.mon.0601-3012.nc"
    fv = addfile(diri+filiv, "r")
    V1 = fv->V(:,{LEV},:,:)
    v1 = month_to_season(V1, season)
    v(i,:,:) = dim_avg_n_Wrap(v1, 0)
end do
z = z/100
z@units = "hPa"
z!0 = "case"
z&case = casename
printVarSummary(z)
u!0 = "case"
u&case = casename
v!0 = "case"
v&case = casename

;----------------------------------------------------------------
tpmaskpath = "/home/yangsong3/wuyt/sysu/scripts/TP_CESM_F/model_need_input_data/"
tpmaskname = "TPmask_1000m_C20190410.nc"
ftp = addfile(tpmaskpath+tpmaskname, "r")

tpmask = ftp->tp
printVarSummary(tpmask)
;----------------------------------------------------------------

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20200102/"

figname = "DJF_U10_SLP"

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
    res@cnFillPalette = "WhiteBlue"

;=====================================================
    vcres                         = True             ; plot mods desired
    vcres@gsnFrame                = False            ; so we can draw time stamp
    vcres@gsnDraw                 = False
    vcres@vcRefAnnoOrthogonalPosF =  -0.7             ; move ref vector up

    ;vcres@vcGlyphStyle            = "CurlyVector"  ;"LineArrow"    ; turn on curly vectors
    ;vcres@vcMinDistanceF          = 0.017
    vcres@vcGlyphStyle            = "LineArrow"
    vcres@vcMinDistanceF          = 0.03
    vcres@vcRefLengthF            = 0.04            ; define length of vec ref
    
    vcres@vcLineArrowThicknessF   = 3.0
    ;vecres@vcLineArrowHeadMinSizeF = 0.008
    ;vecres@vcLineArrowHeadMaxSizeF = 0.3

    vcres@vcRefAnnoString2On = False
    vcres@vcRefAnnoFontHeightF = 0.02

;=====================================================

; res@cnLevelSelectionMode = "ExplicitLevels"
; res@cnLevels = (/12460,12470,12480,12500,12520,12540,12560/)
; ;res@cnFillColors = ()
res@cnLevelSelectionMode = "ManualLevels"
res@cnMaxLevelValF = 1050
res@cnMinLevelValF = 1030
res@cnLevelSpacingF = 5


vcres@vcRefMagnitudeF         = 5.0  ; define vector ref mag


leftsti = (/"Ctrl","NT","NTNH"/)
do i = 0, ncase-1
  res@gsnLeftString = leftsti(i)+"_Usfc_SLP"
  plot(i) = gsn_csm_contour_map(wks, z(i,:,:), res)
  plotc(i) = gsn_csm_vector(wks, u(i,:,:), v(i,:,:), vcres)
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