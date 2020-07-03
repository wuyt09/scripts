
begin

LEV = 200
nlat = 96
nlon = 144

; lonstr = 50
; lonend = 180
; latstr = -10
; latend = 40
lonstr = 50
lonend = 130
latstr = 20 ;-10
latend = 50

casename = (/"CTRL","NOQ1"/) ;"NOTP"
ncase = dimsizes(casename)
year = ispan(6, 30, 1)
nyear = dimsizes(year)
month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
nmonth = dimsizes(month)
season = "JJA"
LEV = 850
z = new((/ncase,96,144/), float)
u = z
v = z

do i = 0, ncase-1
  diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(i)+"_F2000_F19G16_CAM4/pro/"
  zname = "F2000_"+casename(i)+".PREC.mon.0601-3012.nc"
  fz = addfile(diri+zname, "r")
  Z = fz->PREC
  Z = Z*1000*86400
  Z@units = "mm/day"
  Z_sea = month_to_season(Z, season)
  z(i,:,:) = dim_avg_n_Wrap(Z_sea, 0)
end do
z!0 = "case"
z&case = casename
printVarSummary(z)

ratio = z(0,:,:)
ratio = (z(1,:,:)-z(0,:,:))/z(0,:,:)
printVarSummary(ratio)
printMinMax(ratio, False)
; exit
;----------------------------------------------------------------
tpmaskpath = "/home/yangsong3/wuyt/sysu/scripts/TP_CESM_F/model_need_input_data/"
tpmaskname = "TPmask_1000m_C20190410.nc"
ftp = addfile(tpmaskpath+tpmaskname, "r")

tpmask = ftp->tp
printVarSummary(tpmask)
;-------------------------------------------------------------------
;-------------------------------------------------------------------

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20200513/"
figname = season+"_Precip_ratio"

wks = gsn_open_wks("pdf", figpath+figname)

nplot = 1
plot = new(nplot, graphic)
plotv = plot
plotc = plot
plot_topo = plot
;----------------------------------------------------
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
;----------------------------------------------------

res = True
res@gsnDraw  = False
res@gsnFrame = False
res@gsnMaximize = True

res@mpOutlineOn = True
res@mpGeophysicalLineThicknessF = 2
res@mpGeophysicalLineColor = "grey25" ;"red" 
res@mpFillDrawOrder = "PreDraw"
; res@mpFillOn = True

res@mpMinLatF = latstr
res@mpMaxLatF = latend
res@mpMinLonF = lonstr
res@mpMaxLonF = lonend

res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@lbLabelBarOn = True
res@cnInfoLabelOn = False
res@pmLabelBarOrthogonalPosF = 0.15

res@gsnStringFontHeightF = 0.03
res@tmXBLabelFontHeightF = 0.02
res@tmYLLabelFontHeightF = 0.02


res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 0.1 ;0.5
res@cnMaxLevelValF = 1.
res@cnMinLevelValF = -1.
res@cnFillPalette = "CBR_drywet"



res@gsnLeftString = season

plot(0) = gsn_csm_contour_map(wks, ratio, res)
plot_topo(0) = gsn_csm_contour(wks, tpmask, res_topo)
overlay(plot(0), plot_topo(0))

draw(plot)


end