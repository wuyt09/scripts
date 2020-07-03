begin

lonstr = 50
lonend = 120
latstr = 20
latend = 50

casename = (/"CTRL"/) ;"NOTP",,"NOTP_NOQ1"
month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)

Cpair = 1004.07
G = 9.8 

i = 0
diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(i)+"_F2000_F19G16_CAM4/pro/"
filiu = "F2000_"+casename(i)+".plevel.DTCOND.mon.0601-3012.nc"
filiv = "F2000_"+casename(i)+".plevel.DTV.mon.0601-3012.nc"
filis = "F2000_"+casename(i)+".plevel.QRS.mon.0601-3012.nc"
filil = "F2000_"+casename(i)+".plevel.QRL.mon.0601-3012.nc"
filips = "F2000_"+casename(i)+".PS.mon.0601-3012.nc"

fu = addfile(diri+filiu, "r")
fv = addfile(diri+filiv, "r")
fs = addfile(diri+filis, "r")
fl = addfile(diri+filil, "r")
fps = addfile(diri+filips, "r")

U = fu->DTCOND
V = fv->DTV
S = fs->QRS
L = fl->QRL
PS = fps->PS

Q1 = U
Q1 = (U+V+S+L)*Cpair/G

plev = U&lev*100
plev!0 = "lev"
plev&lev = U&lev
plev@units = "Pa"
psfc = PS
ptop = min(plev)
dp = dpres_plevel(plev, psfc, ptop, 0)

Q1dp = Q1
Q1dp = Q1*dp
q1dp = dim_sum_n_Wrap(Q1dp, 1)
copy_VarCoords(Q1(:,0,:,:), q1dp)

q1 = clmMonTLL(q1dp)
printVarSummary(q1)

;==================================================================

tpmaskpath = "/home/yangsong3/wuyt/sysu/scripts/TP_CESM_F/model_need_input_data/"
tpmaskname = "TPmask_1000m_C20190410.nc"
ftp = addfile(tpmaskpath+tpmaskname, "r")

tpmask = ftp->tp
printVarSummary(tpmask)


q1 = mask(q1, tpmask.ge.0.9, True)
q1_series = wgt_areaave_Wrap(q1(:,{latstr:latend},{lonstr:lonend}), 1, 1, 0)
q1_std = dim_stddev(q1_series)

;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20191229/"
figname = "Clim_monthly_Q1_Ctrl"

; lonstr = 0
; lonend = 360
; latstr = -30
; latend = 80

wks = gsn_open_wks("pdf", figpath+figname)
nplot = 12
plot = new(nplot, graphic)
plot_topo = plot

 res_topo            = True
 res_topo@cnFillOn   = False
 res_topo@cnLinesOn  = True
 res_topo@gsnFrame   = False
 res_topo@gsnDraw    = False
 res_topo@cnLineLabelsOn         = False
 res_topo@cnLineThicknessF       = 2
 res_topo@cnLineColor            = "grey80"
 res_topo@cnLineDashPattern      = 12
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
res@cnLevelSpacingF = 30.
res@cnMaxLevelValF = 300.
res@cnMinLevelValF = -300.
res@cnFillPalette = "BlueWhiteOrangeRed"

do i = 0, 11
  res@gsnLeftString = month(i)
  plot(i) = gsn_csm_contour_map(wks, q1(i,:,:), res)
  plot_topo(i) = gsn_csm_contour(wks, tpmask, res_topo)
  overlay(plot(i), plot_topo(i))
end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelMainString = "Q1_climatology"
pres@gsnPanelLabelBar = True
;pres@pmLabelBarOrthogonalPosF = -0.01
pres@lbLabelFontHeightF = 0.01

gsn_panel(wks, plot, (/3,4/), pres)

;=========================================================
delete(res)
Month = ispan(1, 12, 1)
month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)

res = True
res@gsnDraw  = False
res@gsnFrame = False
res@gsnMaximize = True

res@vpWidthF = 0.7
res@vpHeightF = 0.4

res@trXMinF = 1
res@trXMaxF = 12
res@tmXBMode = "Explicit"
res@tmXBValues = Month
res@tmXBLabels = month
res@tmXTOn = False
;res@tiXAxisString = "Month"

; res@trYMinF = -0.08
; res@trYMaxF = 0.04
; res@tmYLMode = "Explicit"
; res@tmYLValues = ispan(-, 0, 1)
; res@tmYLLabels = ispan(-5, 0, 1)
res@tmYROn = False
res@tiYAxisString = " "
res@gsnYRefLine = (/-q1_std,0,q1_std/)
res@gsnYRefLineDashPattern = 16
res@gsnYRefLineColor = "gray20"
res@gsnYRefLineThicknessF = 5

Colors = (/"black","red","blue","forestgreen","purple","goldenrod3","lightblue"/)
Patterns = (/0,0,0,0,16,16,16/)
res@xyLineThicknessF = 8
res@xyLineColors = Colors
res@xyDashPatterns = Patterns

res@gsnLeftString = "Q1"

plot(0) = gsn_csm_xy(wks, Month, q1_series, res)
draw(plot(0))
end