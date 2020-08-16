
begin

LEV = 200
nlat = 96
nlon = 144

; lonstr = 50
; lonend = 180
; latstr = -10
; latend = 40
lonstr = 50
lonend = 150
latstr = 0 ;-10
latend = 50

casename = (/"CTRL","NOTPQ1","NOTPQ1_MAR","NOTPQ1_APR","NOTPQ1_MAY"/) ;
ncase = dimsizes(casename)
year = ispan(6, 30, 1)
nyear = dimsizes(year)
month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
nmonth = dimsizes(month)

t = new((/ncase,nyear,nmonth,96,144/), float)

do i = 0, ncase-1
  diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(i)+"_F2000_F19G16_CAM4/pro/"
  tname = "F2000_"+casename(i)+".plevel.T.mon.0601-3012.nc"
  pname = "F2000_"+casename(i)+".PS.mon.0601-3012.nc"
  ft = addfile(diri+tname, "r")
  fp = addfile(diri+pname, "r")
  T = ft->T ;(:,{LEV},:,:)
  PS = fp->PS

  levs = 500
  levt = 200

  plev = T&lev
  plev = plev*100
  plev@units = "Pa"
  print(min(plev))
  dp = dpres_plevel_Wrap(plev, PS, min(plev), 0)
  Tdp = T
  Tdp = T*dp
  tdp = dim_sum_n_Wrap(Tdp(:,{levs:levt},:,:), 1)/dim_sum_n_Wrap(dp(:,{levs:levt},:,:), 1)
  copy_VarCoords(T(:,0,:,:), tdp)

  t(i,:,:,:,:) = reshape(tdp, (/nyear,nmonth,96,144/))
  copy_VarCoords(tdp(0,:,:), t(i,0,0,:,:))
end do
t!0 = "case"
t&case = casename
t!1 = "year"
t&year = year
t!2 = "month"
t&month = month
printVarSummary(t)

;------------ttest--------------------------------------

t_var = dim_variance_n_Wrap(t, 1)
t_ave = dim_avg_n_Wrap(t, 1)
prob = t_ave(0:ncase-2,:,:,:)
t_diff = t_ave(0:ncase-2,:,:,:)
do i = 0, ncase-2
    prob(i,:,:,:) = ttest(t_ave(0,:,:,:), t_var(0,:,:,:), nyear, t_ave(i+1,:,:,:), t_var(i+1,:,:,:), nyear, False, False)
    t_diff(i,:,:,:) = t_ave(0,:,:,:) - t_ave(i+1,:,:,:)
end do
printVarSummary(prob)
; ; print(t_diff(2,2,:,:))
; print(t_ave(2,2,:,:))
; exit
case2 = "Ctrl -"+casename(1:ncase-1)
prob!0 = "case"
prob&case = case2
copy_VarCoords(prob, t_diff)
printVarSummary(t_diff)


;----------------------------------------------------------------
tpmaskpath = "/home/yangsong3/wuyt/sysu/scripts/TP_CESM_F/model_need_input_data/"
; tpmaskname = "TPmask_1000m_C20190410.nc"
tpmaskname = "TPmask_1500m_C20200713.nc"
ftp = addfile(tpmaskpath+tpmaskname, "r")

tpmask = ftp->tp
printVarSummary(tpmask)
;-------------------------------------------------------------------
;-------------------------------------------------------------------
Month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20200704/"
figname = "TTT_diff_multimodel"

wks = gsn_open_wks("pdf", figpath+figname)

nplot = 4*(ncase-1)
plot = new(nplot, graphic)
plot_prob = plot
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
 res_topo@cnLineColor            = "blue" ; gray25
 res_topo@cnInfoLabelOn          = False
 res_topo@cnLevelSelectionMode   ="ExplicitLevels"
 res_topo@cnLevels               = (/0.5/)
 res_topo@gsnLeftString = ""
 res_topo@gsnRightString = ""
;----------------------------------------------------

res = True
res@gsnDraw  = False
res@gsnFrame = False

res@mpOutlineOn = True
res@mpGeophysicalLineThicknessF = 2
res@mpGeophysicalLineColor = "grey25" ;"red" 
res@mpFillDrawOrder = "PreDraw"
res@mpFillOn = True

res@mpMinLatF = latstr
res@mpMaxLatF = latend
res@mpMinLonF = lonstr
res@mpMaxLonF = lonend

res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False

res@lbLabelBarOn = False
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 10
res@cnMaxLevelValF = 100
res@cnMinLevelValF = -100
res@cnFillPalette = "CBR_coldhot"

res@cnInfoLabelOn = False
;res@cnLineDrawOrder = "PostDraw"
; res@cnLabelDrawOrder = "PostDraw"
;res@cnLevelFlag = "LineAndLabel"
;res@cnLineLabelInterval = 1

res@gsnStringFontHeightF = 0.035
res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025

;------------------------
; draw prob
;------------------------
res3  = True
res3@gsnDraw   = False
res3@gsnFrame  = False

;res3@cnFillOn  = False
res3@cnLinesOn = False
res3@cnLineLabelsOn = False
res3@cnInfoLabelOn  = False
;res3@lbLabelBarOn   = False

res3@gsnRightString = ""

res3@cnLevelSelectionMode = "ManualLevels"
res3@cnMinLevelValF = 0.00
res3@cnMaxLevelValF = 0.10
res3@cnLevelSpacingF = 0.01

;res3@cnFillPattern = 17 ;(/17/)
res3@cnFillColor  =  "black" ;"green4"
res3@gsnLeftString = ""
res3@cnFillDotSizeF = 0.003
res3@cnFillScaleF  = 1 ;1.2


do i = 0, ncase-2
    do j = 0, 3
        res@gsnLeftString = Month(j+2)
        res@gsnRightString = case2(i)
        plot(i*4+j) = gsn_csm_contour_map_ce(wks, t_diff(i,j+2,:,:), res)
        plot_prob(i*4+j) = gsn_csm_contour(wks, prob(i,j+2,:,:), res3)
        overlay(plot(i*4+j), plot_prob(i*4+j))
        plot_topo(i*4+j) = gsn_csm_contour(wks, tpmask, res_topo)
        overlay(plot(i*4+j), plot_topo(i*4+j))
    end do
end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = "Tropo Temp ("+levs+"hPa - "+levt+"hPa)"
;pres@pmLabelBarOrthogonalPosF = -0.01

gsn_panel(wks, plot, (/4,4/), pres)

end
