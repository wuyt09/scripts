begin

latstr = -10
latend = 70
lonstr = 10
lonend = 150
yrStrt = 1979
yrLast = 2016
Season = (/"MAM","JJA","SON","DJF"/)
ns = dimsizes(Season)

diri1 = "~/data-observation/GPCP/"
oname = "precip.mon.mean.nc"
diri2 = "/home/yangsong3/data-model/wuyt/TP_CTRL_F2000_F19G16_CAM4/pro/"
cname = "F2000_CTRL.PRECC.mon.0601-3012.nc"
lname = "F2000_CTRL.PRECL.mon.0601-3012.nc"

f1 = addfile(diri1+oname, "r")
f2 = addfile(diri2+cname, "r")
f3 = addfile(diri2+lname, "r")

Time = f1->time
YYYYMM = cd_calendar(Time, -1)/100
iYYYYMM = ind(YYYYMM.ge.yrStrt .and. YYYYMM.le.yrLast)

x1 = f1->precip(iYYYYMM,{latstr:latend},{lonstr:lonend}) ; mm/day
y1 = f2->PRECC(:,{latstr:latend},{lonstr:lonend})
y2 = f3->PRECL(:,{latstr:latend},{lonstr:lonend})
x2  = y1 + y2  ; m/s
x2  = x2*1000*86400
copy_VarMeta(y1, x2)

x10 = month_to_seasonN(x1, Season)
x20 = month_to_seasonN(x2, Season)
x11 = dim_avg_n_Wrap(x10, 1)   ; year
x22 = dim_avg_n_Wrap(x20, 1)   ; year
printVarSummary(x11)
printVarSummary(x22)

x12 = f2fsh(x22, (/dimsizes(x11&lat),dimsizes(x11&lon)/))  ;dimsizes(x22)
printVarSummary(x12)
diff = x12 - x11
copy_VarMeta(x11, diff)
printVarSummary(diff)

pp = new((/3,ns,dimsizes(x11&lat),dimsizes(x11&lon)/), float)
pp(0,:,:,:) = x11
pp(1,:,:,:) = x12
pp(2,:,:,:) = diff
copy_VarMeta(x11, pp(0,:,:,:))
printVarSummary(pp)

;**************************************************************
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/CTRL_RUN/"
figname = "draw_diff_obs_ctrlrun_seasonal_prec"

wks = gsn_open_wks("pdf",figpath+figname)
nplot = 12
plot = new(nplot,graphic)
plot_topo = plot
gsn_define_colormap(wks,"precip_11lev")

 ;------------------------------------------------------------
 ; TOPO
 ;------------------------------------------------------------
 
 ff = "~/wuyt/data/TOPO_1500.nc"
 topo_in = addfile(ff,"r")
 topo = topo_in->TP

 xlat = topo_in->LATIXY(:,0)
 xlon = topo_in->LONGXY(0,:)
 topo!0 = "lsmlat"
 topo!1 = "lsmlon"
 topo&lsmlat = xlat
 topo&lsmlon = xlon
 
 res_topo            = True
 res_topo@cnFillOn   = False
 res_topo@cnLinesOn  = True
 res_topo@gsnFrame   = False
 res_topo@gsnDraw    = False
 res_topo@cnLineLabelsOn         = False
 res_topo@cnLineThicknessF       = 5
 res_topo@cnLineColor            = "black"
 res_topo@cnInfoLabelOn          = False
 res_topo@cnLevelSelectionMode   ="ExplicitLevels"
 res_topo@cnLevels               = (/1500/)

;---------------------------------------

res                 =  True
res@gsnDraw         = False
res@gsnFrame        = False
res@gsnAddCyclic    = False
;res@gsnSpreadColors = True

res@mpOutlineOn  = True
res@mpFillOn     = True
;res@mpCenterLonF = 0.
res@mpMinLatF    = latstr
res@mpMaxLatF    = latend
res@mpMinLonF    = lonstr
res@mpMaxLonF    = lonend
;res@mpCenterLonF = 180
res@mpGeophysicalLineThicknessF = 2

res@cnFillOn     = True
res@cnLinesOn    = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn  = False
;res@lbLabelBarOn   = False

;----------------------------------------

res@gsnRightString = "mm/day"
res@gsnStringFontHeightF = 0.03
;res@tmXBLabelFontHeightF = 0.018
;res@tmYLLabelFontHeightF = 0.018

do i = 0, ns-1

res@gsnRightString = Season(i)
res@cnFillPalette = "precip_11lev"
res@cnLevelSelectionMode = "ExplicitLevels"
;res@cnLevels = (/3,6,9,12,15,18,21,24,27,30/)
res@cnLevels = (/0.2,1,2,4,7,12,15,18,20/)

res@gsnLeftString = "Obs"
;plot(0) = gsn_csm_contour_map(wks,x11,res)
plot(i) = gsn_csm_contour_map(wks, pp(0,i,:,:), res)
plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(i), plot_topo(i))

res@gsnLeftString = "Ctrl_run"
;plot(1) = gsn_csm_contour_map(wks,x22,res)
plot(i+4) = gsn_csm_contour_map(wks, pp(1,i,:,:), res)
plot_topo(i+4) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(i+4), plot_topo(i+4))

res@gsnLeftString = "Ctrl_run - Obs"
res@cnFillPalette = "precip_diff_12lev"
;delete(res@cnLevelSpacingF)
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF  = -5
res@cnMaxLevelValF  = 5
res@cnLevelSpacingF = 0.5
;plot(2) = gsn_csm_contour_map(wks,diff,res)
plot(i+8) = gsn_csm_contour_map(wks, pp(2,i,:,:), res)
plot_topo(i+8) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(i+8), plot_topo(i+8))
end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelMainString = "Total_precip"
;pres@gsnPanelLabelBar = True
;pres@gsnPanelLeft
gsn_panel(wks, plot, (/3,4/), pres)

end