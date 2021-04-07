begin

lonstr = 75 ; 72 ;75
lonend = 105 ;136 ;105
latstr = 25 ;18 ;25
latend = 40 ;54 ; 40


;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)
printVarSummary(topo)

;===================================================================
; OBS grid data
;===================================================================
yrStrt = 1979
yrLast = 2018
year = ispan(yrStrt, yrLast, 1)
nyear = dimsizes(year)
month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
nm = dimsizes(month)

;===================================================================
ts_obs = new((/nyear*nm,163,283/), float)

diriobs = "/home/yangsong3/data-observation/CN05.1/"

filit = "CN05.1_Tm_1961_2018_month_025x025.nc"
filip = "CN05.1_Pre_1961_2018_month_025x025.nc"

ft = addfile(diriobs+filit, "r")
fp = addfile(diriobs+filip, "r")
Time = ft->time
YYYY = ut_calendar(Time, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
t_time = ft->time(iYYYY)
ts_obs = ft->tm(iYYYY,:,:)
pr_obs = fp->pre(iYYYY,:,:)
printVarSummary(ts_obs)

ts_obs = dtrend_msg_n(t_time, ts_obs, False, False, 0)
pr_obs = dtrend_msg_n(t_time, pr_obs, False, False, 0)

ts_jja = month_to_season(ts_obs, "JJA")
ts_JJA = dim_avg_n_Wrap(ts_jja, 0)
pr_jja = month_to_season(pr_obs, "JJA")
pr_JJA = dim_avg_n_Wrap(pr_jja, 0)
printVarSummary(ts_JJA)


;======================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_surf_temp_EA_prec/20210316/"
figname = "check_t2m_pr_CN051"
wks = gsn_open_wks("pdf", figpath+figname)

plot = new(2, graphic)
plot_topo = plot

 res_topo            = True
 res_topo@cnFillOn   = False
 res_topo@cnLinesOn  = True
 res_topo@gsnFrame   = False
 res_topo@gsnDraw    = False
 res_topo@cnLineLabelsOn         = False
 res_topo@cnLineThicknessF       = 3
 res_topo@cnLineColor            = "grey40"
 res_topo@cnInfoLabelOn          = False
 res_topo@cnLevelSelectionMode   ="ExplicitLevels"
 res_topo@cnLevels               = (/2000/)
 res_topo@gsnLeftString = ""
 res_topo@gsnRightString = ""
;---------------------------------------------------------
res = True
res@gsnDraw = False
res@gsnFrame = False
; res@gsnMaximize = True
res@gsnAddCyclic = False

res@mpOutlineOn = True
res@mpFillOn = True
res@mpGeophysicalLineThicknessF = 2
res@mpGeophysicalLineColor = "grey25" ;"red" 
; res@mpFillDrawOrder = "PostDraw"
; res@mpCenterLonF = 180

res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False
; res@cnFillMode = "CellFill"

res@gsnStringFontHeightF = 0.025
res@tmXBLabelFontHeightF = 0.02
res@tmYLLabelFontHeightF = 0.02

res@mpMinLatF = 22 ;latstr
res@mpMaxLatF = 40 ;latend
res@mpMinLonF = 106 ;lonstr
res@mpMaxLonF = 122 ;lonend
res@cnFillPalette = "WhiteYellowOrangeRed"
; res@cnLevelSelectionMode = "ManualLevels"
; res@cnLevelSpacingF = 0.01
; res@cnMaxLevelValF = 0.05
; res@cnMinLevelValF = -0.05
res@gsnLeftString = "JJA mean t2m"
plot(0) = gsn_csm_contour_map(wks, ts_JJA, res)
plot_topo(0) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(0), plot_topo(0))

res@cnFillPalette = "CBR_wet"
res@gsnLeftString = "JJA mean prec"
plot(1) = gsn_csm_contour_map(wks, pr_JJA, res)
plot_topo(1) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(1), plot_topo(1))

pres = True
pres@gsnMaximize = True
gsn_panel(wks, plot, (/1,2/), pres)

end