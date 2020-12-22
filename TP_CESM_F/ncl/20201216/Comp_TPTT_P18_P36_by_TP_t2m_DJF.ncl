
begin

yrStrt = 1979
yrLast = 2019
mnStrt = 1
mnLast = 12

year = ispan(yrStrt, yrLast, 1)
month = ispan(mnStrt, mnLast, 1)
penta = ispan(1, 73, 1)
nyear = dimsizes(year)
nmonth = dimsizes(month)
npenta = dimsizes(penta)
season = (/"MAM"/)

fx = "time_series_TP_t2m_DJF_MAM_OBS.txt"
x = asciiread(fx, (/nyear,3/), "float")
inx = x(0:nyear-3,1)

hig_ind = ind(inx.ge.1.)
low_ind = ind(inx.le.-1.)
print(year(hig_ind))
print(year(low_ind))

hig_ind = hig_ind+1
low_ind = low_ind+1

;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)
printVarSummary(topo)
;====================================================================
latstr = 25
latend = 40
lonstr = 75
lonend = 105
pens = 18
pene = 36

diri = "/home/yangsong3/wuyt/sysu/scripts/monsoon_onset_pentad/"
fili = "1_Tropo_Temp_pentad_ERA5_1979_2018.nc"
f = addfile(diri+fili, "r")

tropot = f->tropoT
printVarSummary(tropot)
tropot = mask(tropot, topo.ge.1500., True)

TPTT = wgt_areaave_Wrap(tropot(:,:,{latstr:latend},{lonstr:lonend}), 1., 1., 0)

var_hig = dim_avg_n_Wrap(TPTT(hig_ind,pens-1:pene-1), 0)
var_low = dim_avg_n_Wrap(TPTT(low_ind,pens-1:pene-1), 0)
printVarSummary(var_hig)

lines = new((/2,pene-pens+1/), float)
lines(0,:) = var_hig
lines(1,:) = var_low
;====================================================================
; plot
;====================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20201216/"
figname = "Comp_TPTT_P18_P36_by_TP_t2m_DJF"
wks = gsn_open_wks("pdf", figpath+figname)
plot = new(1, graphic)

Colors = (/"red","blue"/)
Patterns = (/0,0/)

res = True
res@gsnFrame = False
res@gsnDraw  = False
res@gsnMaximize = True

res@vpWidthF = 0.7
res@vpHeightF = 0.45

res@xyLineThicknessF = 10
res@xyLineColors = Colors
res@xyDashPatterns = Patterns

res@trXMinF = pens
res@trXMaxF = pene
res@tmXBPrecision = 2

res@gsnYRefLine = (/0,1,-1/)
res@gsnYRefLineColor = "grey40"
res@gsnYRefLineThicknessF = 5
res@gsnYRefLineDashPattern = 2

res@gsnRightString = ""
res@gsnStringFontHeightF = 0.025

plot(0) = gsn_csm_xy(wks, penta(pens-1:pene-1), lines, res)
draw(plot)
end