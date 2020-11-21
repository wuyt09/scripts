
begin

ism_onset = (/33,31,30,30,33,30,30,31,31,29,31,27,32,33,30,30,33,31,33,32,28,30,32,32,32,27,32,29,31,30,31,29,30,32,30,32,32,31,31/)
printVarSummary(ism_onset)

yrStrt = 1979
yrLast = 2017
year = ispan(yrStrt, yrLast, 1)
nyear = dimsizes(year)
month = ispan(1, 12, 1)
nmonth = dimsizes(month)

season = (/"DJF","MAM"/)
nlat = 181
nlon = 360

sdiri = "/home/yangsong3/data-observation/ERA5-daily/Single-level/"
T2M = new((/nyear*nmonth,nlat,nlon/), float)
do j = 0, nyear-1
    YEAR = sprinti("%0.4i", year(j))
    ; print(YEAR)
    do i = 0, nmonth-1
        MONTH = sprinti("%0.2i", month(i))
        ; print(MONTH)
        filist = systemfunc("ls "+sdiri+"t2m/t2m."+YEAR+"-"+MONTH+"*.daily.nc")
        ; print(filist)
        fst = addfile(filist, "r")
        T2MX = short2flt(fst->t2m(:,::-1,:))
        T2M(j*12+i,:,:) = dim_avg_n_Wrap(T2MX, 0)
        copy_VarCoords(T2MX(0,:,:), T2M(j*12+i,:,:))
        delete([/T2MX/])
        print("end of "+YEAR+"-"+MONTH)
    end do
end do
T2M!0 = "time"
printVarSummary(T2M)

t2m_sea = month_to_seasonN(T2M, season)
printVarSummary(t2m_sea)

rc = regCoef_n(ism_onset, t2m_sea, 0, 1)
copy_VarCoords(t2m_sea(:,0,:,:), rc)
printVarSummary(rc)

df = onedtond(rc@nptxy, dimsizes(rc))-2
tval = onedtond(rc@tval, dimsizes(rc))
b = tval
b = 0.5
prob = betainc(df/(df+tval^2), df/2.0, b)
copy_VarCoords(rc, prob)
printVarSummary(prob)

;------------------------------------
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)
;------------------------------------
;=============================================================================

lonstr = 60
lonend = 120
latstr = 10
latend = 50


figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20201116/"
figname = "TP_t2m_reg_on_ISM_onset_pentad"

wks = gsn_open_wks("pdf", figpath+figname) ;
plot = new(2, graphic)
plot_topo = plot
plot_prob = plot

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

res                       = True     ; plot mods desired
res@gsnDraw               = False
res@gsnFrame              = False


res@gsnMaximize           = True     ; Make plot fit the frame
res@cnFillOn              = True     ; turn on color fill
res@cnLinesOn             = False    ; turn of contour lines
res@cnLineLabelsOn        = False
res@cnInfoLabelOn         = False
;res@lbLabelBarOn          = False

res@mpMinLatF            = latstr      ; range to zoom in on
res@mpMaxLatF            = latend
res@mpMinLonF            = lonstr
res@mpMaxLonF            = lonend
; res@mpCenterLonF         = 180.
res@tmXBLabelFontHeightF = 0.02
res@tmYLLabelFontHeightF = 0.02

res@gsnStringFontHeightF = 0.025
; res@gsnRightString       = "K"

res@pmLabelBarOrthogonalPosF = 0.2
res@pmLabelBarWidthF = 0.7
res@lbLabelFontHeightF = 0.015
res@cnFillPalette = "CBR_coldhot"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMaxLevelValF = 0.6
res@cnMinLevelValF = -0.6
res@cnLevelSpacingF = 0.1
res@lbLabelBarOn = False

;----------------------------------------------
; draw prob
;----------------------------------------------
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

 ; res3@cnFillPattern = 17 ;(/17/)
 res3@cnFillColor  =  "black"
 res3@gsnLeftString = ""
 res3@cnFillDotSizeF = 0.003
 res3@cnFillScaleF  = 1.2
 ;----------------------------------------------
 ;----------------------------------------------

do i = 0, 1
    res@gsnLeftString = "t2m: "+season(i)
    plot(i) = gsn_csm_contour_map(wks, rc(i,:,:), res)
    plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(i), plot_topo(i))

    plot_prob(i) = gsn_csm_contour(wks,prob(i,:,:),res3)
    plot_prob(i) = ShadeLtContour(plot_prob(i),0.05,6)
    overlay(plot(i), plot_prob(i))
end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
pres@pmLabelBarWidthF = 0.7

gsn_panel(wks, plot, (/1,2/), pres)
end