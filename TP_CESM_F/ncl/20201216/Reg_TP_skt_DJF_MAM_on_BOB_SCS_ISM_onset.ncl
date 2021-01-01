
begin

lonstr = 75 ; 72 ;75
lonend = 105 ;136 ;105
latstr = 25 ;18 ;25
latend = 40 ;54 ; 40

yrStrt = 1979
yrLast = 2017
mnStrt = 1
mnLast = 12

year = ispan(yrStrt, yrLast, 1)
month = ispan(mnStrt, mnLast, 1)
penta = ispan(1, 73, 1)
nyear = dimsizes(year)
nmonth = dimsizes(month)
npenta = dimsizes(penta)
season = (/"DJF","MAM"/)
monsoon = (/"BOB","SCS","ISM"/)

fx = "/home/yangsong3/wuyt/sysu/scripts/monsoon_onset_pentad/LBQ_onset_pentad_1979_2017.txt"
onset_time = asciiread(fx, (/nyear,4/), "integer")
inx = onset_time(:,1:3)


;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)
printVarSummary(topo)

;===================================================================
diri = "/home/yangsong3/data-observation/ERA5-monthly/surface/"
fili1 = "skt_mon_1x1_1979_2019.nc"

f1 = addfile(diri+fili1, "r")

Time = f1->time
YYYY = cd_calendar(Time, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
var = short2flt(f1->skt(iYYYY,:,:))
printVarSummary(var)

t2m_sea = month_to_seasonN(var, season)
printVarSummary(t2m_sea)

rc = new((/2,3,dimsizes(t2m_sea(0,0,:,0)),dimsizes(t2m_sea(0,0,0,:))/), float)
prob = rc
do i = 0,1
    do j = 0,2
        rc_med = regCoef_n(inx(:,j), t2m_sea(i,:,:,:), 0, 0)
        rc(i,j,:,:) = rc_med

        df = onedtond(rc_med@nptxy, dimsizes(rc_med))-2
        tval = onedtond(rc_med@tval, dimsizes(rc_med))
        b = tval
        b = 0.5
        prob(i,j,:,:) = betainc(df/(df+tval^2), df/2.0, b)
    end do
end do
copy_VarCoords(t2m_sea(:,0,:,:), rc(:,0,:,:))
copy_VarCoords(rc, prob)
printVarSummary(rc)
printVarSummary(prob)


;==============================================================
;plot
;==============================================================
; lonstr = 75
; lonend = 105
; latstr = 25
; latend = 40
lonstr = 60
lonend = 130
latstr = 10
latend = 50


figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20201216/"
figname = "Reg_TP_skt_DJF_MAM_on_BOB_SCS_ISM_onset"

wks = gsn_open_wks("pdf", figpath+figname) ;
plot = new(6, graphic)
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
res@gsnAddCyclic          = False


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
    res@gsnLeftString = season(i)+"_TP_SKT"
    do j = 0,2
        res@gsnRightString = monsoon(j)+"_onset"
        plot(i*3+j) = gsn_csm_contour_map(wks, rc(i,j,:,:), res)
        plot_topo(i*3+j) = gsn_csm_contour(wks, topo, res_topo)
        overlay(plot(i*3+j), plot_topo(i*3+j))

        plot_prob(i*3+j) = gsn_csm_contour(wks,prob(i,j,:,:),res3)
        plot_prob(i*3+j) = ShadeLtContour(plot_prob(i*3+j),0.05,6)
        overlay(plot(i*3+j), plot_prob(i*3+j))
    end do
end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
pres@pmLabelBarWidthF = 0.7
pres@lbLabelFontHeightF = 0.01
; pres@gsnPanelMainString = "BOB onset"

gsn_panel(wks, plot, (/2,3/), pres)
end