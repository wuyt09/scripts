begin

lonstr = 75
lonend = 105
latstr = 25
latend = 40
yrStrt = 198102
yrLast = 200601

casename = (/"Ctrl_FAMIP","Hist_FAMIP"/)
ncase = dimsizes(casename)

season = "DJF"
u = new((/ncase,96,144/), float)
v = u
speed = u
LEV = 200

do i = 0, ncase-1
    diri = "/home/yangsong3/data-model/wuyt/TPbias/TPbias_"+casename(i)+"/a20191206/"
    filiu = casename(i)+".cam.h0.U.1979-2005.nc"
    filiv = casename(i)+".cam.h0.V.1979-2005.nc"

    fu = addfile(diri+filiu, "r")
    fv = addfile(diri+filiv, "r")

    TIME = fu->time
    YYYY = cd_calendar(TIME, -1)
    iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

    UX = fu->U(iYYYY,{LEV},:,:)
    VX = fv->V(iYYYY,{LEV},:,:)
    ;SPEEDX = wind_speed(UX, VX)
    SPEEDX = UX
    U = month_to_season(UX, season)
    V = month_to_season(VX, season)
    u(i,:,:) = dim_avg_n_Wrap(U, 0)
    v(i,:,:) = dim_avg_n_Wrap(V, 0)
    SPEED = month_to_season(SPEEDX, season)
    speed(i,:,:) = dim_avg_n_Wrap(SPEED, 0)

delete([/diri,filiu,filiv,fu,fv,UX,VX,SPEEDX/])
end do

printVarSummary(U)

prec = speed

udiff = u(0:ncase-2,:,:)
vdiff = v(0:ncase-2,:,:)
precdiff = prec(0:ncase-2,:,:)

do i = 0, ncase-2
    udiff(i,:,:) = u(i+1,:,:) - u(0,:,:)
    vdiff(i,:,:) = v(i+1,:,:) - v(0,:,:)
    precdiff(i,:,:) = prec(i+1,:,:) - prec(0,:,:)
end do
;precdiff(0,:,:) = prec(0,:,:)
printVarSummary(udiff)

;=========================================================
; outline of tp 
;=========================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo

;=========================================================
; plot 
;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/a20200207/"
figname = season+"_westerly_jet"

lonstr = 40
lonend = 150
latstr = 10
latend = 70
; lonstr = 0
; lonend = 360
; latstr = -30
; latend = 90

wks = gsn_open_wks("pdf", figpath+figname)
plot = new(ncase, graphic)
plotv = plot
plot_topo = plot

 res_topo            = True
 res_topo@cnFillOn   = False
 res_topo@cnLinesOn  = True
 res_topo@gsnFrame   = False
 res_topo@gsnDraw    = False
 res_topo@cnLineLabelsOn         = False
 res_topo@cnLineThicknessF       = 3
 res_topo@cnLineColor            = "grey"
 res_topo@cnInfoLabelOn          = False
 res_topo@cnLevelSelectionMode   ="ExplicitLevels"
 res_topo@cnLevels               = (/1500/)
 res_topo@gsnLeftString = ""
 res_topo@gsnRightString = ""

res = True
res@gsnMaximize = True
res@gsnDraw  = False
res@gsnFrame = False

res@mpFillOn = True
res@mpMinLonF = lonstr
res@mpMaxLonF = lonend
res@mpMinLatF = latstr
res@mpMaxLatF = latend

res@cnFillOn = True
res@cnFillPalette = "BlueWhiteOrangeRed";WhiteBlueGreenYellowRed
;res@lbLabelBarOn = False
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False


res@gsnRightString       = "m/s"
res@gsnStringFontHeightF = 0.03
res@tmXBLabelFontHeightF = 0.02
res@tmYLLabelFontHeightF = 0.02

vcres = True
vcres@gsnDraw  = False
vcres@gsnFrame = False
vcres@vcMinDistanceF            = 0.025
vcres@vcRefAnnoOrthogonalPosF   = -1.0            ; move ref vector
vcres@vcRefAnnoArrowLineColor   = "black"         ; change ref vector color
;vcres@vcRefAnnoArrowUseVecColor = False
vcres@vcRefAnnoString2On = False
vcres@vcRefAnnoFontHeightF = 0.02

vcres@gsnLeftString = ""
vcres@gsnRightString = ""

vcres@vcGlyphStyle            = "LineArrow"     ; turn on curly vectors CurlyVector
;vcres@vcLineArrowColor        = "white"           ; change vector color
vcres@vcLineArrowThicknessF   = 3.0               ; change vector thickness
;vcres@vcVectorDrawOrder       = "PostDraw"        ; draw vectors last

res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
res@cnMinLevelValF       = 0               ; set min contour level
res@cnMaxLevelValF       = 60                ; set max contour level
res@cnLevelSpacingF      = 5.                ; set contour spacing
vcres@vcRefMagnitudeF           = 40.0             ; define vector ref mag
vcres@vcRefLengthF              = 0.045           ; define length of vec ref
res@gsnLeftString = "Ctrl"
plot(0) = gsn_csm_contour_map(wks, prec(0,:,:), res)
plotv(0) = gsn_csm_vector(wks, u(0,:,:), v(0,:,:), vcres)
overlay(plot(0), plotv(0))
plot_topo(0) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(0), plot_topo(0))

res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
res@cnMinLevelValF       = -5               ; set min contour level
res@cnMaxLevelValF       = 5                 ; set max contour level
res@cnLevelSpacingF      = 0.5              ; set contour spacing
vcres@vcRefMagnitudeF           = 5.0             ; define vector ref mag
vcres@vcRefLengthF              = 0.045           ; define length of vec ref
res@gsnLeftString = "Diff"
plot(1) = gsn_csm_contour_map(wks, precdiff(0,:,:), res)
plotv(1) = gsn_csm_vector(wks, udiff(0,:,:), vdiff(0,:,:), vcres)
overlay(plot(1), plotv(1))
plot_topo(1) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(1), plot_topo(1))

pres = True
pres@gsnDraw = True
pres@gsnFrame = True
pres@gsnMaximize = True
;pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = season+" mean: UV200 & Jetstream"

gsn_panel(wks, plot, (/1,2/), pres)
end 