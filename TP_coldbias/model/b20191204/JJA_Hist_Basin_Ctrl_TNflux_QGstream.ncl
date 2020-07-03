begin

casename = (/"Ctrl","Hist","Hist_IO","Hist_PO","Hist_AO","Hist_PIO","Hist_PAO","Hist_IAO"/)

ncase = dimsizes(casename)

;U = new((/ncase,324,96,144/), float)
U = new((/ncase,27,96,144/), float)
V = U
S = U
LEV = 200
season = "JJA"

diri = "/home/yangsong3/wuyt/sysu/scripts/TP_coldbias/model/TNflux/"
do i = 0, ncase-1
    filiu = "20191123_TN2001-Fx.monthly."+casename(i)+".nc"
    filiv = "20191123_TN2001-Fy.monthly."+casename(i)+".nc"
    filis = "20191123_psidev.monthly."+casename(i)+".nc"

    fu = addfile(diri+filiu, "r")
    fv = addfile(diri+filiv, "r")
    fs = addfile(diri+filis, "r")

    UXX = fu->Fx(:,{LEV},:,:)
    VXX = fv->Fy(:,{LEV},:,:)
    SXX = fs->psidev(:,{LEV},:,:)
    UX = tofloat(UXX)
    copy_VarCoords(UXX, UX)
    VX = tofloat(VXX)
    copy_VarCoords(VXX, VX)
    SX = tofloat(SXX)
    copy_VarCoords(SXX, SX)

    U(i,:,:,:) = month_to_season(UX, season)
    V(i,:,:,:) = month_to_season(VX, season)
    S(i,:,:,:) = month_to_season(SX, season)

delete([/filiu,filiv,fu,fv,fs,UX,VX,SX/])
end do

printVarSummary(U)

u = dim_avg_n_Wrap(U, 1)
v = dim_avg_n_Wrap(V, 1)
s = dim_avg_n_Wrap(S, 1)
;print(prec(0,{20:40},{110:120}))
;print(prec(0,{20:40},{75:105}))

udiff = u(0:ncase-2,:,:)
vdiff = v(0:ncase-2,:,:)
sdiff = s(0:ncase-2,:,:)

do i = 0, ncase-2
    udiff(i,:,:) = u(i+1,:,:) - u(0,:,:)
    vdiff(i,:,:) = v(i+1,:,:) - v(0,:,:)
    sdiff(i,:,:) = s(i+1,:,:) - s(0,:,:)
end do
;precdiff(0,:,:) = prec(0,:,:)
printVarSummary(udiff)
; print(sdiff(0,{25:40},{110:120}))
; print(udiff(0,{25:40},{110:120}))
; print(vdiff(0,{25:40},{110:120}))
; exit
;=========================================================
; outline of tp 
;=========================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo

;=========================================================
; plot 
;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/"
figname = season+"_Hist_Basin_Ctrl_TNflux_QGstream"

; lonstr = 40
; lonend = 150
; latstr = 10
; latend = 70
lonstr = 0
lonend = 360
latstr = 10
latend = 80

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
res@lbLabelBarOn = False
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False


res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
; res@cnMinLevelValF       = -0.1               ; set min contour level
; res@cnMaxLevelValF       = 0.1                 ; set max contour level
; res@cnLevelSpacingF      = 0.01              ; set contour spacing
res@cnMinLevelValF       = -10               ; set min contour level
res@cnMaxLevelValF       = 10                 ; set max contour level
res@cnLevelSpacingF      = 4              ; set contour spacing
res@gsnRightString       = "m~S~2~N~/s~S~2~N~"
res@gsnStringFontHeightF = 0.03
res@tmXBLabelFontHeightF = 0.02
res@tmYLLabelFontHeightF = 0.02

vcres = True
vcres@gsnDraw  = False
vcres@gsnFrame = False
vcres@vcRefMagnitudeF           = 1             ; define vector ref mag
vcres@vcRefLengthF              = 0.045           ; define length of vec ref
vcres@vcMinDistanceF            = 0.02
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

leftsti = (/"Total","IndianOnly","PacificOnly","AtlanticOnly",\
           "Pacific&Indian","Pacific&Atlantic","Indian&Atlantic"/)

vcres@vcRefMagnitudeF           = 1             ; define vector ref mag
do i = 0, ncase-2
    res@gsnLeftString = leftsti(i)
    ;plot(i) = gsn_csm_contour_map(wks, sdiff(i,:,:), res)
    plot(i) = gsn_csm_map_ce(wks, res)
    plotv(i) = gsn_csm_vector(wks, udiff(i,:,:), vdiff(i,:,:), vcres)
    overlay(plot(i), plotv(i))
    plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(i), plot_topo(i))
end do

; vcres@vcRefMagnitudeF           = 2             ; define vector ref mag
; do i = 0, ncase-1
;     res@gsnLeftString = casename(i)
;     ;plot(i) = gsn_csm_contour_map(wks, s(i,:,:), res)
;     plot(i) = gsn_csm_map_ce(wks, res)
;     plotv(i) = gsn_csm_vector(wks, u(i,:,:), v(i,:,:), vcres)
;     overlay(plot(i), plotv(i))
;     plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
;     overlay(plot(i), plot_topo(i))
; end do


pres = True
pres@gsnDraw = True
pres@gsnFrame = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = season+" mean: WaveFlux";& StreamFunction

gsn_panel(wks, plot, (/4,2/), pres)
end 