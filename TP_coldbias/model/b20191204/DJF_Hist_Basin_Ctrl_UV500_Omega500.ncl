begin

casename = (/"Ctrl","Hist","Hist_IO","Hist_PO","Hist_AO","Hist_PIO","Hist_PAO","Hist_IAO"/)

ncase = dimsizes(casename)

;U = new((/ncase,324,96,144/), float)
U = new((/ncase,27,96,144/), float)
V = U
W = U
LEV = 500
season = "DJF"

do i = 0, ncase-1
    diri = "/home/yangsong3/data-model/wuyt/TPbias_"+casename(i)+"_FAMIP/"
    filiu = casename(i)+".esm.cam.h0.U.1979-2005.nc"
    filiv = casename(i)+".esm.cam.h0.V.1979-2005.nc"
    filiw = casename(i)+".esm.cam.h0.OMEGA.1979-2005.nc"

    fu = addfile(diri+filiu, "r")
    fv = addfile(diri+filiv, "r")
    fw = addfile(diri+filiw, "r")

    UX = fu->U(:,{LEV},:,:)
    VX = fv->V(:,{LEV},:,:)
    WX = fw->OMEGA(:,{LEV},:,:)
    U(i,:,:,:) = month_to_season(UX, season)
    V(i,:,:,:) = month_to_season(VX, season)
    W(i,:,:,:) = month_to_season(WX, season)
delete([/diri,filiu,filiv,filiw,fu,fv,fw,\
       UX,VX,WX/])
end do

printVarSummary(U)

u = dim_avg_n_Wrap(U, 1)
v = dim_avg_n_Wrap(V, 1)
w = dim_avg_n_Wrap(W, 1)
;print(prec(0,{20:40},{110:120}))
;print(prec(0,{20:40},{75:105}))

udiff = u(0:ncase-2,:,:)
vdiff = v(0:ncase-2,:,:)
wdiff = w(0:ncase-2,:,:)

do i = 0, ncase-2
    udiff(i,:,:) = u(i+1,:,:) - u(0,:,:)
    vdiff(i,:,:) = v(i+1,:,:) - v(0,:,:)
    wdiff(i,:,:) = w(i+1,:,:) - w(0,:,:)
end do
printVarSummary(udiff)

;=========================================================
; outline of tp 
;=========================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo

;=========================================================
; plot 
;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/"
figname = season+"_Hist_Basin_Ctrl_UV500_Omega500"

; lonstr = 40
; lonend = 150
; latstr = 10
; latend = 70
lonstr = 0
lonend = 360
latstr = -30
latend = 90

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
res@cnFillPalette = "BlueWhiteOrangeRed"
res@lbLabelBarOn = False
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False


res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
; res@cnMinLevelValF       = -0.1               ; set min contour level
; res@cnMaxLevelValF       = 0.1                 ; set max contour level
; res@cnLevelSpacingF      = 0.01              ; set contour spacing
res@cnMinLevelValF       = -0.05               ; set min contour level
res@cnMaxLevelValF       = 0.05                 ; set max contour level
res@cnLevelSpacingF      = 0.005              ; set contour spacing
res@gsnRightString       = "Pa/s"
res@gsnStringFontHeightF = 0.03
res@tmXBLabelFontHeightF = 0.02
res@tmYLLabelFontHeightF = 0.02

vcres = True
vcres@gsnDraw  = False
vcres@gsnFrame = False
vcres@vcRefMagnitudeF           = 2.0             ; define vector ref mag
vcres@vcRefLengthF              = 0.045           ; define length of vec ref
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

leftsti = (/"Total","IndianOnly","PacificOnly","AtlanticOnly",\
           "Pacific&Indian","Pacific&Atlantic","Indian&Atlantic"/)
do i = 0, ncase-2
    res@gsnLeftString = leftsti(i)
    plot(i) = gsn_csm_contour_map(wks, wdiff(i,:,:), res)
    plotv(i) = gsn_csm_vector(wks, udiff(i,:,:), vdiff(i,:,:), vcres)
    overlay(plot(i), plotv(i))
    plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(i), plot_topo(i))
end do


pres = True
pres@gsnDraw = True
pres@gsnFrame = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = season+" mean: UV500 & Omega500"

gsn_panel(wks, plot, (/4,2/), pres)
end 