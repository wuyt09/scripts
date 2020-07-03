begin

casename = (/"Ctrl","Hist","Hist_IO","Hist_PO","Hist_AO","Hist_PIO","Hist_PAO","Hist_IAO"/)

ncase = dimsizes(casename)

;U = new((/ncase,324,96,144/), float)
U = new((/ncase,27,17,96,144/), float)
V = U
W = U
;LEV = 500
season = "DJF"

do i = 0, ncase-1
    diri = "/home/yangsong3/data-model/wuyt/TPbias_"+casename(i)+"_FAMIP/"
    ;filiu = casename(i)+".esm.cam.h0.U.1979-2005.nc"
    filiv = casename(i)+".esm.cam.h0.V.1979-2005.nc"
    filiw = casename(i)+".esm.cam.h0.OMEGA.1979-2005.nc"

    ;fu = addfile(diri+filiu, "r")
    fv = addfile(diri+filiv, "r")
    fw = addfile(diri+filiw, "r")

    ;UX = fu->U(:,{LEV},:,:)
    VX = fv->V
    WX = fw->OMEGA
    ;U(i,:,:,:,:) = month_to_annual(UX, 1)
    V(i,:,:,:,:) = month_to_season(VX, season)
    W(i,:,:,:,:) = month_to_season(WX, season)
delete([/diri,filiv,filiw,fv,fw,VX,WX/])
end do

printVarSummary(V)

;u = dim_avg_n_Wrap(U, 1)
v = dim_avg_n_Wrap(V, 1)
w = dim_avg_n_Wrap(W, 1)
;print(prec(0,{20:40},{110:120}))
;print(prec(0,{20:40},{75:105}))

;udiff = u(0:ncase-2,:,:,:)
vdiff = v(0:ncase-2,:,:,:)
wdiff = w(0:ncase-2,:,:,:)

do i = 0, ncase-2
    ;udiff(i,:,:,:) = u(i+1,:,:,:) - u(0,:,:,:)
    vdiff(i,:,:,:) = v(i+1,:,:,:) - v(0,:,:,:)
    wdiff(i,:,:,:) = w(i+1,:,:,:) - w(0,:,:,:)
end do
printVarSummary(vdiff)

wgt = vdiff(0,0,0,{75:105})
wgt = 1
vcell = dim_avg_wgt_n_Wrap(vdiff(:,:,:,{75:105}), wgt, 0, 3)
wcell = dim_avg_wgt_n_Wrap(wdiff(:,:,:,{75:105}), wgt, 0, 3)
wcell = wcell*100.
;=========================================================
; outline of tp 
;=========================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo

;=========================================================
; plot 
;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/"
figname = season+"_Hist_Basin_Ctrl_vertical_cell_meridional"

; lonstr = 40
; lonend = 150
; latstr = 10
; latend = 70
; lonstr = 0
; lonend = 360
latstr = -30
latend = 80

wks = gsn_open_wks("pdf", figpath+figname)
plot = new(ncase, graphic)
; plotv = plot
; plot_topo = plot

;  res_topo            = True
;  res_topo@cnFillOn   = False
;  res_topo@cnLinesOn  = True
;  res_topo@gsnFrame   = False
;  res_topo@gsnDraw    = False
;  res_topo@cnLineLabelsOn         = False
;  res_topo@cnLineThicknessF       = 3
;  res_topo@cnLineColor            = "grey"
;  res_topo@cnInfoLabelOn          = False
;  res_topo@cnLevelSelectionMode   ="ExplicitLevels"
;  res_topo@cnLevels               = (/1500/)
;  res_topo@gsnLeftString = ""
;  res_topo@gsnRightString = ""


res = True
res@gsnMaximize = True
res@gsnDraw  = False
res@gsnFrame = False

; res@mpFillOn = True
; res@mpMinLonF = lonstr
; res@mpMaxLonF = lonend
; res@mpMinLatF = latstr
; res@mpMaxLatF = latend

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
res@cnMinLevelValF       = -0.8               ; set min contour level
res@cnMaxLevelValF       = 0.8                 ; set max contour level
res@cnLevelSpacingF      = 0.1              ; set contour spacing
res@gsnRightString       = "Pa/s"
res@gsnStringFontHeightF = 0.04
res@tmXBLabelFontHeightF = 0.03
res@tmYLLabelFontHeightF = 0.03

; vcres = True
; vcres@gsnDraw  = False
; vcres@gsnFrame = False
res@vcRefMagnitudeF           = 1.0             ; define vector ref mag
res@vcRefLengthF              = 0.045           ; define length of vec ref
res@vcMinDistanceF            = 0.025
res@vcRefAnnoOrthogonalPosF   = -1.0            ; move ref vector
res@vcRefAnnoArrowLineColor   = "black"         ; change ref vector color
;vcres@vcRefAnnoArrowUseVecColor = False
res@vcRefAnnoString2On = False
res@vcRefAnnoFontHeightF = 0.02

res@vcGlyphStyle            = "LineArrow"     ; turn on curly vectors CurlyVector
;vcres@vcLineArrowColor        = "white"           ; change vector color
res@vcLineArrowThicknessF   = 3.0               ; change vector thickness
;vcres@vcVectorDrawOrder       = "PostDraw"        ; draw vectors last

leftsti = (/"Total","IndianOnly","PacificOnly","AtlanticOnly",\
           "Pacific&Indian","Pacific&Atlantic","Indian&Atlantic"/)
; do i = 0, ncase-2
;     res@gsnLeftString = leftsti(i)
;     plot(i) = gsn_csm_contour_map(wks, wdiff(i,:,:), res)
;     plotv(i) = gsn_csm_vector(wks, udiff(i,:,:), vdiff(i,:,:), vcres)
;     overlay(plot(i), plotv(i))
;     plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
;     overlay(plot(i), plot_topo(i))
; end do

do i = 0, ncase-2
    res@gsnLeftString = leftsti(i)
    plot(i) = gsn_csm_pres_hgt_vector(wks, wcell(i,{1000:100},{latstr:latend}), vcell(i,{1000:100},{latstr:latend}), wcell(i,{1000:100},{latstr:latend}), res)
    ;plotv(i) = gsn_csm_vector(wks, udiff(i,:,:), vdiff(i,:,:), vcres)
    ;overlay(plot(i), plotv(i))
    ; plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
    ; overlay(plot(i), plot_topo(i))
end do


pres = True
pres@gsnDraw = True
pres@gsnFrame = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = season+" mean: Vertical Cell (75-105E)"

gsn_panel(wks, plot, (/2,4/), pres)
end 