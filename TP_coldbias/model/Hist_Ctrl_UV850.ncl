begin

casename = (/"Ctrl","Hist","Futr","Fubias"/)

ncase = dimsizes(casename)

U = new((/ncase,324,96,144/), float)
V = U

do i = 0, 1 ;ncase-1
    diri = "/home/yangsong3/data-model/wuyt/TPbias_"+casename(i)+"_FAMIP/"
    filiu = casename(i)+".esm.cam.h0.U.1979-2005.nc"
    filiv = casename(i)+".esm.cam.h0.V.1979-2005.nc"

    fu = addfile(diri+filiu, "r")
    fv = addfile(diri+filiv, "r")

    U(i,:,:,:) = fu->U(:,{850},:,:)
    V(i,:,:,:) = fv->V(:,{850},:,:)


delete([/diri,filiu,filiv,fu,fv/])
end do
printVarSummary(U)

u = dim_avg_n_Wrap(U, 1)
v = dim_avg_n_Wrap(V, 1)

uHisteffect = u(1,:,:)
uHisteffect = u(1,:,:) - u(0,:,:)
vHisteffect = v(1,:,:)
vHisteffect = v(1,:,:) - v(0,:,:)
printVarSummary(uHisteffect)
speed = wind_speed(uHisteffect, vHisteffect)

;=========================================================
; plot 
;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/"
figname = "Hist_Ctrl_UV850"

wks = gsn_open_wks("pdf", figpath+figname)
plot = new(2, graphic)
plotv = plot

res = True
res@gsnMaximize = True
res@gsnDraw  = False
res@gsnFrame = False

res@mpFillOn = True
res@mpMinLonF = 40
res@mpMaxLonF = 150
res@mpMinLatF = 10
res@mpMaxLatF = 70

res@cnFillOn = True
res@cnFillPalette = "wind_17lev"

res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
res@cnMinLevelValF       = 0               ; set min contour level
res@cnMaxLevelValF       = 3                 ; set max contour level
res@cnLevelSpacingF      = 0.10               ; set contour spacing
res@tiMainString         = "Sea Surface Temperatures and 1000 mb Winds"

vcres = True
vcres@gsnDraw  = False
vcres@gsnFrame = False
vcres@vcRefMagnitudeF           = 1.0             ; define vector ref mag
vcres@vcRefLengthF              = 0.045           ; define length of vec ref
vcres@vcRefAnnoOrthogonalPosF   = -1.0            ; move ref vector
vcres@vcRefAnnoArrowLineColor   = "black"         ; change ref vector color
vcres@vcRefAnnoArrowUseVecColor = False

vcres@vcGlyphStyle            = "CurlyVector"     ; turn on curly vectors
;vcres@vcLineArrowColor        = "white"           ; change vector color
vcres@vcLineArrowThicknessF   = 2.0               ; change vector thickness
;vcres@vcVectorDrawOrder       = "PostDraw"        ; draw vectors last

plot(0) = gsn_csm_contour_map(wks, speed, res)
plotv(0) = gsn_csm_vector(wks, uHisteffect, vHisteffect, vcres)
overlay(plot(0), plotv(0))



draw(plot(0))
end 