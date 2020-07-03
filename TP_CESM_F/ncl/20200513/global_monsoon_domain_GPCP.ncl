begin

lonstr = 50
lonend = 150
latstr = -10
latend = 50

z_range = new((/72,144/), float)
z_MPI = z_range

diri = "/home/yangsong3/data-observation/GPCP/"
zname = "precip.mon.mean.nc"
fz = addfile(diri+zname, "r")
TIME = fz->time
YYYY = cd_calendar(TIME, -1)/100
iYYYY = ind(YYYY.ge.1979 .and. YYYY.le.2015)
Z = fz->precip(iYYYY,:,:)

z = clmMonTLL(Z)
z_MJJAS = dim_avg_n_Wrap(z(4:8,:,:), 0)
z_NDJFM = (dim_sum_n_Wrap(z(0:2,:,:), 0) + dim_sum_n_Wrap(z(10:11,:,:), 0))/5
z_range = z_MJJAS - z_NDJFM
copy_VarCoords(Z(0,:,:), z_range)
Z_ann = month_to_annual(Z, 1)
z_ann = dim_avg_n_Wrap(Z_ann, 0)
printVarSummary(z_ann)
z_MPI = z_range/z_ann
copy_VarCoords(z_range, z_MPI)
printVarSummary(z_range)
; print(z_range&lat)
; exit

; nlat = 72
; nlon = 144
; x = new((/72,144/), integer)
; x = 0
; do j = 0, 35
;     do i = 0, nlon-1
;         z_rangeS = -1*z_range(j,i)
;         z_MPIS = -1*z_MPI(j,i)
;         inx = z_rangeS.gt.300/365. .and. z_MPIS.gt.0.5
;         if (inx) then
;           x(j,i) = 1
;         end if
;     end do
; end do

; do j = 36, 71
;     do i = 0, nlon-1
;         inx = z_range(j,i).gt.300/365. .and. z_MPI(j,i).gt.0.5
;         if (inx) then
;           x(j,i) = 1
;         end if
;     end do
; end do
; copy_VarCoords(z_range, x)
; printVarSummary(x)
; print(x({10},{90:120}))

; ondx = ndtooned(x)
; indx = ind(ondx.eq.1)
; xnd = ind_resolve(indx, dimsizes(x))
; print(xnd)

z_rangeX = z_range
z_MPIX = z_MPI
do j = 0,35
  z_rangeX(j,:) = -1*z_range(j,:)
  z_MPIX(j,:) = -1*z_MPI(j,:)
end do
;z_rangeX = z_rangeX*365
printVarSummary(z_rangeX)
printVarSummary(z_MPIX)

;-------------------------------------------------------------------

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20200513/"
figname = "global_monsoon_domain_GPCP"

wks = gsn_open_wks("pdf", figpath+figname)
plot = new(1, graphic)
plotc = plot
;----------------------------------------------------

res = True
res@gsnDraw  = False
res@gsnFrame = False

res@mpOutlineOn = True
res@mpGeophysicalLineThicknessF = 2
res@mpGeophysicalLineColor = "grey25" ;"red" 
res@mpFillDrawOrder = "PreDraw"
res@mpFillOn = True

res@mpMinLatF = -50
res@mpMaxLatF = 50
res@mpCenterLonF = 180
; res@mpMinLonF = lonstr
; res@mpMaxLonF = lonend

res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@lbLabelBarOn = True
res@cnInfoLabelOn = False
;res@cnLineDrawOrder = "PostDraw"
; res@cnLabelDrawOrder = "PostDraw"
;res@cnLevelFlag = "LineAndLabel"
;res@cnLineLabelInterval = 1

res@gsnStringFontHeightF = 0.035
res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025

res@gsnLeftString = "GPCP"

;=====================================================
    cnres = True
    cnres@gsnDraw = False
    cnres@gsnFrame = False

    cnres@cnFillOn = False
    cnres@cnLinesOn = True
    cnres@cnLineLabelsOn = False
    ; cnres@cnLineLabelFontColor = "red"
    cnres@cnLineColor = "red" ;goldenrod2 ;dodgerblue2
    cnres@cnLineThicknessF = 5
    cnres@cnInfoLabelOn = True
    ;cnres@cnLineDrawOrder = "PostDraw"
    cnres@cnLabelDrawOrder = "PostDraw"
    ;cnres@cnLevelFlag = "LineAndLabel"
    cnres@cnLineLabelInterval = 1
    cnres@gsnLeftString = ""
    cnres@gsnRightString = ""
    cnres@cnLevelSelectionMode = "ExplicitLevels"
    cnres@cnLevels = (/2.5/)
;=====================================================

res@cnFillPalette = "CBR_drywet"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -3
res@cnMaxLevelValF = 3
res@cnLevelSpacingF = 0.6
plot(0) = gsn_csm_contour_map(wks, z_MPIX, res)
plotc(0) = gsn_csm_contour(wks, z_rangeX, cnres)
overlay(plot(0), plotc(0))
draw(plot)


end