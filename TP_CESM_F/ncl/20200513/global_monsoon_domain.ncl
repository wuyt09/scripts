begin

lonstr = 50
lonend = 150
latstr = -10
latend = 50

casename = (/"CTRL","NOQ1","NOTP"/) ;
ncase = dimsizes(casename)
LEV = 850
; z = new((/ncase,12,96,144/), float)
z_range = new((/ncase,96,144/), float)
z_MPI = z_range
; u = z
; v = z

do i = 0, ncase-1
  diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(i)+"_F2000_F19G16_CAM4/pro/"
  zname = "F2000_"+casename(i)+".PREC.mon.0601-3012.nc"
  fz = addfile(diri+zname, "r")
  Z = fz->PREC
  Z = Z*1000*86400
  Z@units = "mm/day"
  z = clmMonTLL(Z)
  z_MJJAS = dim_avg_n_Wrap(z(4:8,:,:), 0)
  z_NDJFM = (dim_sum_n_Wrap(z(0:2,:,:), 0) + dim_sum_n_Wrap(z(10:11,:,:), 0))/5
  z_range(i,:,:) = z_MJJAS - z_NDJFM
  copy_VarCoords(Z(0,:,:), z_range(i,:,:))
  Z_ann = month_to_annual(Z, 1)
  z_ann = dim_avg_n_Wrap(Z_ann, 0)
  printVarSummary(z_ann)
  z_MPI(i,:,:) = z_range(i,:,:)/z_ann
end do
copy_VarCoords(z_range, z_MPI)
printVarSummary(z_range)
; print(z_range&lat)
; exit

z_MPIX = z_MPI
z_rangeX = z_range
do j = 0, 47
    z_rangeX(:,j,:) = -1*z_range(:,j,:)
    z_MPIX(:,j,:) = -1*z_MPI(:,j,:)
end do

; nlat = 96
; nlon = 144
; x = new((/ncase,96,144/), integer)
; x = 0
; do k = 0, ncase-1
;     do j = 0, 47
;         do i = 0, nlon-1
;             z_rangeS = -1*z_range(k,j,i)
;             z_MPIS = -1*z_MPI(k,j,i)
;             inx = z_rangeS.gt.300/365. .and. z_MPIS.gt.0.5
;             if (inx) then
;                 x(k,j,i) = 1
;             end if
;         end do
;     end do
; end do

; do k = 0, ncase-1
;     do j = 48, 95
;         do i = 0, nlon-1
;             inx = z_range(k,j,i).gt.300/365. .and. z_MPI(k,j,i).gt.0.5
;             if (inx) then
;                 x(k,j,i) = 1
;             end if
;         end do
;     end do
; end do
; copy_VarCoords(z_range, x)
; printVarSummary(x)

;-------------------------------------------------------------------

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20200513/"
figname = "global_monsoon_domain"

wks = gsn_open_wks("pdf", figpath+figname)
plot = new(ncase, graphic)
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
res@lbLabelBarOn = False
res@cnInfoLabelOn = False
;res@cnLineDrawOrder = "PostDraw"
; res@cnLabelDrawOrder = "PostDraw"
;res@cnLevelFlag = "LineAndLabel"
;res@cnLineLabelInterval = 1

res@gsnStringFontHeightF = 0.035
res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025

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
;=====================================================

res@cnFillPalette = "CBR_drywet"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -3
res@cnMaxLevelValF = 3
res@cnLevelSpacingF = 0.6

do i = 0, ncase-1
    res@gsnLeftString = casename(i)
    plot(i) = gsn_csm_contour_map(wks, z_MPIX(i,:,:), res)
    plotc(i) = gsn_csm_contour(wks, z_rangeX(i,:,:), cnres)
    overlay(plot(i), plotc(i))
end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
gsn_panel(wks, plot, (/3,1/), pres)


end