
begin

nlat = 96
nlon = 144

lonstr = 100
lonend = 150
latstr = 5
latend = 50

casename = (/"CTRL","NOTPQ1","NOTPQ1_MAR","NOTPQ1_APR","NOTPQ1_MAY"/)
ncase = dimsizes(casename)
clim_pentad_z = new((/ncase,73,nlat,nlon/), float)
LEV = 500

do i = 0, ncase-1
    diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(i)+"_F2000_F19G16_CAM4/pro/daily/"
 
    zname = systemfunc("ls "+diri+"F2000_"+casename(i)+".plevel.Z3.1*.daily.nc")
    print(zname)

    fz = addfiles(zname, "r")

    Z = fz[:]->Z3(:,{LEV},:,:)
    z = reshape(Z, (/10,365,nlat,nlon/)) ;25
    copy_VarCoords(Z(0:364,:,:), z(0,:,:,:))

    pentad_z = z(:,0:72,:,:)

    do p = 0,72  ; a year has 73 penta
      pentad_z(:,p,:,:) = dim_sum_n_Wrap(z(:,5*p:5*p+4,:,:), 1)/5.0
    end do

    clim_pentad_Z = dim_avg_n_Wrap(pentad_z, 0)
    clim_pentad_z(i,:,:,:) = clim_pentad_Z ;(18:35,:,:,:)
    delete([/Z,z,pentad_z,clim_pentad_Z/])

end do
printVarSummary(clim_pentad_z)

;-------------------------------------------------------------------

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20200704/"
figname = "H500_WPSH"


wks = gsn_open_wks("pdf", figpath+figname)
nplot = ncase
plot = new((/nplot,12/), graphic)

res = True
res@gsnDraw  = False
res@gsnFrame = False

res@mpOutlineOn = True
res@mpGeophysicalLineThicknessF = 2
res@mpGeophysicalLineColor = "grey25" ;"red" 
res@mpFillDrawOrder = "PreDraw"
res@mpFillOn = True

res@mpMinLatF = latstr
res@mpMaxLatF = latend
res@mpMinLonF = lonstr
res@mpMaxLonF = lonend

res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnLineLabelFontColor = "red"
; res@cnLineColor = "green4"
; res@cnLineThicknessF = 5
res@cnInfoLabelOn = False
;res@cnLineDrawOrder = "PostDraw"
; res@cnLabelDrawOrder = "PostDraw"
;res@cnLevelFlag = "LineAndLabel"
; res@cnLineLabelInterval = 1

res@gsnStringFontHeightF = 0.025
res@tmXBLabelFontHeightF = 0.02
res@tmYLLabelFontHeightF = 0.02


res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = 5860
res@cnMaxLevelValF = 5920
res@cnLevelSpacingF = 5
res@cnFillPalette = "BlueWhiteOrangeRed"

res@cnExplicitLineLabelsOn = True
res@cnLineLabelFontHeightF = 0.02


pens = 24

Colors = (/20,30,40,50,60,70,80,90,100,110/)
text_p = ispan(pens+1, pens+12, 1)
; do i = 0, ncase-1
;     res@gsnLeftString = casename(i)
;     res@cnLineColor = Colors(0)
;     res@cnLineLabelStrings = "Pentad "+text_p(0)
;     plot(i,0) = gsn_csm_contour_map(wks, clim_pentad_z(i,pens,:,:), res)
;     do j = 1, 11
;         delete(res@gsnLeftString)
;         res@cnLineColor = Colors(j)
;         res@cnLineLabelStrings = "Pentad "+text_p(j)
;         plot(i,j) = gsn_csm_contour(wks, clim_pentad_z(i,pens+j,:,:), res)
;         overlay(plot(i,0), plot(i,j))
;     end do
; end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True

do i = 0, ncase-1
    do j = 0, 11
        res@gsnLeftString = "Pentad "+text_p(j)
        plot(i,j) = gsn_csm_contour_map(wks, clim_pentad_z(i,pens+j,:,:), res)
    end do
    pres@gsnPanelMainString = casename(i)+" Z500"
    gsn_panel(wks, plot(i,:), (/3,4/), pres)
end do

end