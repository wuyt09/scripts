
begin

nlat = 96
nlon = 144

lonstr = 70
lonend = 100
latstr = 5
latend = 28

casename = (/"CTRL","NOTPQ1","NOTPQ1_MAR","NOTPQ1_MAY"/)
ncase = dimsizes(casename)
clim_pentad_z = new((/ncase,73,nlat,nlon/), float)
; clim_pentad_v = clim_pentad_u
; clim_pentad_z = clim_pentad_u
LEV = 850

do i = 0, ncase-1
    diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(i)+"_F2000_F19G16_CAM4/pro/daily/"
 
    ; uname = systemfunc("ls "+diri+"F2000_"+casename(i)+".plevel.U.*.daily.nc")
    ; vname = systemfunc("ls "+diri+"F2000_"+casename(i)+".plevel.V.*.daily.nc")
    zname = systemfunc("ls "+diri+"F2000_"+casename(i)+".plevel.Z3.1*.daily.nc")
    print(zname)

    fz = addfiles(zname, "r")
    ; fu = addfiles(uname, "r")
    ; fv = addfiles(vname, "r")

    ; U = fu[:]->U(:,{LEV},:,:)
    ; u = reshape(U, (/25,365,nlat,nlon/))
    ; copy_VarCoords(U(0:364,:,:), u(0,:,:,:))

    ; V = fv[:]->V(:,{LEV},:,:)
    ; v = reshape(V, (/25,365,nlat,nlon/))
    ; copy_VarCoords(V(0:364,:,:), v(0,:,:,:))

    Z = fz[:]->Z3(:,{LEV},:,:)
    z = reshape(Z, (/10,365,nlat,nlon/)) ;25
    copy_VarCoords(Z(0:364,:,:), z(0,:,:,:))

    ; pentad_u = u(:,0:72,:,:)
    ; pentad_v = v(:,0:72,:,:)
    pentad_z = z(:,0:72,:,:)

    do p = 0,72  ; a year has 73 penta
      ; pentad_u(:,p,:,:) = dim_sum_n_Wrap(u(:,5*p:5*p+4,:,:), 1)/5.0
      ; pentad_v(:,p,:,:) = dim_sum_n_Wrap(v(:,5*p:5*p+4,:,:), 1)/5.0
      pentad_z(:,p,:,:) = dim_sum_n_Wrap(z(:,5*p:5*p+4,:,:), 1)/5.0
    end do

    ; clim_pentad_U = dim_avg_n_Wrap(pentad_u, 0)
    ; clim_pentad_u(i,:,:,:) = clim_pentad_U ;(18:35,:,:,:)
    ; delete([/U,u,pentad_u,clim_pentad_U/])

    ; clim_pentad_V = dim_avg_n_Wrap(pentad_v, 0)
    ; clim_pentad_v(i,:,:,:) = clim_pentad_V ;(18:35,:,:,:)
    ; delete([/V,v,pentad_v,clim_pentad_V/])

    clim_pentad_Z = dim_avg_n_Wrap(pentad_z, 0)
    clim_pentad_z(i,:,:,:) = clim_pentad_Z ;(18:35,:,:,:)
    delete([/Z,z,pentad_z,clim_pentad_Z/])

end do
printVarSummary(clim_pentad_z)

;-------------------------------------------------------------------

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20200704/"
figname = "H850_monsoon_trough"


wks = gsn_open_wks("pdf", figpath+figname)
gsn_define_colormap(wks, "MPL_rainbow")
nplot = ncase
plot = new((/nplot,10/), graphic)
; plot2 = plot1
; plot3 = plot1
; plot4 = plot1
; plot5 = plot1
; plot6 = plot1

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

res@cnFillOn = False
res@cnLinesOn = True
; res@cnLineLabelsOn = False
res@cnLineLabelFontColor = "red"
; res@cnLineColor = "green4"
res@cnLineThicknessF = 5
res@cnInfoLabelOn = False
;res@cnLineDrawOrder = "PostDraw"
; res@cnLabelDrawOrder = "PostDraw"
;res@cnLevelFlag = "LineAndLabel"
; res@cnLineLabelInterval = 1

res@gsnStringFontHeightF = 0.025
res@tmXBLabelFontHeightF = 0.02
res@tmYLLabelFontHeightF = 0.02


res@cnLevelSelectionMode = "ExplicitLevels"
res@cnLevels = (/1480/)

res@cnExplicitLineLabelsOn = True
res@cnLineLabelFontHeightF = 0.02


pens = 22

Colors = (/20,30,40,50,60,70,80,90,100,110/)
text_p = ispan(pens, pens+10, 1)
do i = 0, ncase-1
    res@gsnLeftString = casename(i)
    res@cnLineColor = Colors(0)
    res@cnLineLabelStrings = "Pentad "+text_p(0)
    plot(i,0) = gsn_csm_contour_map(wks, clim_pentad_z(i,pens,:,:), res)
    do j = 1, 9
        delete(res@gsnLeftString)
        res@cnLineColor = Colors(j)
        res@cnLineLabelStrings = "Pentad "+text_p(j)
        plot(i,j) = gsn_csm_contour(wks, clim_pentad_z(i,pens+j,:,:), res)
        overlay(plot(i,0), plot(i,j))
    end do
end do


; do i = 0, ncase-1
;     res@gsnLeftString = casename(i)
;     ; res@cnLineColor = "grey49"
;     ; plot1(i) = gsn_csm_contour_map(wks, clim_pentad_z(i,24,:,:), res)
;     ; delete(res@gsnLeftString)
;     ; res@cnLineColor = "grey39"
;     ; plot2(i) = gsn_csm_contour(wks, clim_pentad_z(i,25,:,:), res)
;     ; res@cnLineColor = "grey29"
;     ; plot3(i) = gsn_csm_contour(wks, clim_pentad_z(i,26,:,:), res)
;     ; res@cnLineColor = "grey9"
;     ; plot4(i) = gsn_csm_contour(wks, clim_pentad_z(i,27,:,:), res)
;     ; res@cnLineColor = "grey3"
;     ; plot5(i) = gsn_csm_contour(wks, clim_pentad_z(i,28,:,:), res)
;     ; res@cnLineColor = "black"
;     ; plot6(i) = gsn_csm_contour(wks, clim_pentad_z(i,29,:,:), res)

;     res@cnLineColor = "red"
;     plot1(i) = gsn_csm_contour_map(wks, clim_pentad_z(i,pens,:,:), res)
;     delete(res@gsnLeftString)
;     res@cnLineColor = "orange"
;     plot2(i) = gsn_csm_contour(wks, clim_pentad_z(i,pens+1,:,:), res)
;     res@cnLineColor = "yellow"
;     plot3(i) = gsn_csm_contour(wks, clim_pentad_z(i,pens+2,:,:), res)
;     res@cnLineColor = "green"
;     plot4(i) = gsn_csm_contour(wks, clim_pentad_z(i,pens+3,:,:), res)
;     res@cnLineColor = "blue"
;     plot5(i) = gsn_csm_contour(wks, clim_pentad_z(i,pens+4,:,:), res)
;     res@cnLineColor = "purple"
;     plot6(i) = gsn_csm_contour(wks, clim_pentad_z(i,pens+5,:,:), res)
;     overlay(plot1(i), plot2(i))
;     overlay(plot1(i), plot3(i))
;     overlay(plot1(i), plot4(i))
;     overlay(plot1(i), plot5(i))
;     overlay(plot1(i), plot6(i))
; end do


pres = True
pres@gsnMaximize = True
; pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = "Z850 1480gpm"
gsn_panel(wks, plot(:,0), (/ncase/2,2/), pres)
end