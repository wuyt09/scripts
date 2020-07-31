
begin

nlat = 96
nlon = 144

lonstr = 70
lonend = 120
latstr = 5
latend = 28

casename = (/"CTRL","NOTPQ1","NOTPQ1_MAR","NOTPQ1_APR","NOTPQ1_MAY"/)
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
gsn_define_colormap(wks, "example") ;MPL_gist_ncar
nplot = ncase+1
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
res@cnLineLabelsOn = False
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

; res@cnExplicitLineLabelsOn = True
; res@cnLineLabelFontHeightF = 0.015


pens = 22

Colors = (/20,30,40,50,60,70,80,90,100,110/)
text_p = ispan(pens+1, pens+10, 1)
do i = 0, ncase-1
    res@gsnLeftString = casename(i)
    res@cnLineColor = Colors(0)
    ; res@cnLineLabelStrings = "Pentad "+text_p(0)
    plot(i,0) = gsn_csm_contour_map(wks, clim_pentad_z(i,pens,:,:), res)

    delete(res@gsnLeftString)
    do j = 1, 9
        res@cnLineColor = Colors(j)
        ; res@cnLineLabelStrings = "Pentad "+text_p(j)
        plot(i,j) = gsn_csm_contour(wks, clim_pentad_z(i,pens+j,:,:), res)
        overlay(plot(i,0), plot(i,j))
    end do
end do

;==================================================================================
    res_lines                   = True                  ; polyline mods desired
    res_lines@gsLineThicknessF  = 5                     ; line thicker
    

    res_text                    = True                  ; text mods desired
    res_text@txFontHeightF      = 0.018                 ; change text size
    res_text@txJust             = "CenterLeft"          ; text justification

    Texts = "Pentad "+text_p
; ;============================================================
    interval = 0.03
    gx1 = (/0.55,0.6/)
    text_pos1 = 0.62

    gx2 = (/0.75,0.8/)
    text_pos2 = 0.82

    gy1 = (/0.3,0.3/)
    print(gy1)

do i = 0, 4
    res_lines@gsLineColor = Colors(i)
    gsn_polyline_ndc(wks, gx1, gy1-interval*i, res_lines)
    gsn_text_ndc(wks, Texts(i), text_pos1, gy1-interval*i, res_text)


    res_lines@gsLineColor = Colors(i+5)
    gsn_polyline_ndc(wks, gx2, gy1-interval*i, res_lines)
    gsn_text_ndc(wks, Texts(i+5), text_pos2, gy1-interval*i, res_text)
end do
    ; res_lines@gsLineColor = Colors(1)
    ; ;res_lines@gsLineDashPattern = Patterns(1)
    ; dum21 = gsn_add_polyline(wks, plot(0), gx, gy2, res_lines)
    ; dum22 = gsn_add_text(wks, plot(0), Texts(1), text_pos, gy2, res_text)
    ; res_lines@gsLineColor = Colors(2)
    ; ;res_lines@gsLineDashPattern = Patterns(2)
    ; dum31 = gsn_add_polyline(wks, plot(0), gx, gy3, res_lines)
    ; dum32 = gsn_add_text(wks, plot(0), Texts(2), text_pos, gy3, res_text)

pres = True
pres@gsnMaximize = True
; pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = "Z850 1480gpm"
gsn_panel(wks, plot(:,0), (/3,2/), pres)
end