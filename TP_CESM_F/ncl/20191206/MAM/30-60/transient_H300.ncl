
begin
;==================================================================
ca = 60
cb = 30
index = asciiread("./high_point_index_"+cb+"-"+ca+"day.txt", -1, "integer")
print(index)

;==================================================================

LEV = 500

;casename = (/"CTRL","NOTP","NOQ1","NOTP_NOQ1"/)
casename = (/"CTRL","NOQ1"/)
ncase = dimsizes(casename)
H = new((/ncase,25,92,96,144/), float)
year = ispan(6, 30, 1)
do i = 0, ncase-1
    diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(i)+"_F2000_F19G16_CAM4/pro/daily/"
    do j = 0, dimsizes(year)-1
        yr = sprinti("%0.2i", year(j))
        fili = "F2000_"+casename(i)+".plevel.Z3."+yr+".daily.nc"
        f = addfile(diri+fili, "r")
        TIME = f->time
        YYYYMM = cd_calendar(TIME, -1)
        MM = mod(YYYYMM, 100)
        iMM = ind(MM.ge.3 .and. MM.le.5)
        H(i,j,:,:,:) = f->Z3(iMM,{LEV},:,:)
    end do
end do
printVarSummary(H)

H_clm = dim_avg_n_Wrap(H, 1)
h_clm = conform_dims(dimsizes(H), H_clm, (/0,2,3,4/))
H_ano = H
H_ano = H - h_clm
h = reshape(H_ano, (/ncase,25*92,96,144/))
copy_VarCoords(H(0,0,0,:,:), h(0,0,:,:))
h!0 = "case"
h&case = casename
printVarSummary(h)
h_thermal = h(1,:,:,:)
h_thermal = h(0,:,:,:) - h(1,:,:,:)


; H0 = reshape(H, (/ncase,25*92,96,144/))
; copy_VarCoords(H(0,0,0,:,:), H0(0,0,:,:))

; fca = 1./ca
; fcb = 1./cb
; h_bp = bw_bandpass_filter(H0, fca, fcb, True, 1)
; printVarSummary(h_bp)

; h_thermal = h_bp(1,:,:,:)
; h_thermal = h_bp(0,:,:,:) - h_bp(1,:,:,:)
; copy_VarCoords(H0(0,0,:,:), h_thermal(0,:,:))

HX = new((/7,dimsizes(index),96,144/), float)

days = (/10,20,30/)
HX(0,:,:,:) = h_thermal(index-days(2),:,:)
HX(1,:,:,:) = h_thermal(index-days(1),:,:)
HX(2,:,:,:) = h_thermal(index-days(0),:,:)
HX(3,:,:,:) = h_thermal(index,:,:)
HX(4,:,:,:) = h_thermal(index+days(0),:,:)
HX(5,:,:,:) = h_thermal(index+days(1),:,:)
HX(6,:,:,:) = h_thermal(index+days(2),:,:)

hx = dim_avg_n_Wrap(HX, 1)
printVarSummary(hx)

;==================================================================

tpmaskpath = "/home/yangsong3/wuyt/sysu/scripts/TP_CESM_F/model_need_input_data/"
tpmaskname = "TPmask_1000m_C20190410.nc"
ftp = addfile(tpmaskpath+tpmaskname, "r")

tpmask = ftp->tp
printVarSummary(tpmask)

;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20191206/MAM/30-60/"
figname = "transient_H"+LEV+"_"+cb+"-"+ca+"day"

lonstr = 0 ;20
lonend = 360 ;220
latstr = 10
latend = 80

wks = gsn_open_wks("pdf", figpath+figname)
nplot = 7
plot = new(nplot, graphic)
plotv = plot
plot_topo = plot

 res_topo            = True
 res_topo@cnFillOn   = False
 res_topo@cnLinesOn  = True
 res_topo@gsnFrame   = False
 res_topo@gsnDraw    = False
 res_topo@cnLineLabelsOn         = False
 res_topo@cnLineThicknessF       = 3
 res_topo@cnLineColor            = "blue"
 res_topo@cnInfoLabelOn          = False
 res_topo@cnLevelSelectionMode   ="ExplicitLevels"
 res_topo@cnLevels               = (/0.9/)
 res_topo@gsnLeftString = ""
 res_topo@gsnRightString = ""


res = True
res@gsnDraw  = False
res@gsnFrame = False

res@mpOutlineOn = True
res@mpFillOn = True
res@mpGeophysicalLineThicknessF = 2
res@mpGeophysicalLineColor = "grey45" ;"grey25"
;res@mpFillDrawOrder = "PreDraw"

res@mpMinLatF = latstr
res@mpMaxLatF = latend
res@mpMinLonF = lonstr
res@mpMaxLonF = lonend
res@mpCenterLonF = (lonstr+lonend)/2

res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False

res@gsnStringFontHeightF = 0.03
res@tmXBLabelFontHeightF = 0.02
res@tmYLLabelFontHeightF = 0.02
; ;------------------------------------------

;   vcres                         = True             ; plot mods desired
;   vcres@gsnFrame                = False            ; so we can draw time stamp
;   vcres@gsnDraw                 = False
;   vcres@vcRefAnnoOrthogonalPosF =  -0.4 ;-1.0             ; move ref vector up
;   ;vcres@vcGlyphStyle            = "CurlyVector"  ;"LineArrow"    ; turn on curly vectors
;   ;vcres@vcMinDistanceF          = 0.017
;   vcres@vcGlyphStyle            = "LineArrow"
;   vcres@vcMinDistanceF          = 0.03
;   vcres@vcRefMagnitudeF         = 10.0             ; define vector ref mag
  
;   vcres@vcLineArrowThicknessF   = 3.0
;   ;vecres@vcLineArrowHeadMinSizeF = 0.008
;   ;vecres@vcLineArrowHeadMaxSizeF = 0.3

;   vcres@vcRefAnnoString2On = False
;   vcres@vcRefAnnoFontHeightF = 0.015

; ;================================================================================
res@lbLabelBarOn = False
res@cnLevelSelectionMode = "ManualLevels"
; res@cnLevelSpacingF = 10.
; res@cnMaxLevelValF = 150.
; res@cnMinLevelValF = -150.
res@cnLevelSpacingF = 5
res@cnMaxLevelValF = 50
res@cnMinLevelValF = -50
res@cnFillPalette = "BlueWhiteOrangeRed"

; vcres@vcRefMagnitudeF         = 4.0             ; define vector ref mag
; vcres@vcRefLengthF            = 0.035            ; define length of vec ref

leftsti = (/"Day-"+days(2),"Day-"+days(1),"Day-"+days(0),"Day0","Day+"+days(0),"Day+"+days(1),"Day+"+days(2)/)

do i = 0, 6
res@gsnLeftString = leftsti(i)
plot(i) = gsn_csm_contour_map_ce(wks, hx(i,:,:), res)
plot_topo(i) = gsn_csm_contour(wks, tpmask, res_topo)
overlay(plot(i), plot_topo(i))
end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelMainString = "H"+LEV
pres@gsnPanelLabelBar = True
;pres@pmLabelBarOrthogonalPosF = -0.01
pres@lbLabelFontHeightF = 0.015

gsn_panel(wks, plot, (/7,1/), pres)

; plot(0) = gsn_csm_contour_map_ce(wks, hx, res)
; plot_topo(0) = gsn_csm_contour(wks, tpmask, res_topo)
; overlay(plot(0), plot_topo(0))
; draw(plot(0))



end