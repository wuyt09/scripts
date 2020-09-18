
begin

  lonstr = 60 ;10
  lonend = 110 ;150
  latstr = 20 ;-10
  latend = 50 ;70


nvname = (/"albedo","wv","cloud","cloud_sw","cloud_lw","co2","o3","solar",\
           "dyn","atm_dyn","ocean_dyn","shflx","lhflx","error"/)  ; ,"error"

;================read metadata=============================================
diri  = "/home/yangsong3/data-model/wuyt/TPbias/TPbias_Ctrl_FAMIP/a20191206/" 
file1 = "Ctrl_FAMIP.cam.h0.T.NoInterpolation.1979-2005.nc"
f1 = addfile(diri+file1,"r")
lon = tofloat(f1->lon)
lat = tofloat(f1->lat)
level = tofloat(f1->lev)

nlat = dimsizes(lat)
nlon = dimsizes(lon)
nlev = 27
nvars = dimsizes(nvname)
delete(f1)
;===========================================================================

diri2 = "/home/yangsong3/wuyt/sysu/scripts/TP_coldbias/model/a20200910/cfram20200917/"
; file2 = diri2+"forcing_SR.grd"
file2 = diri2+"partial_T_SR.grd"
file3 = diri2+"data/t2m_warm.dat"
file4 = diri2+"data/t2m_base.dat"


lev = new(nlev,float,"No_FillValue")
lev(0:nlev-2) = level
lev(nlev-1) = 1001

;setfileoption("bin","ReadByteOrder","Native")

UNDEF = -999.
x = new((/nvars,nlev,nlat,nlon/),float,UNDEF)
x!0 = "var"
x!1 = "lev"
x!2 = "lat"
x!3 = "lon"
x&var = nvname
x&lat = lat
x&lon = lon
x&lev = lev
x&lat@units = "degrees_north"
x&lon@units = "degrees_east"
x&lev@units = "hPa"
x@long_name = "forcing"
x@units     = "W/m2"

do nv = 0,nvars-1
   do nl = 0,nlev-1
      x(nv,nl,:,:) = fbindirread(file2,nv*nlev+nl,(/nlat,nlon/),"float")
   end do
end do

if (any(isnan_ieee(x))) then   ; if x got a element is Not A Number, then the Fill_value
  value = -999.
  replace_ieeenan (x,value,0)
  x@_FillValue = value
end if
printVarSummary(x)
printMinMax(x(0,:,:,:),1)


; ;===============================================================================================
obs = new((/2,nlat,nlon/),float,UNDEF)
obs(0,:,:) = fbindirread(file3,0,(/nlat,nlon/),"float")
obs(1,:,:) = fbindirread(file4,0,(/nlat,nlon/),"float")
copy_VarMeta(x(0,0,:,:),obs(0,:,:))
printVarSummary(obs)
print("end of data read")



title = (/"t_wv","t_cloud","t_dyn"/)
tx = new((/3,nlat,nlon/), float)
tx(0,:,:) = x(1,nlev-1,:,:) + obs(1,:,:)
tx(1,:,:) = x(2,nlev-1,:,:) + obs(1,:,:)
tx(2,:,:) = x(9,nlev-1,:,:) + obs(1,:,:)
tx!0 = "var"
tx&var = title
copy_VarCoords(x(0,nlev-1,:,:), tx(0,:,:))
printVarSummary(tx)

; t_wv = x(1,nlev-1,:,:) + obs(1,:,:)

outpath = "/home/yangsong3/wuyt/sysu/scripts/TP_coldbias/model/a20200910/cfram_decomposeT/"
; outname = "t_wv.dat"
; system("rm -f "+outpath+outname)
; fbindirwrite(outpath+outname, t_wv)
; i = 0
; setfileoption ("bin", "WriteByteOrder", "BingEndian")
do i = 0, dimsizes(title)-1
  outname = title(i)+".dat"
  system("rm -f "+outpath+outname)
  fbindirwrite(outpath+outname, tx(i,:,:))
end do
exit
;==========================================================
;****************************************************
; plot
;****************************************************

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/a20200910/" ;diri2
figname = "write_t_cloud-wv-dyn"
wks  = gsn_open_wks("pdf",figpath+figname)
gsn_define_colormap(wks,"temp_diff_18lev")

nplot = dimsizes(title)
plot = new(nplot, graphic)
plot_topo = plot
;=======================================
; this is for outline of TP
;=======================================

ff = "~/wuyt/data/TOPO_1500.nc"
topo_in = addfile(ff,"r")
topo = topo_in ->TP

mlat = topo_in->LATIXY(:,0)
mlon = topo_in->LONGXY(0,:)
topo!0 = "lsmlat"
topo!1 = "lsmlon"
topo&lsmlat = mlat
topo&lsmlon = mlon

;----------------------
res_topo            = True
res_topo@cnFillOn   = False
res_topo@cnLinesOn  = True
res_topo@gsnFrame   = False
res_topo@gsnDraw    = False
res_topo@cnLineLabelsOn         = False
res_topo@cnLineThicknessF       = 2
res_topo@cnLineColor            = "grey32"
;res_topo@cnLineDashPattern      = "1"
res_topo@cnInfoLabelOn          = False
res_topo@cnLevelSelectionMode   ="ExplicitLevels"
res_topo@cnLevels               = (/1500/)
;==========================================================

res = True
res@gsnDraw   = False
res@gsnFrame  = False
;res@gsnSpreadColors = True
res@gsnAddCyclic = False

res@mpFillOn  = True
;res@mpMinLonF = 180.
;res@mpMaxLonF = 360.
res@mpMinLatF = latstr
res@mpMaxLatF = latend
res@mpMinLonF = lonstr
res@mpMaxLonF = lonend
;res@mpCenterLonF = 85.
res@mpGeophysicalLineThicknessF = 2

res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False

;res@lbLabelAutoStride = True
;res@lbLabelFontHeightF = 0.01
;res@pmLabelBarOrthogonalPosF = 0.16  ; 水平位置
;res@pmLabelBarHeightF = 0.02    ;色标高度
;res@pmLabelBarWidthF = 0.7
;res@pmLabelBarHeightF = 0.15
;res@lbLabelFontHeightF = 0.016
res@lbLabelBarOn         = False 
res@lbLabelsOn           = False

res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025
res@tmXTOn = False
res@tmYROn = False

;res@gsnLeftStringFontHeightF  = 0.02
;res@gsnRightStringFontHeightF  = 0.02
res@gsnStringFontHeightF = 0.035
res@gsnRightString = "K"

;-------------------------------------------------------

; res@cnLevelSelectionMode = "ExplicitLevels"
; res@cnLevels = (/-3,-2.5,-2,-1.5,-1,-0.5,-0.1,0.1,0.5,1,1.5,2,2.5,3/)
; res@cnFillColors = (/4,5,6,7,8,9,10,11,12,13,14,15,16,17,18/)

res@cnLevelSelectionMode = "ManualLevels"
res@cnMaxLevelValF = 260
res@cnMinLevelValF = 250
res@cnLevelSpacingF = 1

do i = 0, dimsizes(title)-1
    res@gsnLeftString = title(i)
    plot(i) = gsn_csm_contour_map(wks, tx(i,:,:), res)
    plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(i), plot_topo(i))
end do

; i = dimsizes(title)-1
; res@gsnLeftString = title(i)
; plot(i) = gsn_csm_contour_map(wks, x_diff(i,:,:), res)
; plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
; overlay(plot(i), plot_topo(i))

pres = True
pres@gsnMaximize = True
pres@gsnFrame    = True
pres@gsnDraw     = True
pres@gsnPanelLabelBar = True

pres@lbLabelFontHeightF = 0.015
pres@pmLabelBarWidthF = 0.8    ;plot manager 
pres@pmLabelBarHeightF = 0.06

gsn_panel(wks,plot,(/1,3/),pres)

end 
