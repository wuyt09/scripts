
;****************************************************************************
; Aim: Draw all the base/warm.dat for further checking.
;
;            Wu Yuting 2017-8-24
;****************************************************************************

begin

; latstr = 20
; latend = 50
; lonstr = 60
; lonend = 110

path = "../data/"
diri = "/home/yangsong3/data-observation/ERA5-monthly/pressure/"
 
 f = addfile(diri+"t_mon_1x1_1979_2019.nc","r")

 lat  = f->latitude ;(::-1)
 lon  = f->longitude
 lev  = f->level
 print(lev) 
; cc   = short2flt(f->cc)

 delete(f)

 nlat = dimsizes(lat)
 mlon = dimsizes(lon)
 nlev = dimsizes(lev)
 ;totallev = 6+6*nlev

 nvname1 = (/"solarin","ssrd","ssru","t2m","huss","sp","sshf","slhf","o3","t","hus","cc","clwc","ciwc"/)
 ;nvname2 = (/"baseline","warm","albedo","cloud","o3","t","wv","co2","solar"/)

 n1 = dimsizes(nvname1)
 ;n2 = dimsizes(nvname2)

 setfileoption("bin", "ReadByteOrder", "Native")

 ;---------------------------------------
 ; read dataset (*base.dat & *warm.dat)
 ;---------------------------------------
 
 UNDEF = -32767   ;cc@_FillValue 
 data1  = new((/8,nlat,mlon/),float,UNDEF)
 data2  = new((/8,nlat,mlon/),float,UNDEF)
 data_b = new((/6,nlev,nlat,mlon/),float,UNDEF)
 data_w = new((/6,nlev,nlat,mlon/),float,UNDEF)

 do i = 0, 7 ; only the surface level (nvname1, the first six)
   f1 = path+nvname1(i)+"_base.dat"
   data1(i,:,:) = fbindirread(f1,0,(/nlat,mlon/),"float")
   f2 = path+nvname1(i)+"_warm.dat"
   data2(i,:,:) = fbindirread(f2,0,(/nlat,mlon/),"float")
 end do

 data1!1 = "lat"
 data1!2 = "lon"
 data1&lat = lat
 data1&lon = lon
 data1&lat@units = "degrees_north"
 data1&lon@units = "degrees_east"
 copy_VarMeta(data1,data2)

 delete(f1)
 delete(f2)
 
 dd1 = data2 - data1
 copy_VarCoords(data1, dd1)

 do i = 0, 5 
   f1 = path+nvname1(i+8)+"_base.dat"
   f2 = path+nvname1(i+8)+"_warm.dat"
   data_b(i,:,:,:) = fbindirread(f1,0,(/nlev,nlat,mlon/),"float")
   data_w(i,:,:,:) = fbindirread(f2,0,(/nlev,nlat,mlon/),"float")
 end do 

 data_b!1 = "lev"
 data_b!2 = "lat"
 data_b!3 = "lon"
 data_b&lat = lat
 data_b&lon = lon
 data_b&lev = lev
 data_b&lat@units = "degrees_north"
 data_b&lon@units = "degrees_east"
 
 copy_VarMeta(data_b, data_w)
 
 dd2 = data_w - data_b
 copy_VarCoords(data_b, dd2)
 print("end of *_base/warm.data")
 


 ;-------------------------------------------
 ; draw all *_base/warm.dat
 ;-------------------------------------------
  
  figpath = "/home/yangsong3/wuyt/sysu/figures/TP_surrounding/cfram20210212/ncl/"
  figname = "check_diff_before_fortran"
  wks  = gsn_open_wks("pdf",figpath+figname)
  gsn_define_colormap(wks,"BlueWhiteOrangeRed")
  plot_b = new(dimsizes(nvname1),graphic)
  plot_w = new(dimsizes(nvname1),graphic)
  plot_d = plot_b


  res   = True
  res@gsnDraw  = False
  res@gsnFrame = False
  res@gsnSpreadColors = True
  res@gsnAddCyclic = True

  res@cnFillOn      = True
  res@cnLinesOn     = False
  res@cnInfoLabelOn = False

  res@mpOutlineOn   = True
  res@mpFillOn      = True
  ;res@mpCenterLonF  = 0.
  ; res@mpMinLatF = latstr
  ; res@mpMaxLatF = latend
  ; res@mpMinLonF = lonstr
  ; res@mpMaxLonF = lonend

  res@gsnStringFontHeightF = 0.035
  res@tmXBLabelFontHeightF = 0.025
  res@tmYLLabelFontHeightF = 0.025

  do i = 0, 7
    res@gsnLeftString = nvname1(i)
    plot_b(i) = gsn_csm_contour_map(wks,data1(i,:,:),res)
    plot_w(i) = gsn_csm_contour_map(wks,data2(i,:,:),res); write loop content
    plot_d(i) = gsn_csm_contour_map(wks,dd1(i,:,:),res)
  end do
  
  do i = 8, 13 
        res@gsnLeftString  = nvname1(i)
        res@gsnRightString = "500hPa"
        plot_b(i) = gsn_csm_contour_map(wks, data_b(i-8,{500},:,:), res)
        plot_w(i) = gsn_csm_contour_map(wks, data_w(i-8,{500},:,:), res)  
        plot_d(i) = gsn_csm_contour_map(wks, dd2(i-8,{500},:,:), res)      
  end do

  ;----------------------------------

   pres  = True
   pres@gsnMaximize = True
   
   pres@tiMainString = "base state"
   gsn_panel(wks,plot_b,(/3,5/),pres)

   pres@tiMainString = "warm state"
   gsn_panel(wks,plot_w,(/3,5/),pres)

   pres@tiMainString = "Diff"
   gsn_panel(wks,plot_d,(/3,5/),pres)






end 
