load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin

  yrStrt = 1981
  yrLast = 2014

;************************************************
; EP & CP .dat 
;************************************************
  f1     = addfile("./ssrd_and_ssr_era_mn_197901_201507_2time_4step.nc","r")  	
	
  TIME   = f1->time
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
	
  x      = short2flt(f1->ssrd(iYYYY,:,:))  ;Pa
  y      = short2flt(f1->ssr(iYYYY,:,:))
  print(x(0,:,0)+" "+x(1,:,0)+" "+x(2,:,0)+" "+x(3,:,0))

  ;for ERA-Interim heat flux data special
     x11 = x(0::8,:,:)
     x12 = x(1::8,:,:)
     x13 = x(2::8,:,:)
     x14 = x(3::8,:,:)
     x21 = x(4::8,:,:)
     x22 = x(5::8,:,:)
     x23 = x(6::8,:,:)
     x24 = x(7::8,:,:)
     x1  = x14/(12*3600)+x24/(12*3600)
;     x1  = (x12-x11)/(3*3600)+(x22-x21)/(3*3600)+(x14-x13)/(3*3600)+(x24-x23)/(3*3600)
     x1  = x1/2

     y11 = y(0::8,:,:)
     y12 = y(1::8,:,:)
     y13 = y(2::8,:,:)
     y14 = y(3::8,:,:)
     y21 = y(4::8,:,:)
     y22 = y(5::8,:,:)
     y23 = y(6::8,:,:)
     y24 = y(7::8,:,:)
     y1  = y14/(12*3600)+y24/(12*3600)
;     y1  = (y12-y11)/(3*3600)+(y22-y21)/(3*3600)+(y14-y13)/(3*3600)+(y24-y23)/(3*3600)
     y1  = y1/2 
     y1  = x1-y1  
 
     lat = f1->latitude
     lon = f1->longitude
     Time= f1->time(iYYYY)
     time= ut_calendar(Time(0::8),-1)
     print(lat)

     x1!0 = "time"
     x1!1 = "lat"
     x1!2 = "lon"
     x1&time = time                
     x1&lat  = lat
     x1&lon  = lon
     x1&lat@units = "degrees_north"
     x1&lon@units = "degrees_east"
     copy_VarMeta(x1,y1)
     delete(f1)

  X1      = month_to_annual (x1, 1)     ;for ssrd
  Y1      = month_to_annual (y1, 1)     ;for ssru
  print(X1(:,0,0))

;  X1      = Y1/X1
  year    = ispan(yrStrt,yrLast,1)
  X1&year = year
  Y1&year = year

  base_e = (/1985,1988,1991,1994,1997,2000,2003,2006,2009/)
  base_s = (/1981,1984,1987,1990,1993,1996,1999,2002,2005/)
  warm_e = (/1990,1993,1996,1999,2002,2005,2008,2011,2014/)
  warm_s = (/1986,1989,1992,1995,1998,2001,2004,2007,2010/)

  n      = dimsizes(base_e)                        ;Returns the dimension sizes of the input variable
  nx     = dimsizes(lat)
  ny     = dimsizes(lon)

  X1_b   = new((/n,5,nx,ny/),typeof(X1))
  X1_w   = new((/n,5,nx,ny/),typeof(X1))  
  Y1_b   = new((/n,5,nx,ny/),typeof(X1))
  Y1_w   = new((/n,5,nx,ny/),typeof(X1))

  do   i = 0,n-1
     X1_b(i,:,:,:) = X1({year|base_s(i):base_e(i)},lat|:,lon|:)
     X1_w(i,:,:,:) = X1({year|warm_s(i):warm_e(i)},lat|:,lon|:)
     Y1_b(i,:,:,:) = Y1({year|base_s(i):base_e(i)},lat|:,lon|:)
     Y1_w(i,:,:,:) = Y1({year|warm_s(i):warm_e(i)},lat|:,lon|:)
  end do
  
  x1_b   = dim_avg_n_Wrap(X1_b,1)
  x1_w   = dim_avg_n_Wrap(X1_w,1)
  d_x1   = x1_w-x1_b
  y1_b   = dim_avg_n_Wrap(Y1_b,1)
  y1_w   = dim_avg_n_Wrap(Y1_w,1)
  d_y1   = y1_w-y1_b

  copy_VarCoords(X1(0:n-1,:,:),x1_b)
  copy_VarCoords(X1(0:n-1,:,:),x1_w)
  copy_VarCoords(X1(0:n-1,:,:),d_x1)
  copy_VarCoords(X1(0:n-1,:,:),y1_b)
  copy_VarCoords(X1(0:n-1,:,:),y1_w)
  copy_VarCoords(X1(0:n-1,:,:),d_y1)


  print(x1_b(0,:,11))

;---Remove file just in case
  file1 = "ssrd_base.dat"
  system("rm -f " + file1)
;---Write first variable to file
  fbindirwrite(file1,x1_b)
  system("ls -l " + file1) 
  
;---Remove file just in case
  file2 = "ssrd_warm.dat"
  system("rm -f " + file2)
;---Write first variable to file
  fbindirwrite(file2,x1_w)
  system("ls -l " + file2)    
  
;---Remove file just in case
  file1 = "ssru_base.dat"
  system("rm -f " + file1)
;---Write first variable to file
  fbindirwrite(file1,y1_b)
  system("ls -l " + file1)

;---Remove file just in case
  file2 = "ssru_warm.dat"
  system("rm -f " + file2)
;---Write first variable to file
  fbindirwrite(file2,y1_w)
  system("ls -l " + file2)

;************************************************
; create plots
;************************************************
  wks = gsn_open_wks("pdf","ssrd_and_ssru")               ; open a ps file
  gsn_define_colormap(wks,"MPL_YlGnBu")
  plot = new(n,graphic)                          ; create a plot array

  res                  = True                      ; plot mods desired
  res@gsnFrame         = False
  res@gsnDraw          = False
  res@gsnSpreadColors  = True                      ; use full colormap

  res@mpOutlineOn      = True
  res@mpFillOn         = True
  res@mpMinLatF        = -90.	; zoom in on a subregion
  res@mpMaxLatF        = 90.
  res@mpMinLonF        = 0.
  res@mpCenterLonF     = 180.
  res@mpMaxLonF        = 360.
  res@mpGeophysicalLineThicknessF=1.3

  res@cnFillOn         = True                      ; turn on color
  res@cnLinesOn        = False 
  res@cnLineLabelsOn   = False                     ; no contour line labels
  res@cnInfoLabelOn    = False                     ; no info label
  res@lbLabelBarOn     = False

;  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
;  res@cnMinLevelValF       = 0.
;  res@cnMaxLevelValF       = 1.
;  res@cnLevelSpacingF      = .1

; now change the size of the tickmark labels
;  res@tmXBLabelFontHeightF = 0.025                 ; resize tick labels
;  res@tmYLLabelFontHeightF = 0.025
  res@gsnRightString       =  "1" ; "W/m~S~2~N~" 

  do i = 0,n-1
     res@gsnLeftString     = "SSRD("+base_s(i)+"-"+base_e(i)+")" 
     plot(i) = gsn_csm_contour_map_ce(wks,x1_b(i,:,:),res)
  end do
  
;************************************************
; create panel
;************************************************
  resP                     = True
  resP@gsnFrame            = False
  resP@gsnMaximize         = True
  resP@gsnPanelLabelBar    = True                ; add common colorbar
  gsn_panel(wks,plot,(/3,3/),resP)  
  frame(wks)       

;************************************************
; second page
;************************************************
  do i = 0,n-1
     res@gsnLeftString     = "SSRD("+warm_s(i)+"-"+warm_e(i)+")"
     plot(i) = gsn_csm_contour_map_ce(wks,x1_w(i,:,:),res)
  end do
  gsn_panel(wks,plot,(/3,3/),resP)
  frame(wks)

;************************************************
; third page
;************************************************
   do i = 0,n-1
      res@gsnLeftString     = "SSRU("+base_s(i)+"-"+base_e(i)+")"
      plot(i) = gsn_csm_contour_map_ce(wks,y1_b(i,:,:),res)
   end do
;   gsn_panel(wks,plot,(/3,3/),resP)
   frame(wks)

;************************************************
; forth page
;************************************************
  do i = 0,n-1
     res@gsnLeftString     = "SSRU("+warm_s(i)+"-"+warm_e(i)+")"
     plot(i) = gsn_csm_contour_map_ce(wks,y1_w(i,:,:),res)
  end do
;  gsn_panel(wks,plot,(/3,3/),resP)
  frame(wks)

;************************************************
; third page
;************************************************
  gsn_define_colormap(wks,"BlueWhiteOrangeRed")  

;  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
;  res@cnMinLevelValF       = -.1              
;  res@cnMaxLevelValF       = .1                
;  res@cnLevelSpacingF      =   0.01               

  do i = 0,n-1
     res@gsnLeftString     = "Delta SSRD"+warm_s(i)+"/"+base_s(i)+")"
     plot(i) = gsn_csm_contour_map_ce(wks,d_x1(i,:,:),res)
  end do
;  gsn_panel(wks,plot,(/3,3/),resP)
  frame(wks)

  do i = 0,n-1
     res@gsnLeftString     = "Delta SSRU"+warm_s(i)+"/"+base_s(i)+")"
     plot(i) = gsn_csm_contour_map_ce(wks,d_y1(i,:,:),res)
  end do
;  gsn_panel(wks,plot,(/3,3/),resP)
  frame(wks)


delete(wks)
end 
