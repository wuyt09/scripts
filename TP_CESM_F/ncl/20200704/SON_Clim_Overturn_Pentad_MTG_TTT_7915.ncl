;*************************************************
; NCL tutorial script: vert_1.ncl
;************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;************************************************
begin

   latS = 0.
   latN = 30.
   lonL = 50.
   lonR = 150.

  latS1  = 10
  latN1  = 20
  lonL1  = 60
  lonR1  = 85.0

  latS2  = 5
  latN2  = 15
  lonL2  = 90
  lonR2  = 100.0

  latS3  = 10
  latN3  = 20
  lonL3  = 110
  lonR3  = 120.0

  nlon = 73
  nlat = 25

  yrStrt = 1979
  yrLast = 2015
  year = ispan(yrStrt, yrLast, 1)
  nyr  = dimsizes(year)

   npt  = 18
   npnt = 30
  x      = ispan(1,12,1)
  ; Left = (/"(a) High_TTT","(b) Low_TTT"/)
  Left = (/"(a) High_Q1","(b) Low_Q1"/)

   TITLE = "SON_Clim_Overturn_Pentad_MTG_Q1_7915"
; ==============================================================
; (I) Basic Field ==> VP200
; ==============================================================
   ; f3  = addfile("./data/sp.interim.1979.2017.nc","r")
   f3  = addfile("/home/yangsong3/data-observation/wuyt/ERAI/monthly/1.5x1.5/sp_era_mn_197901_201507.nc","r")
   PS0 = short2flt(f3->sp(0:431,::-1,:))

   fils  = systemfunc("ls /home/yangsong3/data-observation/ERA-interim/pressure/air/air.interim.19*.nc /home/yangsong3/data-observation/ERA-interim/pressure/air/air.interim.200*.nc /home/yangsong3/data-observation/ERA-interim/pressure/air/air.interim.201[0-5].nc")
   print(fils)

   f  = addfiles(fils,"r")
   t0 = short2flt(f[:]->air)
   printVarSummary(t0)

   Ps0 = linint2_Wrap(PS0&longitude, PS0&latitude, PS0, True, t0&lon, t0&lat, 0)
   ps0 = clmMonTLL(Ps0)
   ps0 = clmMonTLL(PS0)
   printVarSummary(ps0)


   yr   = new(nyr,integer,-999)
   mn   = new((/nyr,12/),integer,-999)
   dy   = new((/nyr,31/),integer,-999)

   ; yr   = ispan( 1979,2016,1 )
   yr = year
   do n=0,nyr-1
	   mn(n,:)   = ispan(1,12,1)
	   dy(n,:)   = ispan(1,31,1)
   end do
   ddd4   = day_of_year(yr,mn(:,6),dy(:,29))	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;要减一！！！
   ddd6   = day_of_year(yr,mn(:,11),dy(:,25))   ;19-36pentad 18 ;13-42Pentad 30
   ddd12  = day_of_year(yr,mn(:,11),dy(:,30))
   print(ddd4)
   print(ddd6)
   print(ddd12)
;************************************************

   t_TLLL  = t0(:,{500:200},::-1,:);13149
   printVarSummary(t_TLLL)
   nlev = dimsizes(t_TLLL(1,:,1,1))
   nlat = dimsizes(t_TLLL(1,1,:,1))
   nlon = dimsizes(t_TLLL(1,1,1,:))

   tYPLLL_AMJ = new((/nyr,npnt*5,nlev,nlat,nlon/),float,-9999)
;   t_Pnt     = t_TLLL(0:1079,:,:,:)

   days = 0
   do n=0,nyr-1
		ddd6(n) = ddd6(n)-1+days
		ddd4(n) = ddd4(n)-1+days
		tYPLLL_AMJ(n,:,:,:,:) = t_TLLL(ddd4(n):ddd6(n),:,:,:)
		days    = days+ddd12(n)
   end do
   printVarSummary(tYPLLL_AMJ)
   
   delete(t0)

   t_Pnt     = tYPLLL_AMJ(:,0:29,:,:,:)
   t_Pnt    = 0
   do n=0,nyr-1
		do k=0,npnt-1
			do i=0,4
;				t_Pnt(n*npnt+k,:,:,:)   = t_Pnt(n*npnt+k,:,:,:) + tYPLLL_AMJ(n,k*5+i,:,:,:)/5.0
				t_Pnt(n,k,:,:,:)   = t_Pnt(n,k,:,:,:) + tYPLLL_AMJ(n,k*5+i,:,:,:)/5.0
			end do
		end do
   end do
   printVarSummary(t_Pnt)

	;p_ClmPnt    = ps0(0:72,:,:)
	;p_ClmPnt    = 0
 	;do n=0,72
 	;	do i=0,4
 	;		p_ClmPnt(n,:,:)=p_ClmPnt(n,:,:)+ps0(n*5+i,:,:)/5.0
	;	end do
	;end do
	;printVarSummary(p_ClmPnt)
  ;print(p_ClmPnt(26,:,:))
; ==============================================================
; (II) PC1 ==> index
; ==============================================================

;  indx = asciiread("SAH_MaxLon_P19P30_7914_Feq90E120E.txt",(/36/),"float")   ;
;  indx1 = asciiread("indx_SST_SCS_110120E_0515N_7914_Mar.txt",(/36/),"float")   ;
;  indx2 = asciiread("indx_SST_SCS_110120E_0515N_7914_Apr.txt",(/36/),"float")   ;
;  indx3 = asciiread("indx_SST_SCS_110120E_0515N_7914_May.txt",(/36/),"float")   ;
  indxx = asciiread("Q1_TP_dtrend_standardize_MAM_SON_ERAI.txt",(/nyr,2/),"float")   ;
  indx = indxx(:,1)
;  print(indx)

; ==============================================================
; (III) Composite
; ==============================================================
  ind1 = ind(indx.ge.1)
  ind2 = ind(indx.le.-1)
 ; ind1 = ind(indx.ge.0.8)
 ; ind2 = ind(indx.le.-0.8)
  print(year(ind1))
  print(year(ind2))
  exit

  num1 = dimsizes(ind1)
  num2 = dimsizes(ind2)
  sp  = t_Pnt ;t_Pnt(ind1,:,:,:,:)
  sn  = t_Pnt(ind2,:,:,:,:)

  Cmp_t = t_Pnt(0:1,:,:,:,:)
  Cmp_t(0,:,:,:,:)= dim_avg_n_Wrap(sp,0)
  Cmp_t(1,:,:,:,:)= dim_avg_n_Wrap(sn,0)
  printVarSummary(Cmp_t)
  delete(t_Pnt)

; MTG 有限中央差分
   longitude  = f3->longitude
   dlon = (longitude(2)-longitude(1))*0.0174533
   dY   = 6371229.0*dlon                             ; dY
   tY   = Cmp_t
   tY   = center_finite_diff_n (Cmp_t, dY , True,0,3)
   MTG  = tY
   MTG!0 = "comp"
   printVarSummary(MTG)
   delete(Cmp_t)
; 垂直积分
	linlog =1
	pbot   =50000
	ptop   =20000
	p = t_TLLL&level
    p = p*100
    print(p)

	Cmp_MTG = MTG(:,:,0,:,:)
	Cmp_MTG = 0
  ppppp = ps0(9,:,:)
	do k=0,1
	do n=0,29
		kkkkk=MTG(comp|k,time|n,lat|:,lon|:,level|:)
		;ppppp=p_ClmPnt(26,:,:)
;		printVarSummary(kkkkk)
;		printVarSummary(ppppp)
		Cmp_MTG(k,n,:,:) = vibeta(p,kkkkk,linlog,ppppp,pbot,ptop)
    end do
	end do
	printVarSummary(Cmp_MTG)



; ==============================================================
; (IV) Reversal Pentad
; ==============================================================

  bob = wgt_areaave_Wrap(Cmp_MTG(:,:,{latS2:latN2},{lonL2:lonR2}),1.0,1.0,1)
	scs = wgt_areaave_Wrap(Cmp_MTG(:,:,{latS3:latN3},{lonL3:lonR3}),1.0,1.0,1)
	ism = wgt_areaave_Wrap(Cmp_MTG(:,:,{latS1:latN1},{lonL1:lonR1}),1.0,1.0,1)
	printVarSummary(bob)

	diff  = Cmp_MTG
	last  = Cmp_MTG
	diff1 = bob
	last1 = bob
	diff2 = bob
	last2 = bob
	diff3 = bob
	last3 = bob
	do n=1,npnt-1
		diff(:,n,:,:) = Cmp_MTG(:,n,:,:)*Cmp_MTG(:,n-1,:,:)
		diff1(:,n) = bob(:,n)*bob(:,n-1)
		diff2(:,n) = scs(:,n)*scs(:,n-1)
		diff3(:,n) = ism(:,n)*ism(:,n-1)
    end do
	do n=0,npnt-2 ;3
		last(:,n,:,:) = Cmp_MTG(:,n,:,:)*Cmp_MTG(:,n+1,:,:)
 		last1(:,n) = bob(:,n)*bob(:,n+1)
		last2(:,n) = scs(:,n)*scs(:,n+1)
		last3(:,n) = ism(:,n)*ism(:,n+1)
   end do
;************************************************
	u_peak = Cmp_MTG(:,0,:,:)
	BOB = new( 2 , integer )
	SCS = BOB
	ISM = BOB

	do k=0,1
		do j=0,nlat-1
		do i=0,nlon-1
	;		a= ind(MTG(:,j,i).gt.0)
			a= ind((diff(k,:,j,i).lt.0).and.(Cmp_MTG(k,:,j,i).lt.0).and.(last(k,:,j,i).gt.0))
			u_peak(k,j,i)=(/a(0)/)	+43
			delete(a)
		end do
		end do
		a= ind((diff1(k,:).lt.0).and.(bob(k,:).lt.0).and.(last1(k,:).gt.0))
		BOB(k) =(/a(0)/)	+43
		delete(a)
		a= ind((diff2(k,:).lt.0).and.(scs(k,:).lt.0).and.(last2(k,:).gt.0))
		SCS(k) =(/a(0)/)	+43
		delete(a)
		a= ind((diff3(k,:).lt.0).and.(ism(k,:).lt.0).and.(last3(k,:).gt.0))
		ISM(k) =(/a(0)/)	+43
		delete(a)
	end do
	print(u_peak(:,{10},{50}))


	print(BOB)
	print(SCS)
	print(ISM)

; ==============================================================
;  plot (a)
; ==============================================================
 figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20200704/" 
 wks  = gsn_open_wks("ps",figpath+TITLE)       ; specifies a ps plot
; gsn_define_colormap(wks,"sunshine_9lev") ; choose colormap
; gsn_define_colormap(wks,"BlWhRe") ; choose colormap
; gsn_define_colormap(wks,"testcmap") ; choose colormap
 gsn_define_colormap(wks,"so4_23")  ; choose colormap
 plot = new(2, graphic )

 resR                       = True     ; plot mods desired
 resR@gsnDraw               = False
 resR@gsnFrame              = False
 resR@gsnSpreadColors       = True     ; use full range of color map
 resR@gsnSpreadColorStart   = 24 ;40       ; choose starting color
 resR@gsnSpreadColorEnd     = 2 ;-40      ; choose ending color
 ; resR@gsnLeftStringFontHeightF =0.012
 ; resR@gsnRightStringFontHeightF =0.012
 resR@gsnStringFontHeightF = 0.02
 resR@gsnRightString = "" ;"ERA-40: JJA" ;MTG (P19-P40)
;  resR@gsnMaximize          = True
;  resR@gsnStringFontHeightF     = 18    ;左右角标字体大小
  resR@tmYLLabelFontHeightF = 0.015     ;坐标字体大小
  resR@tmXBLabelFontHeightF = 0.015

 resR@cnFillOn              = True     ; turn on color
 resR@cnFillMode            = "CellFill"
 resR@cnLinesOn             = False    ; turn off contour lines
 resR@cnLineLabelsOn       = False
 resR@cnInfoLabelOn        = False
 resR@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
 resR@cnMinLevelValF        =  47 ;-6e-5 ;1.6e7 ;2.5         ; set min contour level
 resR@cnMaxLevelValF        =  70 ;6e-5 ;1.6e7 ;2.5          ; set max contour level
 resR@cnLevelSpacingF       =  1 ;1e-5 ;0.5 ; set contour interval

 resR@lbLabelBarOn         = False
 resR@lbLabelStride        = 1        ; label bar stride
; resR@lbOrientation        = "Vertical"
 resR@lbLabelAlignment       = "BoxCenters"
 resR@lbTopMarginF 			= 0.2

 resR@tmYLTickSpacingF  = 10
 resR@tmXBTickSpacingF  = 20

; resR@vpXF      = 0.2        ; In page coordinates, where to start
; resR@vpYF      = 0.90        ; the plot
; resR@vpWidthF         = 0.6

  resR@mpFillOn             = False        ; turn off map fill
  resR@mpMinLatF            = latS         ; zoom in on map
  resR@mpMaxLatF            = latN
  resR@mpMinLonF            = lonL
  resR@mpMaxLonF            = lonR
  resR@mpCenterLonF         = 115. ;180 ;

 resR@tiMainString          = " "
 resR@tiMainFontHeightF  = 0.018


; --------------------------- SAH Track
  resDot = True
  resLine = True

; 绘制00时的点
  resDot@gsMarkerColor = "black"
  resDot@gsMarkerIndex = 1
  resDot@gsMarkerSizeF = 0.02

; 绘制线
;  resLine@gsLineColor = "black"
  resLine@gsLineThicknessF = 5
  resLine@gsLineColor = "black" ;"red"

; add the text
   txid = new((/2,3/),graphic)                  ; Nine text strings
   txres               = True
   txres@txFont        = "helvetica"
;   txres@txFontHeightF = 0.01
;   txres@txFont        = "duplex_roman"
   txres@txFontHeightF = 0.02 ;0.014
   txres@txFontColor 	= 1
;   txres@txFuncCode    = "~"

;************************************************
; add the box
;************************************************
 ypts1 = (/latN1, latN1, latS1, latS1, latN1/)
 xpts1 = (/lonL1, lonR1, lonR1, lonL1, lonL1/)
 ypts2 = (/latN2, latN2, latS2, latS2, latN2/)
 xpts2 = (/lonL2, lonR2, lonR2, lonL2, lonL2/)
 ypts3 = (/latN3, latN3, latS3, latS3, latN3/)
 xpts3 = (/lonL3, lonR3, lonR3, lonL3, lonL3/)

  resp                  = True                      ; polyline mods desired
;  resp@gsLineColor      = "darkpink";"white"                     ; color of lines
  resp@gsLineThicknessF = 4.0                       ; thickness of lines
  resp@gsLineDashPattern = 2
;  resp@gsLineLabelString= "test"                    ; adds a line label string

; create array of dummy graphic variables. This is required, b/c each line
; must be associated with a unique dummy variable.

  dum1 = new((/2,4/),graphic)
  dum2 = new((/2,4/),graphic)
  dum3 = new((/2,4/),graphic)

  do k = 0,1 ;4-6月
	  resR@gsnLeftString  = Left(k)
	  plot(k) = gsn_csm_contour_map(wks,u_peak(k,:,:),resR) ; create a default plot

	   txid(k,0) = gsn_add_text (wks,plot(k),tostring(BOB(k)), 95.0,  10.0 ,txres)
	   txid(k,1) = gsn_add_text (wks,plot(k),tostring(SCS(k)), 115.0, 15.0,txres)
	   txid(k,2) = gsn_add_text (wks,plot(k),tostring(ISM(k)), 72.5,  15.0,txres)

	; draw each line separately. Each line must contain two points.
	  do i = 0 , 3
		dum1(k,i)=gsn_add_polyline(wks,plot(k),xpts1(i:i+1),ypts1(i:i+1),resp)
	  end do
	  do i = 0 , 3
		dum2(k,i)=gsn_add_polyline(wks,plot(k),xpts2(i:i+1),ypts2(i:i+1),resp)
	  end do
	  do i = 0 , 3
		dum3(k,i)=gsn_add_polyline(wks,plot(k),xpts3(i:i+1),ypts3(i:i+1),resp)
	  end do
  end do

;************************************************
; create panel
;************************************************
  resP                  = True                   ; modify the panel plot
  resP@gsnFrame         = False                  ; don't advance panel plot
  resP@gsnPanelLabelBar = True                   ; add common colorbar
  resP@gsnPanelBottom   = 0.02
  resP@txString   = "MTG Reversal Pentad"
  resP@gsnPanelYWhiteSpacePercent = 3   ;图间隔大小
  resP@lbLabelAlignment       = "BoxCenters"
  resP@txFontHeightF = 0.03   ;标题字体大小
  gsn_panel(wks,plot,(/2,1/),resP)               ; now draw as one plot
  frame(wks)                                  ; advance the frame

end
