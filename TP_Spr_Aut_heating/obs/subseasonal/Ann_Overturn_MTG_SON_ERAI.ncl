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

  nlon = 240
  nlat = 121
   nyr  = 38
   npt  = 18
   npnt = 30
  x      = ispan(1,12,1)

; ==============================================================
; (I) Basic Field ==> VP200
; ==============================================================
   ; f3  = addfile("./data/sp.interim.1979.2017.nc","r")
   ; PS0 = short2flt(f3->sp(:,:,:))
   f3  = addfile("/home/yangsong3/data-observation/wuyt/ERAI/monthly/1x1/surf_pressure.mon.mean.7916.r1.nc","r")
   PS0 = short2flt(f3->sp(0:455,:,:))

   fils  = systemfunc("ls /home/yangsong3/data-observation/ERA-interim/pressure/air/air.interim.*.nc")
   ;print(fils)
   f  = addfiles(fils,"r")
   t0 = short2flt(f[:]->air)
   printVarSummary(t0)

   Ps0 = linint2(PS0&longitude, PS0&latitude, PS0, True, t0&lon, t0&lat, 0)
   ps0 = clmMonTLL(Ps0)
   printVarSummary(ps0)


   yr   = new(nyr,integer,-999)
   mn   = new((/nyr,12/),integer,-999)
   dy   = new((/nyr,31/),integer,-999)

   yr   = ispan( 1979,2016,1 )
   do n=0,nyr-1
	   mn(n,:)   = ispan(1,12,1)
	   dy(n,:)   = ispan(1,31,1)
   end do
   ddd4   = day_of_year(yr,mn(:,7),dy(:,1))	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;要减一！！！
   ddd6   = day_of_year(yr,mn(:,11),dy(:,28))   ;19-36pentad 18 ;13-42Pentad 30
   ddd12  = day_of_year(yr,mn(:,11),dy(:,30))
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
   t_pnt = dim_avg_n_Wrap(t_Pnt, 0)

; MTG 有限中央差分
   longitude  = f3->longitude
   dlon = (longitude(2)-longitude(1))*0.0174533
   dY   = 6371229.0*dlon                             ; dY
   tY   = t_pnt
   tY   = center_finite_diff_n (t_pnt, dY , True,0,2)
   MTG  = tY
   ;MTG!0 = "comp"
   printVarSummary(MTG)
   ;delete(Cmp_t)


; 垂直积分
  plev = t_TLLL&level
  plev = plev*100
  psfc = max(plev)
  ptop = min(plev)
  dP = dpres_plevel_Wrap(plev, psfc, ptop, 0)
  dp = conform_dims(dimsizes(MTG), dP, 1)
  MTGDP = MTG
  MTGDP = MTG*dp
  mtg = dim_sum_n_Wrap(MTGDP, 1)  ;/sum(dP)
  Cmp_MTG = mtg

; ==============================================================
; (IV) Reversal Pentad
; ==============================================================

  bob = wgt_areaave_Wrap(Cmp_MTG(:,{latS2:latN2},{lonL2:lonR2}),1.0,1.0,1)
	scs = wgt_areaave_Wrap(Cmp_MTG(:,{latS3:latN3},{lonL3:lonR3}),1.0,1.0,1)
	ism = wgt_areaave_Wrap(Cmp_MTG(:,{latS1:latN1},{lonL1:lonR1}),1.0,1.0,1)
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
		diff(n,:,:) = Cmp_MTG(n,:,:)*Cmp_MTG(n-1,:,:)
		diff1(n) = bob(n)*bob(n-1)
		diff2(n) = scs(n)*scs(n-1)
		diff3(n) = ism(n)*ism(n-1)
  end do
	do n=0,npnt-2
		last(n,:,:) = Cmp_MTG(n,:,:)*Cmp_MTG(n+1,:,:)
 		last1(n) = bob(n)*bob(n+1)
		last2(n) = scs(n)*scs(n+1)
		last3(n) = ism(n)*ism(n+1)
  end do
;************************************************
	u_peak = Cmp_MTG(0,:,:)
	BOB = new( 1 , integer )
	SCS = BOB
	ISM = BOB

	;do k=0,1
	do j=0,nlat-1
		do i=0,nlon-1
	;		a= ind(MTG(:,j,i).gt.0)
			a= ind((diff(:,j,i).lt.0).and.(Cmp_MTG(:,j,i).lt.0).and.(last(:,j,i).gt.0))
			u_peak(j,i)=(/a(0)/)	+42
			delete(a)
		end do
	end do
		a= ind((diff1.lt.0).and.(bob.lt.0).and.(last1.gt.0))
		BOB =(/a(0)/)	+42
		delete(a)
		a= ind((diff2.lt.0).and.(scs.lt.0).and.(last2.gt.0))
		SCS =(/a(0)/)	+42
		delete(a)
		a= ind((diff3.lt.0).and.(ism.lt.0).and.(last3.gt.0))
		ISM =(/a(0)/)	+42
		delete(a)

	print(u_peak({10},{50}))


	print(BOB)
	print(SCS)
	print(ISM)

; ==============================================================
;  plot (a)
; ==============================================================
 figpath = "/home/yangsong3/wuyt/sysu/figures/TP_Spr_Aut_heating/obs/subseasonal/"
 figname = "Ann_Overturn_MTG_SON_ERAI"

 wks  = gsn_open_wks("pdf",figpath+figname)       ; specifies a ps plot
; gsn_define_colormap(wks,"sunshine_9lev") ; choose colormap
; gsn_define_colormap(wks,"BlWhRe") ; choose colormap
; gsn_define_colormap(wks,"testcmap") ; choose colormap
 gsn_define_colormap(wks,"so4_23")  ; choose colormap
 plot = new(1, graphic )

 resR                       = True     ; plot mods desired
 resR@gsnDraw               = False
 resR@gsnFrame              = False
 resR@gsnSpreadColors       = True     ; use full range of color map
 resR@gsnSpreadColorStart   = 24 ;40       ; choose starting color
 resR@gsnSpreadColorEnd     = 2 ;-40      ; choose ending color
 resR@gsnLeftStringFontHeightF =0.012
 resR@gsnRightStringFontHeightF =0.012
 resR@gsnRightString = "MTG (P46-P64)" ;"ERA-40: JJA" ;
;  resR@gsnMaximize          = True
;  resR@gsnStringFontHeightF     = 18    ;左右角标字体大小
  resR@tmYLLabelFontHeightF = 0.012     ;坐标字体大小
  resR@tmXBLabelFontHeightF = 0.012

 resR@cnFillOn              = True     ; turn on color
 resR@cnFillMode            = "CellFill"
 resR@cnLinesOn             = False    ; turn off contour lines
 resR@cnLineLabelsOn       = False
 resR@cnInfoLabelOn        = False
 resR@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
 resR@cnMinLevelValF        =  44 ;-6e-5 ;1.6e7 ;2.5         ; set min contour level
 resR@cnMaxLevelValF        =  70 ;6e-5 ;1.6e7 ;2.5          ; set max contour level
 resR@cnLevelSpacingF       =  1 ;1e-5 ;0.5 ; set contour interval

 resR@lbLabelBarOn         = True
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
   txid = new((/3/),graphic)                  ; Nine text strings
   txres               = True
   txres@txFont        = "helvetica"
;   txres@txFontHeightF = 0.01
;   txres@txFont        = "duplex_roman"
   txres@txFontHeightF = 0.02 ;14
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
  resp@gsLineThicknessF = 2.0                       ; thickness of lines
  resp@gsLineDashPattern = 2
;  resp@gsLineLabelString= "test"                    ; adds a line label string

; create array of dummy graphic variables. This is required, b/c each line
; must be associated with a unique dummy variable.

  dum1 = new(4,graphic)
  dum2 = new(4,graphic)
  dum3 = new(4,graphic)
    
    resR@gsnLeftString = "ERAI"
    resR@gsnRightString = "MTG (P46-P64)"
    resR@gsnStringFontHeightF = 0.04
	  plot(0) = gsn_csm_contour_map(wks,u_peak,resR) ; create a default plot

	   txid(0) = gsn_add_text (wks,plot(0),BOB, 95.0,  10.0 ,txres)
	   txid(1) = gsn_add_text (wks,plot(0),SCS, 115.0, 15.0,txres)
	   txid(2) = gsn_add_text (wks,plot(0),ISM, 72.5,  15.0,txres)

	; draw each line separately. Each line must contain two points.
	  do i = 0 , 3
		dum1(i)=gsn_add_polyline(wks,plot(0),xpts1(i:i+1),ypts1(i:i+1),resp)
	  end do
	  do i = 0 , 3
		dum2(i)=gsn_add_polyline(wks,plot(0),xpts2(i:i+1),ypts2(i:i+1),resp)
	  end do
	  do i = 0 , 3
		dum3(i)=gsn_add_polyline(wks,plot(0),xpts3(i:i+1),ypts3(i:i+1),resp)
	  end do
draw(plot)
end
