
begin

lonstr = 50
lonend = 180
latstr = 0
latend = 30

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

casename = (/"CTRL","NOTPQ1_MAY"/)
ncase = dimsizes(casename)
pens = 13 ;19
pene = 42 ;36
npen = pene - pens + 1
nlat = 96
nlon = 144
nlev = 37
clim_pentad_ps = new((/ncase,npen,nlat,nlon/), float)
clim_pentad_t = new((/ncase,npen,nlev,nlat,nlon/), float)

do i = 0, ncase-1
  diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(i)+"_F2000_F19G16_CAM4/pro/daily/"

  pname = "F2000_"+casename(i)+".PS.daily.06-30.nc"
  fp = addfile(diri+pname, "r")

  ps = fp->PS
  printVarSummary(ps)

  pentad_ps = ps(:,0:72,:,:)

  do p = 0,72  ; a year has 73 penta
    pentad_ps(:,p,:,:) = dim_sum_n_Wrap(ps(:,5*p:5*p+4,:,:), 1)/5.0
  end do

  clim_pentad_PS = dim_avg_n_Wrap(pentad_ps, 0)
  clim_pentad_ps(i,:,:,:) = clim_pentad_PS(pens-1:pene-1,:,:)
  delete([/ps,pentad_ps,clim_pentad_PS/])

  ;--------------------------------------------------------------------
  tname = systemfunc("ls "+diri+"F2000_"+casename(i)+".plevel.T.*.daily.nc")
  print(tname)

  ft = addfiles(tname, "r")

  T = ft[:]->T(:,::-1,:,:)
  t = reshape(T, (/25,365,37,nlat,nlon/))
  copy_VarCoords(T(0:364,:,:,:), t(0,:,:,:,:))

  pentad_t = t(:,0:72,:,:,:)

  do p = 0,72  ; a year has 73 penta
    pentad_t(:,p,:,:,:) = dim_sum_n_Wrap(t(:,5*p:5*p+4,:,:,:), 1)/5.0
  end do

  clim_pentad_T = dim_avg_n_Wrap(pentad_t, 0)
  clim_pentad_t(i,:,:,:,:) = clim_pentad_T(pens-1:pene-1,:,:,:)
  delete([/T,t,pentad_t,clim_pentad_T/])
end do
printVarSummary(clim_pentad_ps)
printVarSummary(clim_pentad_t)
;-------------------------------------------------------------------
; partial T / partial y
rearth = 6371229.0 ; meters

lon = fp->lon
dlon = (lon(2) - lon(1))*0.0174533
dY = tofloat(rearth*dlon) 
;dX = cos(lat*0.0174533)*dlon
printVarSummary(dY)

MTG = clim_pentad_t
MTG = center_finite_diff_n(clim_pentad_t, dY, True, 0, 3)
MTG!0 = "case"
MTG&case = casename
printVarSummary(MTG)

; vertical integration

lev = MTG&lev
plev = lev*100.
plev@units = "Pa"
plev!0 = "lev"
plev&lev = lev
;psfc = max(plev)
;ptop = min(plev)
;dP = dpres_plevel_Wrap(plev, psfc, ptop, 0)
;dp = conform_dims(dimsizes(MTG), dP, 1)

;tdp = MTG
;tdp = MTG*dp
;printVarSummary(tdp)

;mtg = dim_sum_n_Wrap(tdp(:,{500:200},:,:), 1)/sum(dP({500:200}))
;copy_VarCoords(MTG(:,0,:,:), mtg)

Cmp_MTG = MTG(:,:,0,:,:)
pbot = 50000
ptop = 20000
do k = 0,1
  do p = 0, npen-1
    xxx  = MTG(case|k,time|p,lat|:,lon|:,lev|:)
    Cmp_MTG(k,p,:,:) = vibeta(plev, xxx, 1, clim_pentad_ps(k,p,:,:), pbot, ptop)
  end do
end do

printVarSummary(Cmp_MTG)


;-------------------------------------------------------------------
; Reversal Pentad

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
  do n=1,npen-1
    diff(:,n,:,:) = Cmp_MTG(:,n,:,:)*Cmp_MTG(:,n-1,:,:)
    diff1(:,n) = bob(:,n)*bob(:,n-1)
    diff2(:,n) = scs(:,n)*scs(:,n-1)
    diff3(:,n) = ism(:,n)*ism(:,n-1)
    end do
  do n=0,npen-3
    last(:,n,:,:) = Cmp_MTG(:,n,:,:)*Cmp_MTG(:,n+1,:,:)
    last1(:,n) = bob(:,n)*bob(:,n+1)
    last2(:,n) = scs(:,n)*scs(:,n+1)
    last3(:,n) = ism(:,n)*ism(:,n+1)
   end do
;************************************************
  u_peak = Cmp_MTG(:,0,:,:)
  BOB = new(2,integer)
  SCS = BOB
  ISM = BOB
  
  do k = 0, 1
    do j=0,nlat-1
      do i=0,nlon-1
        a= ind((diff(k,:,j,i).lt.0).and.(Cmp_MTG(k,:,j,i).gt.0).and.(last(k,:,j,i).gt.0))
        u_peak(k,j,i)=(/a(0)/)  +pens-1
        delete(a)
      end do
    end do

    a= ind((diff1(k,:).lt.0).and.(bob(k,:).gt.0).and.(last1(k,:).gt.0))
    BOB(k) =(/a(0)/)  +pens
    delete(a)
    a= ind((diff2(k,:).lt.0).and.(scs(k,:).gt.0).and.(last2(k,:).gt.0))
    SCS(k) =(/a(0)/)  +pens
    delete(a)
    a= ind((diff3(k,:).lt.0).and.(ism(k,:).gt.0).and.(last3(k,:).gt.0))
    ISM(k) =(/a(0)/)  +pens
    delete(a)
end do

  print(u_peak(:,{10},{100}))

  print(BOB)
  print(SCS)
  print(ISM)

;-------------------------------------------------------------------

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20200704/"
figname = "onset_pentad_model_"+casename(0)+"_"+casename(1)

wks = gsn_open_wks("ps", figpath+figname)
gsn_define_colormap(wks,"so4_23")  ; choose colormap

nplot = 2
plot = new(nplot, graphic)
Left = (/"Ctrl","NoHeat"/)

 resR                       = True     ; plot mods desired
 resR@gsnDraw               = False
 resR@gsnFrame              = False
 resR@gsnSpreadColors       = True     ; use full range of color map
 resR@gsnSpreadColorStart   = 24 ;40       ; choose starting color
 resR@gsnSpreadColorEnd     = 2 ;-40      ; choose ending color
 resR@gsnLeftStringFontHeightF =0.025
 resR@gsnRightStringFontHeightF =0.025
 resR@gsnRightString = " " ;"ERA-40: JJA" ;
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
 resR@cnMinLevelValF        = pens ;-2 ; 17 ;-6e-5 ;1.6e7 ;2.5         ; set min contour level
 resR@cnMaxLevelValF        = pene ;+2 ;40 ;6e-5 ;1.6e7 ;2.5          ; set max contour level
 resR@cnLevelSpacingF       =  1 ;1e-5 ;0.5 ; set contour interval

 resR@lbLabelBarOn         = False
 resR@lbLabelStride        = 1        ; label bar stride
; resR@lbOrientation        = "Vertical"
 resR@lbLabelAlignment       = "BoxCenters"
 resR@lbTopMarginF      = 0.2

 resR@tmYLTickSpacingF  = 10
 resR@tmXBTickSpacingF  = 20

; resR@vpXF      = 0.2        ; In page coordinates, where to start
; resR@vpYF      = 0.90        ; the plot
; resR@vpWidthF         = 0.6

  resR@mpFillOn             = False        ; turn off map fill
  resR@mpMinLatF            = latstr         ; zoom in on map
  resR@mpMaxLatF            = latend
  resR@mpMinLonF            = lonstr
  resR@mpMaxLonF            = lonend
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
   txres@txFontHeightF = 0.02
   txres@txFontColor  = 1
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
  resP@gsnMaximize      = True                  ; don't advance panel plot
  resP@gsnPanelLabelBar = True                   ; add common colorbar
  ; resP@gsnPanelBottom   = 0.02
  resP@txString   = "MTG Reversal Pentad"
  resP@gsnPanelYWhiteSpacePercent = 3   ;图间隔大小
  resP@lbLabelAlignment       = "BoxCenters"
  resP@txFontHeightF = 0.03   ;标题字体大小
  gsn_panel(wks,plot,(/2,1/),resP)               ; now draw as one plot
  frame(wks)   

end