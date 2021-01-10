
begin

yrStrt = 1979
yrLast = 2005
mnStrt = 1
mnLast = 12

year = ispan(yrStrt, yrLast, 1)
month = ispan(mnStrt, mnLast, 1)
penta = ispan(1, 73, 1)
nyear = dimsizes(year)
nmonth = dimsizes(month)
npenta = dimsizes(penta)

nlat = 96
nlon = 144
nlev = 17
;===================================================


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

;===================================================

pens = 18
diri = "/home/yangsong3/wuyt/sysu/scripts/TP_FAMIP/"

casename = (/"Ctrl_FAMIP","Hist_FAMIP","Hist_PO"/) ;,"Hist_FAMIP_CMIP6"
ncase = dimsizes(casename)

u_peak = new((/ncase,nlat,nlon/), float)

do i = 0, ncase-1
    fili = "2_MTG_onset_pentad_"+yrStrt+"_"+yrLast+"_"+casename(i)+".nc"
    f = addfile(diri+fili, "r")
    var_med = f->MTG
    printVarSummary(var_med)
    u_peak(i,:,:) = var_med
end do
ISM = wgt_areaave(u_peak(:,{latS1:latN1},{lonL1:lonR1}), 1., 1., 0)
BOB = wgt_areaave(u_peak(:,{latS2:latN2},{lonL2:lonR2}), 1., 1., 0)
SCS = wgt_areaave(u_peak(:,{latS3:latN3},{lonL3:lonR3}), 1., 1., 0)
; ==============================================================
;  plot (a)
; ==============================================================

 figpath = "/home/yangsong3/wuyt/sysu/figures/TP_FAMIP/" 
 figname = "3_onset_pentad_grid_MTG"
 wks  = gsn_open_wks("ps",figpath+figname)       ; specifies a ps plot
; gsn_define_colormap(wks,"sunshine_9lev") ; choose colormap
; gsn_define_colormap(wks,"BlWhRe") ; choose colormap
; gsn_define_colormap(wks,"testcmap") ; choose colormap
 gsn_define_colormap(wks,"so4_23")  ; choose colormap
 plot = new(ncase, graphic )

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
 resR@cnMinLevelValF        =  17 ;-6e-5 ;1.6e7 ;2.5         ; set min contour level
 resR@cnMaxLevelValF        =  40 ;6e-5 ;1.6e7 ;2.5          ; set max contour level
 resR@cnLevelSpacingF       =  1 ;1e-5 ;0.5 ; set contour interval

 resR@lbLabelBarOn         = False
 resR@lbLabelStride        = 1        ; label bar stride
; resR@lbOrientation        = "Vertical"
 resR@lbLabelAlignment       = "BoxCenters"
 resR@lbTopMarginF          = 0.2

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
   txid = new((/ncase,3/),graphic)                  ; Nine text strings
   txres               = True
   txres@txFont        = "helvetica"
;   txres@txFontHeightF = 0.01
;   txres@txFont        = "duplex_roman"
   txres@txFontHeightF = 0.02 ;0.014
   txres@txFontColor    = 1
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

  dum1 = new((/ncase,4/),graphic)
  dum2 = new((/ncase,4/),graphic)
  dum3 = new((/ncase,4/),graphic)

  do k = 0,ncase-1 ;4-6月
      resR@gsnLeftString  = casename(k)
      plot(k) = gsn_csm_contour_map(wks,u_peak(k,:,:),resR) ; create a default plot

       txid(k,0) = gsn_add_text (wks,plot(k),sprintf("%4.2f", BOB(k)), 95.0,  10.0 ,txres)
       txid(k,1) = gsn_add_text (wks,plot(k),sprintf("%4.2f", SCS(k)), 115.0, 15.0,txres)
       txid(k,2) = gsn_add_text (wks,plot(k),sprintf("%4.2f", ISM(k)), 72.5,  15.0,txres)

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
  gsn_panel(wks,plot,(/ncase,1/),resP)               ; now draw as one plot
  frame(wks)   
end