begin

indexname = (/"Chen","Yang","LYang","Jhun","Sun"/)
casename = (/"CTRL","NOTP","NOQ1","NOTP_NOQ1"/)

diri = "./"
fili1 = "calc_monsoon_index_DJF_CTRL.txt"
fili2 = "calc_monsoon_index_DJF_NOTP.txt"
fili3 = "calc_monsoon_index_DJF_NOQ1.txt"
fili4 = "calc_monsoon_index_DJF_NOTP_NOQ1.txt"

ctrl = asciiread(diri+fili1, (/25,5/), "float")
notp = asciiread(diri+fili2, (/25,5/), "float")
noq1 = asciiread(diri+fili3, (/25,5/), "float")
noei = asciiread(diri+fili4, (/25,5/), "float")
printVarSummary(ctrl)
print(ctrl)

index = new((/4,25,5/), float)
index(0,:,:) = ctrl
index(1,:,:) = notp
index(2,:,:) = noq1
index(3,:,:) = noei

;==================================================================
;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/"
figname = "draw_monsoon_index_DJF"
xaxis = ispan(6, 30, 1)

wks = gsn_open_wks("pdf", figpath+figname)
nplot = dimsizes(indexname)

plot = new(nplot+1, graphic)

res = True
res@gsnDraw  = False
res@gsnFrame = False

res@vpHeightF= 0.4                    ; change aspect ratio of plot
res@vpWidthF = 0.8     

res@gsnLeftStringFontHeightF = 0.03

res@trXMinF = 6
res@trXMaxF = 30
;res@tmYLAutoPrecision = False
;res@tmYLPrecision = 3 ; the number of significant digits is 0
res@xyDashPattern = 0
res@xyLineColors = (/"black","blue","red","green4"/)

do i = 0, dimsizes(indexname)-1
    ;res@trYMinF = -5
    ;res@trYMaxF = 20
    res@gsnLeftString = indexname(i)
    plot(i) = gsn_csm_xy(wks, xaxis, index(:,:,i), res)
end do

;   res_text                    = True                  ; text mods desired
;   res_text@txFontHeightF      = 0.015                 ; change text size
;   res_text@txJust             = "CenterLeft"          ; text justification

;   res_lines                   = True                  ; polyline mods desired
;   res_lines@gsLineDashPattern = 0.                    ; solid line
;   res_lines@gsLineThicknessF  = 5.                    ; line thicker
;   res_lines@gsLineColor       = "red"                 ; line color
;   xx = (/1893,1907/)
;   yy = (/0.705,0.705/)
;   dum1 = gsn_add_polyline(wks,bot_plot,xx,yy,res_lines)              ; add polyline
;   dum2 = gsn_add_text(wks,bot_plot,"Anthropogenic + Natural",1910,0.705,res_text); add text
  
;   yy = (/0.79,0.79/)
;   res_lines@gsLineColor       = "blue"                                 ; change to blue
;   dum3 = gsn_add_polyline(wks,bot_plot,xx,yy,res_lines)                ; add polyline
;   dum4 = gsn_add_text(wks,bot_plot,"Natural",1910,0.79,res_text)       ; add text
  
;   yy = (/0.875,0.875/)
;   res_lines@gsLineColor       = "black"                                ; change to black
;   dum5 = gsn_add_polyline(wks,bot_plot,xx,yy,res_lines)                ; add polyline
;   dum6 = gsn_add_text(wks,bot_plot,"Observations",1910,0.875,res_text) ; add text

; lnres = True
; lnres@gsLineColors = (/"black","blue","red","green4"/)

; x1pos = (/0.1, 0.25, 0.4, 0.55, 0.7, 0.85/)
; x2pos = (/0.15, 0.3, 0.45, 0.6, 0.75, 0.9/)
; gsn_polyline_ndc(wks, x1pos, x2pos, lnres)


pres = True
pres@gsnMaximize = True
gsn_panel(wks, plot, (/3,2/), pres)


end