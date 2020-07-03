begin

indexname = (/"WY","WNP","SA","SA-II","EA","AUS"/)
casename = (/"CTRL","NOTP","NOQ1","NOTP_NOQ1"/)

diri = "./"
fili1 = "calc_monsoon_index_12month_CTRL.txt"
fili2 = "calc_monsoon_index_12month_NOTP.txt"
fili3 = "calc_monsoon_index_12month_NOQ1.txt"
fili4 = "calc_monsoon_index_12month_NOTP_NOQ1.txt"

ctrl = asciiread(diri+fili1, (/12,6/), "float")
notp = asciiread(diri+fili2, (/12,6/), "float")
noq1 = asciiread(diri+fili3, (/12,6/), "float")
noei = asciiread(diri+fili4, (/12,6/), "float")
printVarSummary(ctrl)

index = new((/4,5,6/), float)
index(0,:,:) = ctrl(4:8,:)
index(1,:,:) = notp(4:8,:)
index(2,:,:) = noq1(4:8,:)
index(3,:,:) = noei(4:8,:)
printVarSummary(index)

index2 = new((/4,6,6/), float)
index2(:,0:4,:) = index
index2(:,5,:) = dim_avg_n_Wrap(index, 1)
printVarSummary(index2)
;==================================================================
;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/"
figname = "draw_monsoon_index_MJJAS_CTRL_S1_S2"
xaxis = fspan(1, 6, 6)
xtitle = (/"May","June","July","Aug","Sep","Ave"/)

wks = gsn_open_wks("pdf", figpath+figname)
nplot = dimsizes(indexname)

;==========================================================
;****************************************************
; plot
;****************************************************

plots = new(4, graphic)
plotb = new(4, graphic)
plotq = new(4, graphic)
plotc = new(4, graphic)

res   = True
res@gsnDraw  = False
res@gsnFrame = False
;res@gsnMaximize = True
; these four resources allow the user to stretch the plot size, and
; decide exactly where on the page to draw it.

; res@vpXF      = 0.10        ; In page coordinates, where to start
; res@vpYF      = 0.75        ; the plot
res@vpHeightF = 0.35        ; Changes the aspect ratio
res@vpWidthF  = 0.7

;res@tiYAxisString = "Contribution (~S~o~N~C)" ; y-axis label      
;res@tiXAxisString = ""
;res@tiMainString  = "****deva*****" ; title

res@gsnYRefLine           = 0.              ; reference line   
res@gsnXYBarChart         = True            ; create bar chart 
res@gsnXYBarChartBarWidth = 0.25
;res@gsnAboveYRefLineColor = "red" ;firebrick2           ; above ref line fill red
;res@gsnBelowYRefLineColor = "blue" ;dodgerblue2         ; below ref line fill blue

res@tmXBMode   = "Explicit"
res@tmXBValues = xaxis
print(res@tmXBValues)
res@tmXBLabels = xtitle

res@tmYLLabelFontHeightF = 0.0205
;res@tmXBLabelFontHeightF = 0.0205
res@tmXBLabelFontHeightF = 0.0205
res@gsnLeftStringFontHeightF = 0.025
res@tmXTOn = False
res@tmYROn = False

res@trXMinF               = 0.4
res@trXMaxF               = 6.5
;res@trYMinF               = -30
;res@trYMaxF               = 30.0
;res@tmYLValues = fspan(-0.5, 0.4, 10)
;res@tmYLLabels = (/-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4/)
       
       res@trYMinF               = 0
       res@trYMaxF               = 30
       res@gsnLeftString = indexname(0)
       res@gsnXYBarChartColors = "green4" ;powderblue
       plots(0) = gsn_csm_xy(wks,fspan(.75,5.75,6),index2(0,:,0),res)

       res@gsnXYBarChartColors = "red2"
       plots(1) = gsn_csm_xy(wks,fspan(1.0,6,6),index2(2,:,0),res)

       res@gsnXYBarChartColors = "dodgerblue"
       plots(2) = gsn_csm_xy(wks,fspan(1.25,6.25,6),index2(3,:,0),res)

       ;res@gsnXYBarChartColors = "darkgoldenrod4"
       ;plots(3) = gsn_csm_xy(wks,fspan(1.5,5.5,5),index2(3,:,0),res)

       overlay(plots(0),plots(1))
       overlay(plots(0),plots(2))
       ;overlay(plots(0),plots(3))

       res@trYMinF               = -35
       res@trYMaxF               = 30
       res@gsnLeftString = indexname(4)
       res@gsnXYBarChartColors = "green4" ;powderblue
       plotb(0) = gsn_csm_xy(wks,fspan(.75,5.75,6),index2(0,:,4),res)

       res@gsnXYBarChartColors = "red2"
       plotb(1) = gsn_csm_xy(wks,fspan(1.0,6,6),index2(2,:,4),res)

       res@gsnXYBarChartColors = "dodgerblue"
       plotb(2) = gsn_csm_xy(wks,fspan(1.25,6.25,6),index2(3,:,4),res)

       ;res@gsnXYBarChartColors = "darkgoldenrod4"
       ;plots(3) = gsn_csm_xy(wks,fspan(1.5,5.5,5),index2(3,:,0),res)

       overlay(plotb(0),plotb(1))
       overlay(plotb(0),plotb(2))
       ;overlay(plots(0),plots(3))

       res@trYMinF               = 0
       res@trYMaxF               = 20
       res@gsnLeftString = indexname(2)
       res@gsnXYBarChartColors = "green4" ;powderblue
       plotq(0) = gsn_csm_xy(wks,fspan(.75,5.75,6),index2(0,:,2),res)

       res@gsnXYBarChartColors = "red2"
       plotq(1) = gsn_csm_xy(wks,fspan(1.0,6,6),index2(2,:,2),res)

       res@gsnXYBarChartColors = "dodgerblue"
       plotq(2) = gsn_csm_xy(wks,fspan(1.25,6.25,6),index2(3,:,2),res)

       ;res@gsnXYBarChartColors = "darkgoldenrod4"
       ;plots(3) = gsn_csm_xy(wks,fspan(1.5,5.5,5),index2(3,:,0),res)

       overlay(plotq(0),plotq(1))
       overlay(plotq(0),plotq(2))
       ;overlay(plots(0),plots(3))

       res@trYMinF               = -10
       res@trYMaxF               = 10
       res@gsnLeftString = indexname(3)
       res@gsnXYBarChartColors = "green4" ;powderblue
       plotc(0) = gsn_csm_xy(wks,fspan(.75,5.75,6),index2(0,:,3),res)

       res@gsnXYBarChartColors = "red2"
       plotc(1) = gsn_csm_xy(wks,fspan(1.0,6,6),index2(2,:,3),res)

       res@gsnXYBarChartColors = "dodgerblue"
       plotc(2) = gsn_csm_xy(wks,fspan(1.25,6.25,6),index2(3,:,3),res)

       ;res@gsnXYBarChartColors = "darkgoldenrod4"
       ;plots(3) = gsn_csm_xy(wks,fspan(1.5,5.5,5),index2(3,:,0),res)

       overlay(plotc(0),plotc(1))
       overlay(plotc(0),plotc(2))
       ;overlay(plots(0),plots(3))

       pres = True
       ;pres@gsnDraw  = False
       pres@gsnFrame = False
       pres@gsnMaximize = True
       gsn_panel(wks, (/plots(0), plotb(0), plotq(0), plotc(0)/), (/2,2/), pres)
; Create some labelbars for a legend
     lbres                    = True          ; labelbar only resources
     lbres@vpWidthF           = 0.1           ; labelbar width
     lbres@vpHeightF          = 0.2           ; labelbar height
     lbres@lbBoxMajorExtentF  = 0.15          ; puts space between color boxes
     lbres@lbMonoFillPattern  = True          ; Solid fill pattern
     lbres@lbLabelFontHeightF = 0.015         ; font height. default is small
     lbres@lbLabelJust        = "CenterLeft"  ; left justify labels
     lbres@lbPerimOn          = False

     colors = (/"green4", "red2", "dodgerblue"/) ;,"darkgoldenrod4"
     labels = (/"CTRL", "NOQ1","NOTP_NOQ1"/)  ;"NOTP" ,
     xpos   = (/0.33, 0.5, 0.67/)
     
     ; do i=0,2
     ;   lbres@lbFillColors = colors(i)
     ;   gsn_labelbar_ndc(wks,1,labels(i),xpos(i),0.12,lbres)
     ; end do 

     ;xpos   = (/0.2, 0.3, 0.4/)
     
     do i=0,2
       lbres@lbFillColors = colors(i)
       gsn_labelbar_ndc(wks,1,labels(i),xpos(i),0.25,lbres)
     end do 
     ;drawNDCGrid(plots)
     frame(wks)   ; Advance the frame.


end