
begin

diri = "./"
casename = (/"CTRL","NOQ1","NOTP_NOQ1","NOTP"/) ;
ncase = dimsizes(casename)
indexname = (/"Chen","Yang","LYang","Jhun","Sun"/)
nindex = dimsizes(indexname)
month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
nm = dimsizes(month)

index = new((/ncase,nm,nindex/), float)
do i = 0, ncase-1
    fili = "calc_monsoon_index_12month_"+casename(i)+".txt"
    index(i,:,:) = asciiread(diri+fili, (/nm,nindex/), "float")
end do
printVarSummary(index)

title = (/"Thermal","Dynamical"/)
pp = index(0:1,:,:)
pp(0,:,:) = index(0,:,:) - index(1,:,:)
pp(1,:,:) = index(1,:,:) - index(2,:,:)
pp!0 = "effect"
pp&effect = title

djf_mean = pp(:,0,:)
djf_mean = (pp(:,11,:) + pp(:,0,:) + pp(:,1,:))/3


filipath = "./"
filiname = "winter_monsoon_index_TD.txt"

write_table(filipath+filiname, "w", [/"Chen","Yang","LYang","Jhun","Sun"/], "%10s %10s %10s %10s %10s")
write_table(filipath+filiname, "a", [/djf_mean(:,0),djf_mean(:,1),djf_mean(:,2),djf_mean(:,3),djf_mean(:,4)/], "%10.2f %10.2f %10.2f %10.2f %10.2f")

indexx = index
ppx = pp
monthx = month
do i = 6, 11
    ppx(:,i,:) = pp(:,i-6,:)
    indexx(:,i,:) = index(:,i-6,:)
    monthx(i) = month(i-6)
end do
do i = 0, 5
    ppx(:,i,:) = pp(:,i+6,:)
    indexx(:,i,:) = index(:,i+6,:)
    monthx(i) = month(i+6)
end do

;==================================================================
;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20200102/"
figname = "winter_monsoon_index_evolution"
xaxis = ispan(1, 12, 1)
Colors = (/"black","blue","red","green4"/)
Texts = (/"Ctrl","NH","NTNH","NT"/)

wks = gsn_open_wks("eps", figpath+figname)
nplot = 2*nindex
plot = new(nplot, graphic)

res = True
res@gsnDraw  = False
res@gsnFrame = False

res@vpHeightF= 0.4                    ; change aspect ratio of plot
res@vpWidthF = 0.7     
res@gsnStringFontHeightF = 0.035
res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025

res@trXMinF = 4
res@trXMaxF = 9
res@xyDashPattern = 0
res@xyLineThicknessF = 5
res@xyLineColors = (/"black","blue","red","green4"/)
res@tmXBMode = "Explicit"
res@tmXBValues = xaxis
res@tmXBLabels = monthx
res@gsnYRefLine = 0
res@gsnYRefLineColor = "darkgoldenrod1"
res@gsnYRefLineThicknessF = 3
res@gsnYRefLineDashPattern = 2

str = (/"(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)","(j)"/)
do i = 0, nindex-1
    delete(res@gsnYRefLine)
    res@gsnLeftString = str(i)+" "+indexname(i)
    plot(i) = gsn_csm_xy(wks, xaxis(3:8), indexx(:,3:8,i), res)

    res@gsnYRefLine = 0
    res@gsnLeftString = str(i+nindex)+" "+indexname(i)
    plot(i+nindex) = gsn_csm_xy(wks, xaxis(3:8), ppx(:,3:8,i), res)
end do
;==================================================================================
    res_lines                   = True                  ; polyline mods desired
    res_lines@gsLineThicknessF  = 5                     ; line thicker
    

    res_text                    = True                  ; text mods desired
    res_text@txFontHeightF      = 0.02                 ; change text size
    res_text@txJust             = "CenterLeft"          ; text justification

; ;============================================================
    interval = 0.6
    gx = (/5.5,6.2/)
    text_pos = 6.4
    gy1 = (/4.5,4.5/)
    gy2 = gy1-interval ;(/15,15/)
    gy3 = gy2-interval
    gy4 = gy3-interval
    print(gy1)

    res_lines@gsLineColor = Colors(0)
    ;res_lines@gsLineDashPattern = Patterns(0)
    dum11 = gsn_add_polyline(wks, plot(0), gx, gy1, res_lines)
    dum12 = gsn_add_text(wks, plot(0), Texts(0), text_pos, gy1, res_text)
    res_lines@gsLineColor = Colors(1)
    ;res_lines@gsLineDashPattern = Patterns(1)
    dum21 = gsn_add_polyline(wks, plot(0), gx, gy2, res_lines)
    dum22 = gsn_add_text(wks, plot(0), Texts(1), text_pos, gy2, res_text)
    res_lines@gsLineColor = Colors(2)
    ;res_lines@gsLineDashPattern = Patterns(2)
    dum31 = gsn_add_polyline(wks, plot(0), gx, gy3, res_lines)
    dum32 = gsn_add_text(wks, plot(0), Texts(2), text_pos, gy3, res_text)
    res_lines@gsLineColor = Colors(3)
    ;res_lines@gsLineDashPattern = Patterns(3)
    dum41 = gsn_add_polyline(wks, plot(0), gx, gy4, res_lines)
    dum42 = gsn_add_text(wks, plot(0), Texts(3), text_pos, gy4, res_text) 

    interval = 0.2
    gx = (/4.1,4.8/)
    text_pos = 5
    gy1 = (/-0.3,-0.3/)
    gy2 = gy1-interval ;(/15,15/)
    res_lines@gsLineColor = Colors(0)
    ;res_lines@gsLineDashPattern = Patterns(4)
    dum51 = gsn_add_polyline(wks, plot(nindex), gx, gy1, res_lines)
    dum52 = gsn_add_text(wks, plot(nindex), "Thermal", text_pos, gy1, res_text)
    res_lines@gsLineColor = Colors(1)
    ;res_lines@gsLineDashPattern = Patterns(5)
    dum61 = gsn_add_polyline(wks, plot(nindex), gx, gy2, res_lines)
    dum62 = gsn_add_text(wks, plot(nindex), "Dynamical", text_pos, gy2, res_text)

pres = True
pres@gsnMaximize = True
gsn_panel(wks, plot, (/2,nindex/), pres)

end