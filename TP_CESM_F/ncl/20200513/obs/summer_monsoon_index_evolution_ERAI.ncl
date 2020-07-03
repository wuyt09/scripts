
begin

diri = "./"
casename = (/"high_years","low_years"/) 
ncase = dimsizes(casename)
indexname = (/"WY","SA-Goswami","SA-Wang","SEA-Wang","EA-LKY","EA-Wang"/)
indexname2 = (/"WY","SA-Goswami","SEA-Wang","EA-LKY","EA-Wang"/)
nindex = dimsizes(indexname2)
nindexx = dimsizes(indexname)
month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
nm = dimsizes(month)

indexx = new((/ncase,nm,nindexx/), float)
do i = 0, ncase-1
    fili = "calc_monsoon_index_12month_ERAI_JJA_"+casename(i)+".txt"
    indexx(i,:,:) = asciiread(diri+fili, (/nm,nindexx/), "float")
end do
index = indexx(:,:,0:nindexx-2)
index(:,:,0:1) = indexx(:,:,0:1)
index(:,:,2:nindex-1) = indexx(:,:,3:nindexx-1)
index!0 = "case"
index&case = casename
index!1 = "month"
index&month = month
index!2 = "indexname"
index&indexname = indexname2
printVarSummary(index)

; x0 = dim_standardize_n_Wrap(index, 1, 1)
; x = new((/ncase+1,nm,nindex/), float)
; x(0:ncase-1,:,:) = x0
; x(ncase,:,:) = x0(0,:,:) - x0(1,:,:)
; copy_VarCoords(x0(0,:,:), x(0,:,:))
; printVarSummary(x)

;==================================================================
;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20200513/obs/"
figname = "summer_monsoon_index_evolution_ERAI_JJA_Comp"
xaxis = ispan(1, 12, 1)
Colors = (/"red","blue","black"/) ;,"green4","purple"
Patterns = (/0,0,0/)
Texts = (/"High_TTT","Low_TTT"/) ;,"NTNH","NT"

wks = gsn_open_wks("pdf", figpath+figname)
nplot = 6
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
res@trXMaxF = 10
; res@trYMinF = -1
; res@trYMaxF = 2
res@xyDashPatterns = Patterns
res@xyLineThicknessF = 5
res@xyLineColors = Colors
res@tmXBMode = "Explicit"
res@tmXBValues = xaxis
res@tmXBLabels = month
res@gsnYRefLine = 0
res@gsnYRefLineColor = "black" ;darkgoldenrod1
res@gsnYRefLineThicknessF = 3
res@gsnYRefLineDashPattern = 2
res@gsnYRefLine = 0

str = (/"(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)","(j)"/)
; do i = 0, nindex-1
;     res@gsnLeftString = str(i)+" "+indexname2(i)
;     plot(i) = gsn_csm_xy(wks, xaxis(3:8), index(:,3:8,i), res)
; end do
do i = 0, nindexx-1
    res@gsnLeftString = str(i)+" "+indexname(i)
    plot(i) = gsn_csm_xy(wks, xaxis(3:9), indexx(:,3:9,i), res)
end do


;==================================================================================
    res_lines                   = True                  ; polyline mods desired
    res_lines@gsLineThicknessF  = 3                     ; line thicker
    

    res_text                    = True                  ; text mods desired
    res_text@txFontHeightF      = 0.018                 ; change text size
    res_text@txJust             = "CenterLeft"          ; text justification

; ;============================================================
    interval = 5 ;0.3
    gx = (/4.2,4.7/)
    text_pos = 4.8
    ; gy1 = (/1.7,1.7/)
    gy1 = (/26,26/)
    gy2 = gy1-interval ;(/15,15/)
    gy3 = gy2-interval
    ; gy4 = gy3-interval
    print(gy1)

    res_lines@gsLineColor = Colors(0)
    res_lines@gsLineDashPattern = Patterns(0)
    dum11 = gsn_add_polyline(wks, plot(0), gx, gy1, res_lines)
    dum12 = gsn_add_text(wks, plot(0), Texts(0), text_pos, gy1, res_text)
    res_lines@gsLineColor = Colors(1)
    res_lines@gsLineDashPattern = Patterns(1)
    dum21 = gsn_add_polyline(wks, plot(0), gx, gy2, res_lines)
    dum22 = gsn_add_text(wks, plot(0), Texts(1), text_pos, gy2, res_text)
    ; res_lines@gsLineColor = Colors(2)
    ; res_lines@gsLineDashPattern = Patterns(2)
    ; dum31 = gsn_add_polyline(wks, plot(0), gx, gy3, res_lines)
    ; dum32 = gsn_add_text(wks, plot(0), Texts(2), text_pos, gy3, res_text)
    ; res_lines@gsLineColor = Colors(3)
    ; ;res_lines@gsLineDashPattern = Patterns(3)
    ; dum41 = gsn_add_polyline(wks, plot(0), gx, gy4, res_lines)
    ; dum42 = gsn_add_text(wks, plot(0), Texts(3), text_pos, gy4, res_text) 

    ; interval = 2
    ; gx = (/7,7.7/)
    ; text_pos = 7.9
    ; gy1 = (/8,8/)
    ; gy2 = gy1-interval ;(/15,15/)
    ; res_lines@gsLineColor = Colors(0)
    ; ;res_lines@gsLineDashPattern = Patterns(4)
    ; dum51 = gsn_add_polyline(wks, plot(nindex), gx, gy1, res_lines)
    ; dum52 = gsn_add_text(wks, plot(nindex), "Thermal", text_pos, gy1, res_text)
    ; res_lines@gsLineColor = Colors(1)
    ; ;res_lines@gsLineDashPattern = Patterns(5)
    ; dum61 = gsn_add_polyline(wks, plot(nindex), gx, gy2, res_lines)
    ; dum62 = gsn_add_text(wks, plot(nindex), "Dynamical", text_pos, gy2, res_text)

pres = True
pres@gsnMaximize = True
gsn_panel(wks, plot, (/2,3/), pres)

end