
begin

yrStrt = 1979
yrLast = 2015
year = ispan(yrStrt, yrLast, 1)
Month = ispan(1, 12, 1)

;=================================================================
filitp = "~/wuyt/data/topo.sfc.1x1.nc"
ftp = addfile(filitp, "r")
topo = ftp->topo(::-1,:)
printVarSummary(topo)
;=================================================================
casename = (/"NOQ1","CTRL"/)
index_wy = new((/3,12/), float)

do i = 0, dimsizes(casename)-1
    diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(i)+"_F2000_F19G16_CAM4/pro/"
    filiu = "F2000_"+casename(i)+".plevel.U.mon.0601-3012.nc"
    filiv = "F2000_"+casename(i)+".plevel.V.mon.0601-3012.nc"

    fu = addfile(diri+filiu, "r")
    fv = addfile(diri+filiv, "r")
    U = fu->U
    V = fv->V
    u = clmMonTLLL(U)
    v = clmMonTLLL(V)
    ;===============================================================================
    ; Webster-Yang index (Asian Monsoon)
    ;(U850-U200 averaged over 0-20N, 40E-110E)
    ;===============================================================================
    u_wy = wgt_areaave_Wrap(u(:,:,{0:20},{40:110}), 1, 1, 0)
    printVarSummary(u_wy)
    index_wy(i,:) = u_wy(:,{850}) - u_wy(:,{200})
    copy_VarCoords(u_wy(:,0), index_wy(i,:))
delete([/diri,filiu,filiv,fu,fv,U,V,u,v,u_wy/])
end do
printVarSummary(index_wy)

;=================================================================

diri = "/home/yangsong3/data-observation/ERA-interim/monthly/"
filiu = "uwnd.mon.mean.nc"
filiv = "vwnd.mon.mean.nc"

fu = addfile(diri+filiu, "r")
fv = addfile(diri+filiv, "r")
TIME = fu->time
YYYY = cd_calendar(TIME, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

U = short2flt(fu->u(iYYYY,:,:,:))
V = short2flt(fv->v(iYYYY,:,:,:))
u = clmMonTLLL(U)
v = clmMonTLLL(V)
;===============================================================================
; Webster-Yang index (Asian Monsoon)
;(U850-U200 averaged over 0-20N, 40E-110E)
;===============================================================================
u_wy = wgt_areaave_Wrap(u(:,:,{0:20},{40:110}), 1, 1, 0)
printVarSummary(u_wy)
index_wy(2,:) = u_wy(:,{850}) - u_wy(:,{200})
copy_VarCoords(u_wy(:,0), index_wy(2,:))
printVarSummary(index_wy)
index_wy = index_wy(::-1,:)
dk = index_wy
do i = 1, 11
    dk(:,i) = index_wy(:,i) - index_wy(:,i-1)
end do
dk(:,0) = index_wy(:,0) - index_wy(:,11)

; rc1 = regline(Month(2:5), index_wy(:,2:5))
; rc2 = regline(Month(7:10), index_wy(:,7:10))
; print("rc1 "+rc1+"  rc2 "+rc2) ;
;================================================================
; plot
;================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_Spr_Aut_heating/"
figname = "draw_monthly_WY_index_OBS_CTRL_NOQ1"

wks = gsn_open_wks("pdf", figpath+figname)
plot = new(1, graphic)

res = True
res@gsnDraw  = False
res@gsnFrame = False
res@gsnMaximize = True

res@vpWidthF = 0.7
res@vpHeightF = 0.4

res@trXMinF = 1
res@trXMaxF = 12
res@tmXBMode = "Explicit"
res@tmXBValues = Month
res@tmXBLabels = Month
res@tmXTOn = False
res@tiXAxisString = "Month"

; res@trYMinF = -16
; res@trYMaxF = 20
; res@tmYLMode = "Explicit"
; res@tmYLValues = ispan(-15, 20, 5)
; res@tmYLLabels = ispan(-15, 20, 5)
res@tmYROn = False

res@gsnYRefLine = 0
res@gsnYRefLineThicknessF = 5
res@gsnYRefLineColor = "grey53"
res@gsnYRefLineDashPattern = 2

;Colors = (/"deeppink2", "blue", "purple","green4","orange2", "black" /);
Colors = (/"black","blue","orange2"/)
Patterns = (/0,0,0/)
res@xyLineThicknessF = 5
res@xyLineColors = Colors
res@xyDashPatterns = Patterns
;res@gsnLeftString = "Webster-Yang Index, "+sprintf("%4.2f", rc1)+" / "+sprintf("%4.2f", rc2)
res@gsnLeftString = "Webster-Yang Index"
res@gsnStringFontHeightF = 0.025

plot(0) = gsn_csm_xy(wks, Month, index_wy, res)

;==================================================================================
    res_lines                   = True                  ; polyline mods desired
    res_lines@gsLineThicknessF  = 5                     ; line thicker
    

    res_text                    = True                  ; text mods desired
    res_text@txFontHeightF      = 0.018                 ; change text size
    res_text@txJust             = "CenterLeft"          ; text justification

    ;Texts  = (/"CTRL ("+sprintf("%4.2f", rc11)+"/"+sprintf("%4.2f", rc12)+")",\
    ;          "NOQ1 ("+sprintf("%4.2f", rc21)+"/"+sprintf("%4.2f", rc22)+")"/)
    Texts = (/"OBS","CTRL","NOQ1"/)
;============================================================
    gx = (/1.5,2.2/)
    text_pos = 2.3
    gy1 = (/26,26/)
    gy2 = (/23,23/)
    gy3 = (/20,20/)
    gy4 = (/9,9/)
    ; gy5 = (/300,300/)
    ; gy6 = (/280,280/)

    res_lines@gsLineColor = Colors(0)
    res_lines@gsLineDashPattern = 0
    dum11 = gsn_add_polyline(wks, plot(0), gx, gy1, res_lines)
    dum12 = gsn_add_text(wks, plot(0), Texts(0), text_pos, gy1, res_text)
    res_lines@gsLineColor = Colors(1)
    ;res_lines@gsLineDashPattern = 2
    dum21 = gsn_add_polyline(wks, plot(0), gx, gy2, res_lines)
    dum22 = gsn_add_text(wks, plot(0), Texts(1), text_pos, gy2, res_text)
    res_lines@gsLineColor = Colors(2)
    ;res_lines@gsLineDashPattern = 2
    dum31 = gsn_add_polyline(wks, plot(0), gx, gy3, res_lines)
    dum32 = gsn_add_text(wks, plot(0), Texts(2), text_pos, gy3, res_text)


draw(plot(0))
frame(wks)


end