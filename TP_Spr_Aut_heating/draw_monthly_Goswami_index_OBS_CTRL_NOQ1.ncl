
begin

yrStrt = 1979
yrLast = 2015
year = ispan(yrStrt, yrLast, 1)
Month = ispan(1, 12, 1)

;=================================================================
filitp = "~/wuyt/data/topo.sfc.1x1.nc"
ftp = addfile(filitp, "r")
topo = ftp->topo
printVarSummary(topo)
;=================================================================
Month = ispan(1, 12, 1)
casename = (/"NOQ1","CTRL"/)
index_rm1 = new((/3,12/), float)

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
    ; South Asian monsoon index RM1 Lau Kim Yang (2000)/Goswami
    ; (V850-V200 averaged over 10N-30N, 70E-110E)
    ;===============================================================================
    v_rm1 = wgt_areaave_Wrap(v(:,:,{10:30},{70:110}), 1, 1, 0)
    index_rm1(i,:) = v_rm1(:,{850}) - v_rm1(:,{200})
    copy_VarCoords(v_rm1(:,0), index_rm1(i,:))
delete([/diri,filiu,filiv,fu,fv,U,V,u,v,v_rm1/])
end do

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
; South Asian monsoon index RM1 Lau Kim Yang (2000)/Goswami
; (V850-V200 averaged over 10N-30N, 70E-110E)
;===============================================================================
v_rm1 = wgt_areaave_Wrap(v(:,:,{10:30},{70:110}), 1, 1, 0)
index_rm1(2,:) = v_rm1(:,{850}) - v_rm1(:,{200})
copy_VarCoords(v_rm1(:,0), index_rm1(2,:))
index_rm1 = index_rm1(::-1,:)

dk = index_rm1
do i = 1, 11
    dk(:,i) = index_rm1(:,i) - index_rm1(:,i-1)
end do
dk(:,0) = index_rm1(:,0) - index_rm1(:,11)

rc01 = regline(Month(2:5), index_rm1(0,2:5))
rc02 = regline(Month(7:10), index_rm1(0,7:10))
print("OBS: rc1 "+rc01+"  rc2 "+rc02) ;

rc11 = regline(Month(2:5), index_rm1(1,2:5))
rc12 = regline(Month(7:10), index_rm1(1,7:10))
print("CTRL: rc1 "+rc11+"  rc2 "+rc12) ;

rc21 = regline(Month(2:5), index_rm1(2,2:5))
rc22 = regline(Month(7:10), index_rm1(2,7:10))
print("NOQ1: rc1 "+rc21+"  rc2 "+rc22) ;
;================================================================
; plot
;================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_Spr_Aut_heating/"
figname = "draw_monthly_Goswami_index_OBS_CTRL_NOQ1"

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
res@gsnLeftString = "Goswami Index";, +sprintf("%4.2f", rc1)+" / "+sprintf("%4.2f", rc2)
res@gsnStringFontHeightF = 0.025

plot(0) = gsn_csm_xy(wks, Month, index_rm1, res)

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
    gy1 = (/5,5/)
    gy2 = (/4,4/)
    gy3 = (/3,3/)
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