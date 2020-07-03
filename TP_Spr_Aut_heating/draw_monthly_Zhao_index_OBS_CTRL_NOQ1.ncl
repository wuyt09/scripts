
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
index_zhao = new((/3,12/), float)

do i = 0, dimsizes(casename)-1
    diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(i)+"_F2000_F19G16_CAM4/pro/"
    fili = "F2000_"+casename(i)+".PSL.mon.0601-3012.nc"

    fu = addfile(diri+fili, "r")
    U = fu->PSL
    u = clmMonTLL(U)
    ;===============================================================================
    ; Zhao&Zhou monsoon index (For east asian monsoon)
    ; (SLP (30N-40N, 160E) - (SLP 40N-50N, 110E))
    ;===============================================================================
    u_zhao = u
    u1 = dim_avg_n_Wrap(u_zhao(:,{40:50},{110}), 1)
    u2 = dim_avg_n_Wrap(u_zhao(:,{25:35},{160}), 1)
    u1 = dim_standardize_Wrap(u1, 0)
    u2 = dim_standardize_Wrap(u2, 0)
    index_zhao(i,:) = u2 - u1
    copy_VarCoords(u1, index_zhao(i,:))
delete([/diri,fili,fu,U,u,u1,u2,u_zhao/])
end do

diri = "/home/yangsong3/data-observation/ERA-interim/monthly/"
fili = "slp.mon.mean.nc"

fu = addfile(diri+fili, "r")
TIME = fu->time
YYYY = cd_calendar(TIME, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

U = short2flt(fu->msl(iYYYY,:,:))
;U = dim_standardize_n_Wrap(U, 0, 0)
u = clmMonTLL(U)

;===============================================================================
; Zhao&Zhou monsoon index (For east asian monsoon)
; (SLP (30N-40N, 160E) - (SLP 40N-50N, 110E))
;===============================================================================
u_zhao = u
u1 = dim_avg_n_Wrap(u_zhao(:,{40:50},{110}), 1)
u2 = dim_avg_n_Wrap(u_zhao(:,{25:35},{160}), 1)
u1 = dim_standardize_Wrap(u1, 0)
u2 = dim_standardize_Wrap(u2, 0)
index_zhao(2,:) = u2 - u1
copy_VarCoords(u1, index_zhao(2,:))
index_zhao = index_zhao(::-1,:)
printVarSummary(index_zhao)

; ;===============================================================================
; ; Shi monsoon index (For east asian monsoon)
; ; (SLP (30N-40N, 160E) - (SLP 40N-50N, 110E))
; ;===============================================================================
; u_zhao = u
; u1 = dim_avg_n_Wrap(u_zhao(:,{20:50},{110}), 1)
; u2 = dim_avg_n_Wrap(u_zhao(:,{20:50},{160}), 1)
; u1 = dim_standardize_Wrap(u1, 0)
; u2 = dim_standardize_Wrap(u2, 0)
; index_zhao = u1 - u2
; copy_VarCoords(u1, index_zhao)
; printVarSummary(index_zhao)


dk = index_zhao
do i = 1, 11
    dk(:,i) = index_zhao(:,i) - index_zhao(:,i-1)
end do
dk(:,0) = index_zhao(:,0) - index_zhao(:,11)

rc01 = regline(Month(2:5), index_zhao(0,2:5))
rc02 = regline(Month(7:10), index_zhao(0,7:10))
print("OBS: rc1 "+rc01+"  rc2 "+rc02) ;

rc11 = regline(Month(2:5), index_zhao(1,2:5))
rc12 = regline(Month(7:10), index_zhao(1,7:10))
print("CTRL: rc1 "+rc11+"  rc2 "+rc12) ;

rc21 = regline(Month(2:5), index_zhao(2,2:5))
rc22 = regline(Month(7:10), index_zhao(2,7:10))
print("NOQ1: rc1 "+rc21+"  rc2 "+rc22) ;
;================================================================
; plot
;================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_Spr_Aut_heating/"
figname = "draw_monthly_ZhaoI_index_OBS_CTRL_NOQ1"

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
res@gsnLeftString = "Zhao Index";, +sprintf("%4.2f", rc1)+" / "+sprintf("%4.2f", rc2)
res@gsnStringFontHeightF = 0.025

plot(0) = gsn_csm_xy(wks, Month, index_zhao, res)
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
    gy1 = (/2.6,2.6/)
    gy2 = (/2.3,2.3/)
    gy3 = (/2,2/)
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