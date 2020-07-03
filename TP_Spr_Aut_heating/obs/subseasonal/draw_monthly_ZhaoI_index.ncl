
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
index_zhao = u2 - u1
copy_VarCoords(u1, index_zhao)
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
    dk(i) = index_zhao(i) - index_zhao(i-1)
end do
dk(0) = index_zhao(0) - index_zhao(11)

rc1 = regline(Month(2:5), index_zhao(2:5))
rc2 = regline(Month(7:10), index_zhao(7:10))
print("rc1 "+rc1+"  rc2 "+rc2) ;
;================================================================
; plot
;================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_Spr_Aut_heating/obs/subseasonal/"
figname = "draw_monthly_ZhaoI_index"

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
Colors = (/"black","blue","forestgreen"/)
Patterns = (/0,0,0/)
res@xyLineThicknessF = 5
res@xyLineColors = Colors
res@xyDashPatterns = Patterns
res@gsnLeftString = "Zhao Index, "+sprintf("%4.2f", rc1)+" / "+sprintf("%4.2f", rc2)
res@gsnStringFontHeightF = 0.025

plot(0) = gsn_csm_xy(wks, Month, index_zhao, res)
draw(plot(0))
frame(wks)


end