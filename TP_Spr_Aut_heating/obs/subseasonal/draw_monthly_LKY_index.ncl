
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
; LKY monsoon index (For east asian monsoon) RM2
; (U200 (40N-50N, 110-150E) - (U200 25N-35N, 110E-150E))
;===============================================================================
u_rm2 = u(:,{200},:,:)
u1_rm2 = wgt_areaave_Wrap(u_rm2(:,{40:50},{110:150}), 1, 1, 0)
u2_rm2 = wgt_areaave_Wrap(u_rm2(:,{25:35},{110:150}), 1, 1, 0)
index_rm2 = u1_rm2 - u2_rm2
copy_VarCoords(u1_rm2, index_rm2)
printVarSummary(index_rm2)


dk = index_rm2
do i = 1, 11
    dk(i) = index_rm2(i) - index_rm2(i-1)
end do
dk(0) = index_rm2(0) - index_rm2(11)

rc1 = regline(Month(2:5), index_rm2(2:5))
rc2 = regline(Month(7:10), index_rm2(7:10))
print("rc1 "+rc1+"  rc2 "+rc2) ;
;================================================================
; plot
;================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_Spr_Aut_heating/obs/"
figname = "draw_monthly_LKY_index"

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
res@gsnLeftString = "LKY Index, "+sprintf("%4.2f", rc1)+" / "+sprintf("%4.2f", rc2)
res@gsnStringFontHeightF = 0.025

plot(0) = gsn_csm_xy(wks, Month, index_rm2, res)
draw(plot(0))
frame(wks)


end