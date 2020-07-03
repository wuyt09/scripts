
begin

yrStrt = 1979
yrLast = 2015
year = ispan(yrStrt, yrLast, 1)
Month = ispan(1, 12, 1)
LEV = 500

latstr = 10
latend = 90
lonstr = 110
lonend = 150
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

U = short2flt(fu->u(iYYYY,{LEV},{latstr:latend},{lonstr:lonend}))
V = short2flt(fv->v(iYYYY,{LEV},{latstr:latend},{lonstr:lonend}))
u = clmMonTLL(U)
v = clmMonTLL(V)

filiz = "hgt.mon.mean.nc"
fz = addfile(diri+filiz, "r")
Z = short2flt(fz->z(iYYYY,{LEV},{latstr:latend},{lonstr:lonend}))
Z = Z/9.8
z = clmMonTLL(Z)
zx = ndtooned(z)
ind_height = ind_resolve(ind(zx.ge.5880), dimsizes(z))
printVarSummary(ind_height)
exit

lon = tofloat(u&longitude)
lat = tofloat(u&latitude)
dlat = (lat(2)-lat(1))*0.0174533
dy = 6378388.*dlat
dudx = u
dudx = center_finite_diff_n(u, dy, True, 0, 1)
printVarSummary(dudx)

UU = ndtooned(u)
XX = ndtooned(dudx)
error = 0.5
index1 = ind_resolve(ind(XX.ge.0), dimsizes(dudx)) 
index2 = ind_resolve(ind(abs(UU).le.error), dimsizes(u)) 
ind_final = ind(index1.eq.True .and. index2.eq.True .and. ind_height.eq.True)
print(ind_final)
;print(index2)
exit
;===============================================================================
; Webster-Yang index (Asian Monsoon)
;(U850-U200 averaged over 0-20N, 40E-110E)
;===============================================================================
u_wy = wgt_areaave_Wrap(u(:,:,{0:20},{40:110}), 1, 1, 0)
printVarSummary(u_wy)
index_wy = u_wy(:,{850}) - u_wy(:,{200})
copy_VarCoords(u_wy(:,0), index_wy)
printVarSummary(index_wy)

dk = index_wy
do i = 1, 11
    dk(i) = index_wy(i) - index_wy(i-1)
end do
dk(0) = index_wy(0) - index_wy(11)

rc1 = regline(Month(2:5), index_wy(2:5))
rc2 = regline(Month(7:10), index_wy(7:10))
print("rc1 "+rc1+"  rc2 "+rc2) ;
;================================================================
; plot
;================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_Spr_Aut_heating/obs/"
figname = "draw_monthly_WY_index"

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
res@gsnLeftString = "Webster-Yang Index, "+sprintf("%4.2f", rc1)+" / "+sprintf("%4.2f", rc2)
res@gsnStringFontHeightF = 0.025

plot(0) = gsn_csm_xy(wks, Month, index_wy, res)
draw(plot(0))
frame(wks)


end