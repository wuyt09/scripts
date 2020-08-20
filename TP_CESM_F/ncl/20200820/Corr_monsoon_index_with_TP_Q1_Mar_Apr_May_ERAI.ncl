
begin

lonstr = 40
lonend = 150
latstr = -15
latend = 50

yrStrt = 1979
yrLast = 2015
year = ispan(yrStrt, yrLast, 1)
nyear = dimsizes(year)
;season = (/"MAM","SON"/)
; season = (/"MAM","JJA","SON"/)

;================================================================
filitp = "~/wuyt/data/topo.sfc.1x1.nc"
ftp = addfile(filitp, "r")
TOPO = ftp->topo(::-1,:)
printVarSummary(TOPO)
;================================================================

diri = "/home/yangsong3/wuyt/data/ERA-Interim/monthly/"
fili = "q1.mon.mean.1979.2016.ERAI.nc"

f = addfile(diri+fili, "r")
TIME = f->time
YYYY = cd_calendar(TIME, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

T = short2flt(f->q1(iYYYY,:,::-1,:))
printVarSummary(T)
printMinMax(T, True)

diri2 = "/home/yangsong3/data-observation/ERA-interim/monthly/"
filip = "Surface_pressure.mon.mean.nc"
fp = addfile(diri2+filip, "r")
PS0 = short2flt(fp->sp(iYYYY,::-1,:))
PS = linint2_Wrap(PS0&longitude, PS0&latitude, PS0, True, T&lon, T&lat, 0)
; print(T&lat)


topo = linint2_Wrap(TOPO&lon, TOPO&lat, TOPO, True, T&lon, T&lat, 0)

nlat = dimsizes(T&lat)
nlon = dimsizes(T&lon)

;===========integration========================================
; levs = 600
; levt = 300

lev = T&lev
plev = lev*100.
psfc = PS
ptop = min(plev)

dP = dpres_plevel_Wrap(plev, psfc, ptop, 0)
dP!1 = "lev"
dP&lev = lev
printVarSummary(dP)

Tdp = T*dP
copy_VarCoords(T, Tdp)
TTT = dim_sum_n_Wrap(Tdp, 1)
copy_VarCoords(T(:,0,:,:), TTT)
printVarSummary(TTT)
printMinMax(TTT, True)

;===============================================================

; Tclm = clmMonTLL(TTT)
; Tano = calcMonAnomTLL(TTT, Tclm)
; ;print(Tano(20,{25:40},{75:105}))
; printVarSummary(Tano)


; t = new((/3,nyear,nlat,nlon/), float)
; t_Mar = TTT(2::12,:,:)
; t_Apr = TTT(3::12,:,:)
; ; t_MarApr = (t_Mar+t_Apr)/2
; t_May = TTT(4::12,:,:)
; t(0,:,:,:) = t_Mar
; t(1,:,:,:) = t_Apr
; t(2,:,:,:) = t_May ;t_MarApr
; copy_VarCoords(TTT(0,:,:), t(0,0,:,:))
; ; t = month_to_seasonN(TTT, season)
; ; printVarSummary(t)
t = mask(TTT, topo.ge.1500, True)
copy_VarCoords(TTT, t)

ttt = wgt_areaave_Wrap(t(:,{25:40},{75:105}), 1, 1, 0)
printVarSummary(ttt)

;===============================================================
diri3 = "/home/yangsong3/data-observation/ERA-interim/monthly/"
filiu = "uwnd.mon.mean.nc"
filiv = "vwnd.mon.mean.nc"

fu = addfile(diri3+filiu, "r")
fv = addfile(diri3+filiv, "r")

TIMEo = fu->time
YYYYo = cd_calendar(TIMEo, -1)/100
iYYYYo = ind(YYYYo.ge.yrStrt.and.YYYYo.le.yrLast)

u = short2flt(fu->u(iYYYYo,:,{latstr:latend},{lonstr:lonend}))
v = short2flt(fv->v(iYYYYo,:,{latstr:latend},{lonstr:lonend}))
printVarSummary(u)

;===============================================================================
; Webster-Yang index (Asian Monsoon) 
;(U850-U200 averaged over 0-20N, 40E-110E)
;===============================================================================
u_wy = wgt_areaave_Wrap(u(:,:,{0:20},{40:110}), 1, 1, 0)
printVarSummary(u_wy)
index_wy = u_wy(:,{850}) - u_wy(:,{200})
copy_VarCoords(u_wy(:,0), index_wy)
printVarSummary(index_wy)
;===============================================================================
; East Asian - Western North Pacific monsoon index (southeast monsoon)Wang etal (2001)
;(U850 (5N-15N, 90E-130E) - U850 (22.5N-32.5N, 110E-140E))
;===============================================================================
u_wnp = u(:,{850},:,:)
u1_wnp = wgt_areaave_Wrap(u_wnp(:,{5:15},{90:130}), 1, 1, 0)
u2_wnp = wgt_areaave_Wrap(u_wnp(:,{22.5:32.5},{110:140}), 1, 1, 0)
index_wnp = u1_wnp - u2_wnp
; index_wnp = (-1.) * index_wnp
copy_VarCoords(u1_wnp, index_wnp)
printVarSummary(index_wnp)

;===============================================================================
; Dynamic Indian monsoon index Wang etal (2001)
; (U850 (5N-15N, 40E-80E) - (U850 20N-30N, 70E-90E))
;===============================================================================
; u_sa = u(:,:,{850},:,:)
; u1_sa = wgt_areaave_Wrap(u_sa(:,:,{5:15},{40:80}), 1, 1, 0)
; u2_sa = wgt_areaave_Wrap(u_sa(:,:,{20:30},{70:90}), 1, 1, 0)
; index_sa = u1_sa - u2_sa
; copy_VarCoords(u1_sa, index_sa)
; printVarSummary(index_sa)

;===============================================================================
; South Asian monsoon index RM1 Lau Kim Yang (2000)
; (V850-V200 averaged over 10N-30N, 70E-110E)
;===============================================================================
v_rm1 = wgt_areaave_Wrap(v(:,:,{10:30},{70:110}), 1, 1, 0)
index_rm1 = v_rm1(:,{850}) - v_rm1(:,{200})
copy_VarCoords(v_rm1(:,0), index_rm1)
printVarSummary(index_rm1)

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

;===============================================================================
indexname = (/"WY","WNP","SA","EA"/)
indx = new((/4,444/), float)
indx(0,:) = index_wy
indx(1,:) = index_rm1
indx(2,:) = index_wnp
indx(3,:) = index_rm2
indx!0 = "index"
indx&index = indexname
printVarSummary(indx)


;===============================================================================
; wtite data to txt file
;===============================================================================
Month = (/"March","April","May"/)
monS = 5
filipath = "./"
filiname = "Corr_monsoon_index_with_TPQ1_ERAI_May.txt"
corr = new((/3,4/), float)
do i = 0, 2
    do j = 0, 3
        qq = ttt(i+2::12)
        pp = indx(j,monS-1::12)
        corr(i,j) = escorc(qq, pp)
    end do
end do
write_table(filipath+filiname, "w", [/"Month","WY","WNP","SA","EA"/], "%10s %10s %10s %10s %10s")
write_table(filipath+filiname, "a", [/Month,corr(:,0),corr(:,1),corr(:,2),corr(:,3)/], "%10s %10.2f %10.2f %10.2f %10.2f")

monS = 6
filipath = "./"
filiname = "Corr_monsoon_index_with_TPQ1_ERAI_June.txt"
corr = new((/3,4/), float)
do i = 0, 2
    do j = 0, 3
        qq = ttt(i+2::12)
        pp = indx(j,monS-1::12)
        corr(i,j) = escorc(qq, pp)
    end do
end do
write_table(filipath+filiname, "w", [/"Month","WY","WNP","SA","EA"/], "%10s %10s %10s %10s %10s")
write_table(filipath+filiname, "a", [/Month,corr(:,0),corr(:,1),corr(:,2),corr(:,3)/], "%10s %10.2f %10.2f %10.2f %10.2f")

exit

;==================================================================
;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20200820/"
figname = "calc_monsoon_index_multi_model"
xaxis = ispan(6, 30, 1)
wks = gsn_open_wks("pdf", figpath+figname)
nplot = 3
plot = new(nplot, graphic)

res = True
res@gsnDraw  = False
res@gsnFrame = False

res@vpHeightF= 0.4                    ; change aspect ratio of plot
res@vpWidthF = 0.8     

res@trXMinF = 6
res@trXMaxF = 30           

res@gsnLeftString = "CTRL"
plot(0) = gsn_csm_xy(wks, xaxis, index1, res)
res@gsnLeftString = "NOQ1"
plot(1) = gsn_csm_xy(wks, xaxis, index2, res)
res@gsnLeftString = "NOTP_NOQ1"
plot(2) = gsn_csm_xy(wks, xaxis, index3, res)


pres = True
pres@gsnMaximize = True
pres@gsnPanelMainString = "Webster-Yang Index"

gsn_panel(wks, plot, (/3,1/), pres)

end