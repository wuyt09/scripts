
begin

lonstr = 40
lonend = 150
latstr = -15
latend = 50

yrStrt = 1979
yrLast = 2015
year = ispan(yrStrt, yrLast, 1)
nyear = dimsizes(year)

x = asciiread("../20200704/TTT_TP_dtrend_standardize_Mar_Apr_May_ERAI.txt",(/nyear,3/),"float")
; x = indxx(:,0)
; print(x)
; ind1 = ind(x.ge.1)
; ind2 = ind(x.le.-1)
; print(ind1)
; print(ind2)
; exit
; num1 = dimsizes(ind1)
; num2 = dimsizes(ind2)

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
Month = (/"May","June"/)
monS = 6
filipath = "./"
filiname = "Comp_"+Month(monS-5)+"_monsoon_index_with_TTT_ERAI.txt"
inx_hig = new((/3,4/), float)
inx_low = inx_hig
inx_may = indx(:,monS-1::12)
printVarSummary(inx_may)

do i = 0, 2
    ind1 := ind(x(:,i).ge.1.)
    ind2 := ind(x(:,i).le.-1.)
    inx_hig(i,:) = dim_avg_n_Wrap(inx_may(:,ind1), 1)
    inx_low(i,:) = dim_avg_n_Wrap(inx_may(:,ind2), 1)
end do
printVarSummary(inx_low)

month_hig = (/"March_hig","April_hig","May_hig"/)
month_low = (/"March_low","April_low","May_low"/)
write_table(filipath+filiname, "w", [/"Comp","WY","WNP","SA","EA"/], "%10s %10s %10s %10s %10s")
write_table(filipath+filiname, "a", [/month_hig,inx_hig(:,0),inx_hig(:,1),inx_hig(:,2),inx_hig(:,3)/], "%10s %10.2f %10.2f %10.2f %10.2f")
write_table(filipath+filiname, "a", [/month_low,inx_low(:,0),inx_low(:,1),inx_low(:,2),inx_low(:,3)/], "%10s %10.2f %10.2f %10.2f %10.2f")

; monS = 6
; filipath = "./"
; filiname = "Corr_monsoon_index_with_TPQ1_ERAI_June.txt"
; corr = new((/3,4/), float)
; do i = 0, 2
;     do j = 0, 3
;         qq = ttt(i+2::12)
;         pp = indx(j,monS-1::12)
;         corr(i,j) = escorc(qq, pp)
;     end do
; end do
; write_table(filipath+filiname, "w", [/"Month","WY","WNP","SA","EA"/], "%10s %10s %10s %10s %10s")
; write_table(filipath+filiname, "a", [/Month,corr(:,0),corr(:,1),corr(:,2),corr(:,3)/], "%10s %10.2f %10.2f %10.2f %10.2f")

; exit


end