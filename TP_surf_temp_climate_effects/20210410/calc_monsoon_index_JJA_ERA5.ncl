
begin

yrStrt = 1979
yrLast = 2019
year = ispan(yrStrt, yrLast, 1)

lonstr = 40
lonend = 150
latstr = -15
latend = 50

diri = "/home/yangsong3/data-observation/ERA5-monthly/pressure/"
filiu = "u_mon_1x1_1979_2019.nc"
filiv = "v_mon_1x1_1979_2019.nc"


fu = addfile(diri+filiu, "r")
fv = addfile(diri+filiv, "r")
Time = fu->time
YYYY = cd_calendar(Time, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
u1 = short2flt(fu->u(iYYYY,:,{latstr:latend},{lonstr:lonend}))
v1 = short2flt(fv->v(iYYYY,:,{latstr:latend},{lonstr:lonend}))
printVarSummary(u1)

u = month_to_season(u1, "JJA")
v = month_to_season(v1, "JJA")
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
copy_VarCoords(u1_wnp, index_wnp)
printVarSummary(index_wnp)

;===============================================================================
; Dynamic Indian monsoon index Wang etal (2001)
; (U850 (5N-15N, 40E-80E) - (U850 20N-30N, 70E-90E))
;===============================================================================
u_sa = u(:,{850},:,:)
u1_sa = wgt_areaave_Wrap(u_sa(:,{5:15},{40:80}), 1, 1, 0)
u2_sa = wgt_areaave_Wrap(u_sa(:,{20:30},{70:90}), 1, 1, 0)
index_sa = u1_sa - u2_sa
copy_VarCoords(u1_sa, index_sa)
printVarSummary(index_sa)

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
; wtite data to txt file
;===============================================================================
filipath = "./"
filiname = "calc_monsoon_index_JJA_ERA5.txt"

write_table(filipath+filiname, "w", [/"Year","WY","WNP","SA","Gosiwami","EA"/], "%4s %10s %10s %10s %10s %10s")
write_table(filipath+filiname, "a", [/year,index_wy,index_wnp,index_sa,index_rm1,index_rm2/], "%4.2i %10.2f %10.2f %10.2f %10.2f %10.2f")

exit
;==================================================================
;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/"
figname = "calc_monsoon_index_JJA_WY_CTRL_S1_S2"
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