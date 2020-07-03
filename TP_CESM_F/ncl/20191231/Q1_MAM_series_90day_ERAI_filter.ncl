
begin

yrStrt = 1979
yrLast = 2016
year = ispan(yrStrt, yrLast, 1)
nyr = dimsizes(year)

;==================================================================

tpmaskpath = "/home/yangsong3/wuyt/data/"
tpmaskname = "topo.sfc.1x1.nc"
ftp = addfile(tpmaskpath+tpmaskname, "r")

tp = ftp->topo(::-1,:)
printVarSummary(tp)

;==================================================================

diri = "/home/yangsong3/wuyt/data/ERA-Interim/daily/Q1/"
Q1 = new((/nyr,92/), float)

do i = 0, nyr-1
    yr = sprinti("%4i", year(i))
    printVarSummary(yr)
    fili = "q1."+yr+".nc"
    f = addfile(diri+fili, "r")
    TIME = f->time
    YYYYMM = cd_calendar(TIME, -1)
    MM = mod(YYYYMM, 100)
    iMM = ind(MM.ge.3 .and. MM.le.5)
    Q1X = f->q1(iMM,:,::-1,:)
    tpmask = linint2_Wrap(tp&lon, tp&lat, tp, True, Q1X&lon, Q1X&lat, 0)
    Q1X = mask(Q1X, tpmask.ge.1500, True)
    ;printVarSummary(Q1)
    Q1X = Q1X*1004.07/9.8

    plev = Q1X&lev
    plev = plev*100
    psfc = max(plev)
    ptop = min(plev)
    dP = dpres_plevel_Wrap(plev, psfc, ptop, 0)
    dp = conform_dims(dimsizes(Q1X), dP, 1)
    Q1dp = Q1X*dp
    q1dp = dim_sum_n_Wrap(Q1dp, 1)/sum(dP)
    copy_VarCoords(Q1X(:,0,:,:), q1dp)
    Q1_area = wgt_areaave_Wrap(q1dp(:,{25:40},{75:105}), 1, 1, 0)
    Q1(i,:) = Q1_area
delete([/f,TIME,YYYYMM,MM,iMM,Q1X/])
end do
printVarSummary(Q1)

;=========================================================
; Q1_area = wgt_areaave_Wrap(Q1(:,:,{25:40},{75:105}), 1, 1, 0)
; time = Q1_area&time
; printVarSummary(Q1_area)
; ;Q1_area = dtrend_msg_n(time, Q1_area, False, False, 1)
; ;dtrend(Q1_area, False)
; ;printVarSummary(Q1_area)

q1 = ndtooned(Q1)

;=========================================================
; q1_specx = Q1_area
; spcx = new((/25,46/), float)
; do i = 0, 24
;     q1_specx(i,:) = specx_anal(Q1_area(i,:), 1, 3, 0.1)
;     spcx(i,:) = q1_specx@spcx
; end do

iopt  = 0
jave  = 1
pct   = 0.1
qdof = specx_anal(q1, iopt, jave, pct)
splt = specx_ci(qdof,0.05,0.95)
printVarSummary(qdof)
printVarSummary(splt)
spcx = qdof@spcx

f = qdof@frq
printVarSummary(f)
p = 1/f
p!0 = "f"
p&f = f
p@long_name = "period"
p@units     = "day"
ip = ind(p.le.80) 

;=========================================================

ca = 8
cb = 2.5
; ca = 30
; cb = 10
fca = 1./ca
fcb = 1./cb
q1_bf = bw_bandpass_filter(q1, fca, fcb, True, 0)
copy_VarCoords(q1, q1_bf)
q1_bf@long_name = "Band Pass: "+cb+"-"+ca+" day"
q1_bf = dim_standardize_n_Wrap(q1_bf, 1, 0)
printVarSummary(q1_bf)

x = new(12, float)
index = new(25*92, logical)
do i = 3, 25*92-4
    x(0) = q1_bf(i) - q1_bf(i-3)
    x(1) = q1_bf(i) - q1_bf(i-2)
    x(2) = q1_bf(i) - q1_bf(i-1)
    x(3) = q1_bf(i) - q1_bf(i+3)
    x(4) = q1_bf(i) - q1_bf(i+2)
    x(5) = q1_bf(i) - q1_bf(i+1)

    ; x(6) = q1_bf(i) - q1_bf(i-4)
    ; x(7) = q1_bf(i) - q1_bf(i-5)
    ; x(8) = q1_bf(i) - q1_bf(i-6)
    ; x(9) = q1_bf(i) - q1_bf(i+4)
    ; x(10) = q1_bf(i) - q1_bf(i+5)
    ; x(11) = q1_bf(i) - q1_bf(i+6)
    index(i) = .not.any(x.lt.0) .and. q1_bf(i).ge.1 ;
    ;print(i+" "+index(i))
end do
ind1 = ind(index)
printVarSummary(ind1)
q_highpoint = q1_bf(ind1)
printVarSummary(q_highpoint)
q_higavg = avg(q_highpoint)
print(q_higavg)

ind2 = ind(q1_bf.ge.q_higavg)
printVarSummary(ind2)

; print(ind1)
; print("=======================")
; print(ind2)


; filename = "./high_point_index_"+cb+"-"+ca+"day.txt"
; write_table(filename, "w", [/ind1/], "%10i")
; exit

Q1_bf = runave_Wrap(q1_bf, 5, 0)

qq = new((/2,92*38/), float)
qq(0,:) = q1_bf
qq(1,0:dimsizes(Q1_bf)-1) = Q1_bf



;=======================================================
; plot 
;=======================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20191231/"
figname = "Q1_MAM_series_90day_ERAI_filter"
wks = gsn_open_wks("pdf", figpath+figname)
plot = new(1, graphic)

; res = True
; res@gsnDraw  = False
; res@gsnFrame = False
; res@gsnMaximize = True

; res@vpWidthF = 0.7
; res@vpHeightF = 0.4

; res@trXMinF = 999
; res@trXMaxF = 1100
; ; res@tmXBMode = "Explicit"
; ; res@tmXBValues = year
; ; res@tmXBLabels = Month
; ; res@tmXTOn = False
; ; res@tiXAxisString = "Month"

; ;res@xyLineThicknessF = 5
; res@xyLineColors = (/"black", "blue", "red", "green4"/)
; ;(/"CTRL", "HIST", "FUTR", "FUBIAS"/)
; res@xyDashPattern = 0
; res@gsnYRefLine = q_higavg
; res@gsnYRefLineDashPattern = 2
; res@gsnYRefLineThicknesses = 4
; res@gsnYRefLineColor = "blue"

; res@tiYAxisString = "Q1 over the TP in Spring"
; plot(0) = gsn_csm_xy(wks, ispan(1, 92*38, 1), qq, res)
; draw(plot(0))


res2 = True
res2@gsnDraw  = False
res2@gsnFrame = False
res2@gsnMaximize = True
res2@vpWidthF = 0.7
res2@vpHeightF = 0.4

; res2@trXMinF = 0
; res2@trXMaxF = 50

res2@xyLineColors = (/"black", "blue", "red", "green4"/)
res2@xyDashPatterns = (/0,2,2,2/)
res2@xyLineThicknessF = 5
; res@gsnXRefLine = (/6,24/)
; res@gsnXRefLineDashPattern = 2
; res@gsnXRefLineThicknesses = 4
; res@gsnXRefLineColor = "blue"
; do i = 0, 24
;     plot(i) = gsn_csm_y(wks, spcx(i,:), res)
; end do

; pres = True
; pres@gsnMaximize = True
; gsn_panel(wks, plot, (/6,4/), pres)

; plot(1) = gsn_csm_y(wks, spcx, res2)
; draw(plot(1))

;plot(0) = gsn_csm_xy(wks, p(ip), spcx(ip), res2)

plot(0) = gsn_csm_xy(wks, p(ip), splt(:,ip), res2)
draw(plot(0))
;gsn_panel(wks, plot, (/2,1/), True)
end
