
begin

fq = addfile("./Q1.daily.6-30.nc", "r")


TIME = fq->time
YYYYMM = cd_calendar(TIME, -1)
MM = mod(YYYYMM, 100)

iMM = ind(MM.ge.3 .and. MM.le.5)
Q1 = fq->q1(:,iMM,:,:)
q1 = dim_avg_n_Wrap(Q1, 0) ;Q1(20,:,:,:) ;
printVarSummary(q1)

;==================================================================

tpmaskpath = "/home/yangsong3/wuyt/sysu/scripts/TP_CESM_F/model_need_input_data/"
tpmaskname = "TPmask_1000m_C20190410.nc"
ftp = addfile(tpmaskpath+tpmaskname, "r")

tpmask = ftp->tp
printVarSummary(tpmask)

;=========================================================
q1 = mask(q1, tpmask.ge.0.9, True)
q1_area = wgt_areaave_Wrap(q1, 1, 1, 0)
printVarSummary(q1_area)
;q1_area = dtrend(q1_area, False)
;=========================================================

q1_spex = specx_anal(q1_area, 0, 3, 0.1)
splt = specx_ci(q1_spex,0.05,0.95)
freq = q1_spex@frq
period = freq
period = 1/freq
printVarSummary(splt)
ip = ind(period.le.90)

;=========================================================
ca = 8
cb = 2.5
fca = 1./ca
fcb = 1./cb
q1_bf = bw_bandpass_filter(q1_area, fca, fcb, True, 0)
copy_VarCoords(q1_area, q1_bf)
q1_bf@long_name = "Band Pass: "+cb+"-"+ca+" day"
printVarSummary(q1_bf)

x = new(6, float)
index = new(92, logical)
do i = 3, 88
    x(0) = q1_bf(i) - q1_bf(i-3)
    x(1) = q1_bf(i) - q1_bf(i-2)
    x(2) = q1_bf(i) - q1_bf(i-1)
    x(3) = q1_bf(i) - q1_bf(i+3)
    x(4) = q1_bf(i) - q1_bf(i+2)
    x(5) = q1_bf(i) - q1_bf(i+1)
    index(i) = .not.any(x.lt.0) ;
    ;print(i+" "+index(i))
end do

; indx = ind(index)

; filename = "./high_point_index.txt"
; write_table(filename, "w", [/indx/], "%10i")


Q1_bf = runave_Wrap(q1_bf, 5, 0)

qq = new((/2,92/), float)
qq(0,:) = q1_bf
qq(1,0:dimsizes(Q1_bf)-1) = Q1_bf

;=======================================================
; plot 
;=======================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20191130/"
figname = "Q1_MAM_series_90day_CTRL_filter"
wks = gsn_open_wks("pdf", figpath+figname)
plot = new(2, graphic)

res = True
res@gsnDraw  = False
res@gsnFrame = False
res@gsnMaximize = True

res@vpWidthF = 0.7
res@vpHeightF = 0.4

; res@trXMinF = yrStrt
; res@trXMaxF = yrLast
; res@tmXBMode = "Explicit"
; res@tmXBValues = year
; res@tmXBLabels = Month
; res@tmXTOn = False
; res@tiXAxisString = "Month"

res@xyLineThicknessF = 5
res@xyLineColors = (/"black", "blue", "red", "green4"/)
;(/"CTRL", "HIST", "FUTR", "FUBIAS"/)
res@xyDashPattern = 0
res@gsnYRefLine = 0
res@gsnYRefLineDashPattern = 2
res@gsnYRefLineThicknesses = 4
res@gsnYRefLineColor = "blue"

;res@tiYAxisString = "Q1 over the TP in Spring"
;plot(0) = gsn_csm_xy(wks, ispan(1, 92, 1), qq, res)

plot(0) = gsn_csm_xy(wks, period(ip), q1_spex@spcx(ip), res)
draw(plot(0))



end
