
begin

fq = addfile("./Q1.daily.6-30.nc", "r")


TIME = fq->time
YYYYMM = cd_calendar(TIME, -1)
MM = mod(YYYYMM, 100)

iMM = ind(MM.ge.3 .and. MM.le.5)
Q1 = fq->q1(:,iMM,:,:)
q1 = dim_avg_n_Wrap(Q1, 0)
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

Q1_area = runave_Wrap(q1_area, 9, 0)

qq = new((/2,92/), float)
qq(0,:) = q1_area
qq(1,0:dimsizes(Q1_area)-1) = Q1_area
;=========================================================
q1_specx = specx_anal(q1_area, 1, 3, 0.1)
printVarSummary(q1_specx)
freq = q1_specx@frq
spcx = q1_specx@spcx
;=======================================================
; plot 
;=======================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20191130/"
figname = "Q1_MAM_series_90day_CTRL"
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

res@tiYAxisString = "Q1 over the TP in Spring"
plot(0) = gsn_csm_xy(wks, ispan(1, 92, 1), qq, res)

; res@tiYAxisString = "freqency"
; plot(1) = gsn_csm_xy(wks, ispan(1, 46, 1), freq, res)

res@trYMaxF = 400
res@tiYAxisString = "Spectral analysis"
plot(1) = gsn_csm_y(wks, spcx, res)
; draw(plot(0))
; frame(wks)

pres = True
pres@gsnMaximize = True

gsn_panel(wks, plot, (/2,1/), pres)



end