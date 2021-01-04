
begin

LEV = 200
nlat = 96
nlon = 144

lonstr = 106
lonend = 122
latstr = 20
latend = 35

casename = (/"CTRL","TPQ1_double"/) ;"NOTP"
ncase = dimsizes(casename)
year = ispan(6, 30, 1)
nyear = dimsizes(year)
month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
nmonth = dimsizes(month)
season = "MAM"

t = new((/ncase,nmonth,96,144/), float)
u = t
v = t

do i = 0, ncase-1
  diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(i)+"_F2000_F19G16_CAM4/pro/"
  zname = "F2000_"+casename(i)+".PS.mon.0601-3012.nc"
  tname = "F2000_"+casename(i)+".T.mon.0601-3012.nc"
  fz = addfile(diri+zname, "r")
  ft = addfile(diri+tname, "r")
  PS = fz->PS
  T = ft->T

  plev = T&lev
  plev = plev*100
  plev@units = "Pa"
  plev!0 = "lev"
  plev&lev = T&lev
  ; psfc = max(plev)
  ptop = min(plev)
  dp = dpres_plevel_Wrap(plev, PS, ptop, 0)

  Tdp = T
  Tdp = T*dp
  tdp = dim_sum_n_Wrap(Tdp(:,{500:200},:,:), 1)/dim_sum_n_Wrap(dp(:,{500:200},:,:),1)
  copy_VarCoords(T(:,0,:,:), tdp)
  t(i,:,:,:) = clmMonTLL(tdp)
  delete([/PS,T,dp,Tdp,tdp/])
end do
t!0 = "case"
t&case = casename
printVarSummary(t)

; ;------------ttest--------------------------------------
; z_var = dim_variance_n_Wrap(z, 1)
; z_ave = dim_avg_n_Wrap(z, 1)

; probz = ttest(z_ave(0,:,:), z_var(0,:,:), nyear, z_ave(1,:,:), z_var(1,:,:), nyear, False, False)

; diff_z = z_ave(0,:,:)
; diff_z = z_ave(0,:,:) - z_ave(1,:,:); diff_z = mask(diff_z, probz.lt.sig, True)
; printVarSummary(diff_z)
; ;----------------------------------------------------------------
tpmaskpath = "/home/yangsong3/wuyt/sysu/scripts/TP_CESM_F/model_need_input_data/"
tpmaskname = "TPmask_1500m_C20200713.nc.nc"
ftp = addfile(tpmaskpath+tpmaskname, "r")

tpmask = ftp->tp
printVarSummary(tpmask)


t_mask = mask(t, tpmask.ge.0.5, True)
t_TP = wgt_areaave_Wrap(t_mask, 1., 1., 0)
printVarSummary(t_TP)

;===================================================================
; Plot
;===================================================================
filipath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20210102/"
filiname = "model_12month_TP_TropoT"
wks = gsn_open_wks("pdf", filipath+filiname)
plot = new(1, graphic)
plot_topo = plot

res = True
res@gsnDraw  = False
res@gsnFrame = False
res@gsnMaximize = True

res@vpWidthF = 0.7
res@vpHeightF = 0.4

res@trXMinF = 1
res@trXMaxF = 12
res@tmXBMode = "Explicit"
res@tmXBValues = ispan(1, 12, 1)
res@tmXBLabels = month
res@tmXTOn = False
; res@tiXAxisString = "Month"

; res@trYMinF = -20
; res@trYMaxF = 20
; res@tmYLMode = "Explicit"
; res@tmYLValues = ispan(-20, 20, 5)
; res@tmYLLabels = ispan(-20, 20, 5)
; res@tmYROn = False
; res@tiYAxisString = "Temperature (~S~o~N~C)"
;res@gsnLeftString = "Tas: grid2station"

;Colors = (/"deeppink2", "blue", "purple","green4","orange2", "black" /);
Colors = (/"black","red","black","red"/)
Patterns = (/0,0,16,16/)
res@xyLineThicknessF = 15
res@xyLineColors = Colors
res@xyDashPatterns = Patterns

plot(0) = gsn_csm_xy(wks, ispan(1, 12, 1), t_TP, res)

;==================================================================================
    res_lines                   = True                  ; polyline mods desired
    res_lines@gsLineThicknessF  = 5                     ; line thicker
    

    res_text                    = True                  ; text mods desired
    res_text@txFontHeightF      = 0.018                 ; change text size
    res_text@txJust             = "CenterLeft"          ; text justification

    Texts  = casename ;(/"DJF TP Ts bias","MAM EA prec bias"/)

;============================================================
    gx = (/1.5,2.2/)
    text_pos = 2.3
    gy1 = (/292,292/)
    interval = 3
    gy2 = gy1 - interval

    res_lines@gsLineThicknessF = 15
    res_lines@gsLineColor = Colors(0)
    res_lines@gsLineDashPattern = Patterns(0)
    dum11 = gsn_add_polyline(wks, plot(0), gx, gy1, res_lines)
    dum12 = gsn_add_text(wks, plot(0), Texts(0), text_pos, gy1, res_text)
    res_lines@gsLineColor = Colors(1)
    res_lines@gsLineDashPattern = Patterns(1)
    dum21 = gsn_add_polyline(wks, plot(0), gx, gy2, res_lines)
    dum22 = gsn_add_text(wks, plot(0), Texts(1), text_pos, gy2, res_text)


draw(plot(0))
end