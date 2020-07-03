
begin

casename = (/"CTRL","NOTP","NOQ1","NOTP_NOQ1"/) ;
ncase = dimsizes(casename)

varname = (/"SH","LH","netS","netR"/)
nvar = dimsizes(varname)

month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
nm = dimsizes(month)

flux = new((/ncase,nm,96,144/), float)

do i = 0, ncase-1
  diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(i)+"_F2000_F19G16_CAM4/pro/"
  filiu = "F2000_"+casename(i)+".SHFLX.mon.0601-3012.nc"
  filiv = "F2000_"+casename(i)+".LHFLX.mon.0601-3012.nc"
  filis = "F2000_"+casename(i)+".FSNS.mon.0601-3012.nc"
  filil = "F2000_"+casename(i)+".FLNS.mon.0601-3012.nc"

  fu = addfile(diri+filiu, "r")
  fv = addfile(diri+filiv, "r")
  fs = addfile(diri+filis, "r")
  fl = addfile(diri+filil, "r")

  SHFLX = fu->SHFLX
  LHFLX = fv->LHFLX
  FSNS  = fs->FSNS
  FLNS  = fl->FLNS

  FLUX  = SHFLX
  FLUX  = SHFLX+LHFLX+FSNS+FLNS

  flux(i,:,:,:) = clmMonTLL(FLUX)

  ; sh = clmMonTLL(SHFLX)
  ; lh = clmMonTLL(LHFLX)
  ; netS = clmMonTLL(FSNS)
  ; netR = clmMonTLL(FLNS)

  ; flux(i,0,:,:,:) = sh
  ; flux(i,1,:,:,:) = lh
  ; flux(i,2,:,:,:) = netS
  ; flux(i,3,:,:,:) = netR

end do
printVarSummary(flux)


;flux_series = flux(:,:,{30},{90})
flux_series = wgt_areaave_Wrap(flux(:,:,{28:40},{75:105}), 1, 1, 0)


x = new((/6,12/), float)
x(0:1,:) = flux_series(0:1,:)
x(2,:) = flux_series(0,:) - flux_series(1,:)
x(3:4,:) = flux_series(2:3,:)
x(5,:) = flux_series(2,:) - flux_series(3,:)


;===================================================================
; Plot
;===================================================================

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20200106/"
figname = "check_4_experiments_surface_heat_flux"
Month= ispan(1, 12, 1)

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
res@tmXBLabels = month
res@tmXTOn = False
;res@tiXAxisString = "Month"

; res@trYMinF = -2
; res@trYMaxF = 1
; res@tmYLMode = "Explicit"
; res@tmYLValues = ispan(-, 0, 1)
; res@tmYLLabels = ispan(-5, 0, 1)
res@tmYROn = False
res@tiYAxisString = "W/m~S~2~N~"
res@gsnYRefLine = 0
res@gsnYRefLineDashPattern = 16
res@gsnYRefLineColor = "gray20"
res@gsnYRefLineThicknessF = 5

Colors = (/"red","blue","black","goldenrod3","lightblue","forestgreen"/)  ;,"purple",
Patterns = (/0,0,0,16,16,16/) 
res@xyLineThicknessF = 8
res@xyLineColors = Colors
res@xyDashPatterns = Patterns

res@gsnLeftString = "surface heat flux"

plot(0) = gsn_csm_xy(wks, Month, x, res)


;==================================================================================
    res_lines                   = True                  ; polyline mods desired
    res_lines@gsLineThicknessF  = 5                     ; line thicker
    

    res_text                    = True                  ; text mods desired
    res_text@txFontHeightF      = 0.018                 ; change text size
    res_text@txJust             = "CenterLeft"          ; text justification

    Texts  = (/"Ctrl","NT","Ctrl-NT","NH","NTNH","NH-NTNH"/)

; ;============================================================
    interval = 50
    gx = (/1.5,2.2/)
    text_pos = 2.3
    gy1 = (/550,550/)
    gy2 = gy1-interval ;(/15,15/)
    gy3 = gy2-interval
    gy4 = gy1
    gy5 = gy2 ;gy4-interval
    gy6 = gy3 ;gy5-interval

    res_lines@gsLineColor = Colors(0)
    res_lines@gsLineDashPattern = Patterns(0)
    dum11 = gsn_add_polyline(wks, plot(0), gx, gy1, res_lines)
    dum12 = gsn_add_text(wks, plot(0), Texts(0), text_pos, gy1, res_text)
    res_lines@gsLineColor = Colors(1)
    res_lines@gsLineDashPattern = Patterns(1)
    dum21 = gsn_add_polyline(wks, plot(0), gx, gy2, res_lines)
    dum22 = gsn_add_text(wks, plot(0), Texts(1), text_pos, gy2, res_text)
    res_lines@gsLineColor = Colors(2)
    res_lines@gsLineDashPattern = Patterns(2)
    dum31 = gsn_add_polyline(wks, plot(0), gx, gy3, res_lines)
    dum32 = gsn_add_text(wks, plot(0), Texts(2), text_pos, gy3, res_text)
    

    text_pos = 9.8
    gx = (/9,9.7/)
    res_lines@gsLineColor = Colors(3)
    res_lines@gsLineDashPattern = Patterns(3)
    dum41 = gsn_add_polyline(wks, plot(0), gx, gy4, res_lines)
    dum42 = gsn_add_text(wks, plot(0), Texts(3), text_pos, gy4, res_text)
    res_lines@gsLineColor = Colors(4)
    res_lines@gsLineDashPattern = Patterns(4)
    dum51 = gsn_add_polyline(wks, plot(0), gx, gy5, res_lines)
    dum52 = gsn_add_text(wks, plot(0), Texts(4), text_pos, gy5, res_text)
    res_lines@gsLineColor = Colors(5)
    res_lines@gsLineDashPattern = Patterns(5)
    dum61 = gsn_add_polyline(wks, plot(0), gx, gy6, res_lines)
    dum62 = gsn_add_text(wks, plot(0), Texts(5), text_pos, gy6, res_text)


draw(plot(0))
;draw(plot(1))
frame(wks)




end