
begin

diri = "/home/yangsong3/data-observation/GPCP/"
fili = "precip.pentad.clim.mean.nc"

f = addfile(diri+fili, "r")
pentad_prec = f->precip

prec_tlat = dim_avg_n_Wrap(pentad_prec(:,{10:20},:), 1)
; prec_tlat&time = ispan(1, 73, 1)
printVarSummary(prec_tlat)


;-------------------------------------------------------------------
Month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20200513/"
figname = "precip_evolution_10-20N_GPCP"

yaxis = ispan(1, 73, 2)

wks = gsn_open_wks("pdf", figpath+figname)

nplot = 1
plot = new(nplot, graphic)


res          = True                ; plot mods desired
res@gsnDraw  = False
res@gsnFrame = False

res@cnFillOn         = True                ; color on
res@cnFillPalette    = "CBR_wet"    ; set color map
;res@lbLabelBarOn = False
res@cnInfoLabelOn = False

res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = 5
res@cnMaxLevelValF = 20
res@cnLevelSpacingF  = 2                   ; contour spacing

; res@tmYLMode = "Manual"
; ; res@tmYLLabels = yaxis
; ; res@tmYLValues = yaxis
; res@tmYLTickEndF = 73
; res@tmYLTickStartF = 1
; res@tmYLTickSpacingF = 6
; ; res@tmYLMinorValues = 6
; ; res@tmYMajorGrid = True
res@trYMinF = 1
res@trYMaxF = 73
res@tmYLMode = "Explicit"
res@tmYLLabels = yaxis
res@tmYLValues = yaxis
res@tmYLLabelStride = 3
res@tmYLMinorOn = True

res@tmXBMode = "Explicit"
res@tmXBValues = (/40,60,80,100,120,140,160,180/)
res@tmXBLabels = (/"40E","60E","80E","100E","120E","140E","160E","180E"/)
res@tmXBMinorOn = False

; res@tmXBLabelStride  = 2                   ; tick mark label stride
res@tiYAxisString    = "Pentad"              ; y axis title
res@tmYLPrecision = 2

; res@lbLabelStride    = 2                   ; every other label
; res@cnLinesOn        = False               ; turn off contour lines
res@gsnStringFontHeightF = 0.03
res@gsnLeftString = "GPCP"
plot(0) = gsn_csm_time_lat(wks, prec_tlat(:,{40:180}), res)

draw(plot)
end