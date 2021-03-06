
begin
yrStrt = 1979
yrLast = 2016
lonstr = 80
lonend = 100

casename = (/"CTRL","NOTP"/) ;,"NOQ1","NOTP_NOQ1"
ncase = dimsizes(casename)
month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
nm = dimsizes(month)

LEV = (/200/)
u = new((/ncase,12,96,144/), float)


do i = 0, dimsizes(casename)-1
    diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(i)+"_F2000_F19G16_CAM4/pro/"
    filiu = "F2000_"+casename(i)+".plevel.U.mon.0601-3012.nc"
    filiv = "F2000_"+casename(i)+".plevel.V.mon.0601-3012.nc"


    fu = addfile(diri+filiu, "r")
    fv = addfile(diri+filiv, "r")
    u1 = fu->U(:,{LEV},:,:)
    ; v1 = fv->V(:,{LEV},:,:)
    ; v_zonal = dim_avg_n_Wrap(v1, 1)
    ; v1 = v1 - conform_dims(dimsizes(v1), v_zonal, (/0,2/))
    ; printVarSummary(u1)
    
    u(i,:,:,:) = clmMonTLL(u1)
end do

P_band = dim_avg_n_Wrap(u(:,:,:,{lonstr:lonend}), 3)
P_band!0 = "case"
P_band&case = casename
p_band = P_band(case|:,lat|:,month|:)
printVarSummary(P_band)


; ;==================================================================
; dirip = "/home/yangsong3/data-observation/GPCP/"
; filip = "precip.mon.mean.nc"

; fp = addfile(dirip+filip, "r")
; TIME = fp->time
; YYYY = cd_calendar(TIME, -1)/100
; iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
; precip = fp->precip(iYYYY,:,:) ;{latstr:latend},{lonstr:lonend}
; printVarSummary(precip)

; gpcp = clmMonTLL(precip)
; gpcp_band = dim_avg_n_Wrap(gpcp(:,:,{lonstr:lonend}), 2)
; GPCP_band = gpcp_band(lat|:,month|:)

;==================================================================
;==================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20200106/"


  figname = "Ctrl_NT_Jetstream_Evolution"

  wks = gsn_open_wks("pdf", figpath+figname)
  nplot = ncase
  plot = new(nplot, graphic)

  res = True
  res@gsnDraw  = False
  res@gsnFrame = False
  res@vpHeightF = 0.4
  res@vpWidthF = 0.7

  res@cnFillOn = True
  res@cnLinesOn = False
  res@cnLineLabelsOn = False
  res@cnInfoLabelOn = False
  res@lbLabelBarOn = False

  res@gsnStringFontHeightF = 0.035
  res@tmXBLabelFontHeightF = 0.025
  res@tmYLLabelFontHeightF = 0.025
  res@pmLabelBarOrthogonalPosF = 0.2
  ;------------------------------------------
      
  res@gsnRightString = "m/s"

  res@cnLevelSelectionMode = "ManualLevels"
  res@cnLevelSpacingF = 5
  res@cnMaxLevelValF = 40.
  res@cnMinLevelValF = 10.
  ; res@cnLevelSelectionMode = "ExplicitLevels"
  ; res@cnLevels = (/0.5,1,2,4,6,8,10,13/)
  ; res@cnFillColors = (/2,3,4,6,7,9,10,11,12/)
  res@cnFillPalette = "WhiteBlue"
  

  ; res@tmYLMode = "Manual"
  ; res@tmYLTickStartF = 10
  ; res@tmYLTickEndF = 50
  ; res@tmYLTickSpacingF = 10
  res@trYMaxF = 50
  res@trYMinF = 15

  res@tmXBMode = "Explicit"
  res@tmXBValues = ispan(0, 11, 1)
  res@tmXBLabels = month ;(2:9)

  ; res@gsnLeftString = "GPCP"
  ; plot(0) = gsn_csm_lat_time(wks, GPCP_band, res)

  ; res@gsnLeftString = "Ctrl"
  ; plot(1) = gsn_csm_lat_time(wks, p_band(0,:,:), res)

  ; res@gsnLeftString = "NH"
  ; plot(2) = gsn_csm_lat_time(wks, p_band(1,:,:), res)

  ; res@gsnLeftString = "NTNH"
  ; plot(3) = gsn_csm_lat_time(wks, p_band(2,:,:), res)
  
  leftsti = (/"Ctrl","NT"/)
  do i = 0, ncase -1
    res@gsnLeftString = leftsti(i)
      plot(i) = gsn_csm_lat_time(wks, p_band(i,:,:), res)
  end do

  pres = True
  pres@gsnMaximize = True
  pres@gsnPanelMainString = lonstr+"E - "+lonend+"E"
  pres@gsnPanelLabelBar = True
  ;pres@pmLabelBarOrthogonalPosF = -0.01

  gsn_panel(wks, plot, (/2,1/), pres)



end