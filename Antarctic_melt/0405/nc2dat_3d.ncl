begin

yrStrt = 200412
yrLast = 200501

fin = addfile("/home/yangsong3/huxm/antarctic_melt/199012_201001_ta_4time_daily.nc", "r")
lat = fin->latitude
lon = fin->longitude
printVarSummary(lat)

filiname = (/"cc","cli","clw","t","q","o3"/)
varname = (/"cc","ciwc","clwc","t","q","o3"/)
nvar = dimsizes(varname)

outpath = "./"
outname = (/"cc_warm.dat","ciwc_warm.dat","clwc_warm.dat","ta_warm.dat","hus_warm.dat","o3_warm.dat"/)

diri = "/home/yangsong3/wuyt/sysu/scripts/Antarctic_melt/data/"

nlat = 181
nlon = 360
nlev = 37
x_med = new((/62,nlev,nlat,nlon/), float)
do i = nvar-1, nvar-1
    fili1 = diri+filiname(i)+".1982_1991_2004_2015_12.daily.nc"
    fili2 = diri+filiname(i)+".1983_1992_2005_2016_01.daily.nc"
    f1 = addfile(fili1, "r")
    f2 = addfile(fili2, "r")
    TIME1 := f1->time
    YYYY1 := cd_calendar(TIME1, -1)
    iYYYY1 := ind(YYYY1.eq.yrStrt)
    printVarSummary(iYYYY1)
    TIME2 := f2->time
    YYYY2 := cd_calendar(TIME2, -1)
    iYYYY2 := ind(YYYY2.eq.yrLast)
    printVarSummary(iYYYY2)
    x1 = short2flt(f1->$varname(i)$(iYYYY1,:,:,:))
    x2 = short2flt(f2->$varname(i)$(iYYYY2,:,:,:))
    ; printVarSummary(x1)
    ; printVarSummary(x2)
    x_med(0:30,:,:,:) = x1
    x_med(31:61,:,:,:) = x2
    copy_VarCoords(x1(0,:,:,:), x_med(0,:,:,:))

    ; x = linint2_Wrap(x_med&longitude, x_med&latitude, x_med, True, lon, lat, 0)
    x = f2fsh_Wrap(x_med, (/dimsizes(lat),dimsizes(lon)/))
    printVarSummary(x)

    system("rm -f "+outpath+outname(i))
    fbindirwrite(outpath+outname(i),x)
    ;=======================================================================
    ; x_zonal = dim_avg_n_Wrap(x, 3)
    ; wks = gsn_open_wks("pdf", varname(i))

    ; gsn_define_colormap(wks,"BlueDarkRed18")
    ; plot = new(20,graphic)

  ;   res                  = True                      ; plot mods desired
  ;   res@gsnFrame         = False
  ;   res@gsnDraw          = False
  ;   res@gsnSpreadColors  = True

  ;   res@cnFillOn         = True                      ; turn on color
  ;   res@cnLinesOn        = False
  ;   res@cnLineLabelsOn   = False                     ; no contour line labels
  ;   res@cnInfoLabelOn    = False                     ; no info label
  ;   res@lbLabelBarOn     = False

  ;   res@tmYRMode   = "Automatic"
  ;   res@tmYLMode   = "Explicit"
  ;   res@tmYLValues = (/"300","100","70","50","30","20","10","7","5"/)
  ;   res@tmYLLabels = (/"300","100","70","50","30","20","10","7","5"/)
  ;   res@tmYUseLeft = False

  ;   ; res@cnLevelSelectionMode = "ManualLevels"
  ;   ; res@cnMinLevelValF       = 4e-7
  ;   ; res@cnMaxLevelValF       = 4e-5
  ;   ; res@cnLevelSpacingF      = 8e-7

  ;   do j = 0,19
  ;      res@gsnLeftString     = "day "+(j*3+1)
  ;      plot(j) = gsn_csm_pres_hgt(wks, x_zonal(j*3,{300:5},:), res)
  ;   end do

  ; ;************************************************
  ; ; create panel
  ; ;************************************************
  ;   resP                     = True
  ;   resP@gsnFrame            = True
  ;   resP@gsnMaximize         = True
  ;   resP@gsnPanelLabelBar    = True                ; add common colorbar
  ;   resP@txString            = outname(i)
  ;   gsn_panel(wks,plot,(/4,5/),resP)

;===========================================
wks = gsn_open_wks("pdf", varname(i))

plot = new(1, graphic)

res = True
res@gsnDraw  = False
res@gsnFrame = False
res@gsnAddCyclic = False
res@gsnMaximize = True

; res@mpMinLatF = latstr
; res@mpMaxLatF = latend
; res@mpMinLonF = lonstr
; res@mpMaxLonF = lonend

res@mpOutlineOn = True
res@mpFillOn = True
res@mpGeophysicalLineThicknessF = 2
res@mpGeophysicalLineColor = "grey25" ;"red"
; res@mpFillDrawOrder = "PostDraw"
res@mpCenterLonF = 180

res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False

res@gsnStringFontHeightF = 0.025
res@tmXBLabelFontHeightF = 0.018
res@tmYLLabelFontHeightF = 0.018

plot(0) = gsn_csm_contour_map(wks, x(22,10,:,:), res)
draw(plot)
end do
end
