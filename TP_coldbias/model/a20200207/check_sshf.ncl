
begin
yrStrt = 1979
yrLast = 2016
year = ispan(yrStrt, yrLast, 1)
season = "JJA"

;fili = "/home/yangsong3/data-observation/ERA-interim/monthly/sshf.1979-2016.2time_4step.nc"
;fili = "/home/yangsong3/data-observation/wuyt/ERAI/monthly/2.5x2.5/shfx.lhfx.mon0.mean.7917.nc"
;fili = "/home/yangsong3/data-observation/wuyt/ERAI/monthly/1.5x1.5/sshf_and_slhf_era_mn_197901_201507_2time_4step.nc"
;fili = "/home/yangsong3/data-observation/wuyt/ERAI/monthly/1x1/sshf.mon.mean.7916.r1.nc"
;fili = "/home/yangsong3/data-model/wuyt/TPbias/TPbias_Ctrl_FAMIP/a20191206/Ctrl_FAMIP.cam.h0.SHFLX.1979-2005.nc"
fili = "/home/yangsong3/data-model/wuyt/TPbias/TPbias_Ctrl_FAMIP/a20191206/Ctrl_FAMIP.cam.h0.LHFLX.1979-2005.nc"
f = addfile(fili, "r")

;x = short2flt(f->sshf(0:455,:,:)) ;(0:3455,:,:)
x = f->LHFLX
printVarSummary(x)

; ; var1 = x(0::8,:,:)
; ; var2 = x(1::8,:,:)
; ; var3 = x(2::8,:,:)
; var4 = x(3::8,:,:)
; ; var5 = x(4::8,:,:)
; ; var6 = x(5::8,:,:)
; ; var7 = x(6::8,:,:)
; var8 = x(7::8,:,:)
; x1 = var4/(12*3600)+var8/(12*3600)
; x1 = x1/2
; x1 = -1*x1
; copy_VarCoords(x(::8,:,:), x1)
; printVarSummary(x1)
x1 = x ;-1*x
copy_VarCoords(x, x1)

; x_ann = month_to_annual(x1, 1)
; x_ann&year = year
; printVarSummary(x_ann)

; x_base = dim_avg_n_Wrap(x_ann({year|yrStrt:1997 },latitude|:,longitude|:), 0)
; x_warm = dim_avg_n_Wrap(x_ann({year|1998:yrLast},latitude|:,longitude|:), 0)

; diff = x_base
; diff = x_warm - x_base
; printVarSummary(diff)


; x_ann = month_to_season(x1, season)
; x_ann&time = year
; x_base = dim_avg_n_Wrap(x_ann({time|yrStrt:1997 },latitude|:,longitude|:), 0)
; x_warm = dim_avg_n_Wrap(x_ann({time|1998:yrLast},latitude|:,longitude|:), 0)
; diff = x_base
; diff = x_warm - x_base
; printVarSummary(diff)

x_ann = month_to_season(x1, season)
x_base = dim_avg_n_Wrap(x_ann, 0)
diff = x_base


;=========================================================
; plot 
;=========================================================
figpath = "./"
figname = "check_sshf"

lonstr = 75
lonend = 105
latstr = 25
latend = 40

wks = gsn_open_wks("pdf", figpath+figname)
plot = new(1, graphic)
plot_topo = plot

 res_topo            = True
 res_topo@cnFillOn   = False
 res_topo@cnLinesOn  = True
 res_topo@gsnFrame   = False
 res_topo@gsnDraw    = False
 res_topo@cnLineLabelsOn         = False
 res_topo@cnLineThicknessF       = 3
 res_topo@cnLineColor            = "grey"
 res_topo@cnInfoLabelOn          = False
 res_topo@cnLevelSelectionMode   ="ExplicitLevels"
 res_topo@cnLevels               = (/1500/)
 res_topo@gsnLeftString = ""
 res_topo@gsnRightString = ""


res                       = True     ; plot mods desired
res@gsnDraw               = False
res@gsnFrame              = False

res@gsnMaximize           = True     ; Make plot fit the frame
res@cnFillOn              = True     ; turn on color fill
res@cnLinesOn             = False    ; turn of contour lines
res@cnLineLabelsOn        = False
res@cnInfoLabelOn         = False
res@cnFillPalette         = "CBR_coldhot" ;MPL_RdBu
;res@lbLabelBarOn          = False

res@mpMinLatF            = latstr      ; range to zoom in on
res@mpMaxLatF            = latend
res@mpMinLonF            = lonstr
res@mpMaxLonF            = lonend
; res@mpCenterLonF         = 180.
res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025

; res@cnLevelSelectionMode = "ManualLevels"
; res@cnMinLevelValF = -6
; res@cnMaxLevelValF = 6
; res@cnLevelSpacingF = 1
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -100
res@cnMaxLevelValF = 100
res@cnLevelSpacingF = 10
res@gsnStringFontHeightF = 0.035


plot(0) = gsn_csm_contour_map(wks, diff, res)

draw(plot)

end