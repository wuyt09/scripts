
begin

lonstr = 75
lonend = 105
latstr = 25
latend = 40
yrStrt = 198102
yrLast = 200601

casename = (/"Ctrl_FAMIP","Hist_FAMIP"/)
ncase = dimsizes(casename)

season = "DJF"
SC = new((/ncase,96,144/), float)
AL = SC

do i = 0, ncase-1
    diri = "/home/yangsong3/data-model/wuyt/TPbias/TPbias_"+casename(i)+"/a20191206/"
    fili1 = casename(i)+".cam.h0.FSDS.1979-2005.nc"
    fili2 = casename(i)+".cam.h0.FSNS.1979-2005.nc"
    fili3 = casename(i)+".clm2.h0.FSNO.1979-2005.nc"

    f1 = addfile(diri+fili1, "r")
    f2 = addfile(diri+fili2, "r")
    f3 = addfile(diri+fili3, "r")

    TIME = f1->time
    YYYY = cd_calendar(TIME, -1)
    iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
    DS = f1->FSDS(iYYYY,:,:)
    NS = f2->FSNS(iYYYY,:,:)
    SCX = f3->FSNO(iYYYY,:,:)

    ; time = f1->time(iYYYY)
    ; lat = f1->lat
    ; lon = f1->lon
    ; oned = ndtooned(DS)
    ; inx = ind_resolve(ind(oned.eq.0), dimsizes(DS))
    ; print(time(inx(:,0)) +"  "+lat(inx(:,1))+"  "+lon(inx(:,2)))
    ; exit


    ALBEDO = NS
    DS = where(DS.eq.0, 0.00000001, DS)

    ALBEDO = (DS-NS)/DS
    AL0 = month_to_season(ALBEDO, season)
    SC0 = month_to_season(SCX, season)
    AL(i,:,:) = dim_avg_n_Wrap(AL0, 0)
    SC(i,:,:) = dim_avg_n_Wrap(SC0, 0)
delete([/diri,fili1,fili2,fili3,f1,f2,f3,SCX,DS,NS,ALBEDO,AL0,SC0/])
end do
printVarSummary(AL)

al = AL(0,:,:)
al = AL(1,:,:) - AL(0,:,:)
sc = SC(0,:,:)
sc = SC(1,:,:) - SC(0,:,:)

;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)

;=========================================================
; plot 
;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/a20200207/"
figname = season+"_snowcover_albedo"

; lonstr = 60
; lonend = 120
; latstr = 10
; latend = 50

wks = gsn_open_wks("pdf", figpath+figname)
plot = new(6, graphic)
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
res@cnFillPalette         = "MPL_RdBu"
;res@lbLabelBarOn          = False

res@mpMinLatF            = latstr      ; range to zoom in on
res@mpMaxLatF            = latend
res@mpMinLonF            = lonstr
res@mpMaxLonF            = lonend
; res@mpCenterLonF         = 180.
res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025

res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -0.1
res@cnMaxLevelValF = 0.1
res@cnLevelSpacingF = 0.01
res@gsnStringFontHeightF = 0.035
;res@gsnRightString       = "K"


res@gsnLeftString = "ALBEDO"
plot(0) = gsn_csm_contour_map(wks, al, res)
plot_topo(0) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(0), plot_topo(0))

res@gsnLeftString = "SnowCover"
plot(1) = gsn_csm_contour_map(wks, sc, res)
plot_topo(1) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(1), plot_topo(1))


res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = 0.1
res@cnMaxLevelValF = 0.8
res@cnLevelSpacingF = 0.1
res@cnFillPalette = "WhiteBlue"
res@gsnLeftString = "Ctrl: ALBEDO"
plot(2) = gsn_csm_contour_map(wks, AL(0,:,:), res)
plot_topo(2) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(2), plot_topo(2))

res@gsnLeftString = "Ctrl: SnowCover"
plot(3) = gsn_csm_contour_map(wks, SC(0,:,:), res)
plot_topo(3) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(3), plot_topo(3))

res@gsnLeftString = "Hist: ALBEDO"
plot(4) = gsn_csm_contour_map(wks, AL(1,:,:), res)
plot_topo(4) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(4), plot_topo(4))

res@gsnLeftString = "Hist: SnowCover"
plot(5) = gsn_csm_contour_map(wks, SC(1,:,:), res)
plot_topo(5) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(5), plot_topo(5))

pres = True
pres@gsnMaximize = True
;pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = season+" mean"
; pres@pmLabelBarOrthogonalPosF = -0.02
; pres@lbLabelFontHeightF = 0.015

gsn_panel(wks, plot, (/3,2/), pres)

end