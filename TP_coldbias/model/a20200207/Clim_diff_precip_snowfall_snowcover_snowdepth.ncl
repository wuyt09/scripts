
begin

lonstr = 75
lonend = 105
latstr = 25
latend = 40
yrStrt = 198102
yrLast = 200601

casename = (/"Ctrl_FAMIP","Hist_FAMIP"/)
ncase = dimsizes(casename)

pf = new((/ncase,12,96,144/), float)
;pf = new((/12,96,144/), float)
sf = pf
sc = pf
sd = pf

do i = 0, ncase-1
    diri = "/home/yangsong3/data-model/wuyt/TPbias/TPbias_"+casename(i)+"/a20191206/"
    fili1 = casename(i)+".cam.h0.PRECC.1979-2005.nc"
    fili2 = casename(i)+".cam.h0.PRECL.1979-2005.nc"
    fili3 = casename(i)+".clm2.h0.FSNO.1979-2005.nc"
    fili4 = casename(i)+".clm2.h0.H2OSNO.1979-2005.nc"
    fili5 = casename(i)+".cam.h0.PRECSC.1979-2005.nc"
    fili6 = casename(i)+".cam.h0.PRECSL.1979-2005.nc"

    f1 = addfile(diri+fili1, "r")
    f2 = addfile(diri+fili2, "r")
    f3 = addfile(diri+fili3, "r")
    f4 = addfile(diri+fili4, "r")
    f5 = addfile(diri+fili5, "r")
    f6 = addfile(diri+fili6, "r")

    TIME = f1->time
    YYYY = cd_calendar(TIME, -1)
    iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

    PF1 = f1->PRECC(iYYYY,:,:)
    PF2 = f2->PRECL(iYYYY,:,:)
    PF = PF1
    PF = (PF1+PF2)*1000*86400.
    PF@units = "mm/day"
    pf(i,:,:,:) = clmMonTLL(PF)

    SC = f3->FSNO(iYYYY,:,:)
    sc(i,:,:,:) = clmMonTLL(SC)

    SD = f4->H2OSNO(iYYYY,:,:)
    sd(i,:,:,:) = clmMonTLL(SD)

    SF1 = f5->PRECSC(iYYYY,:,:)
    SF2 = f6->PRECSL(iYYYY,:,:)
    SF = SF1 
    SF = (SF1+SF2)*1000*86400.
    SF@units = "mm/day"
    sf(i,:,:,:) = clmMonTLL(SF)
end do
printVarSummary(sf)


pf0 = pf(0,:,:,:)
pf0 = pf(1,:,:,:) - pf(0,:,:,:)

sc0 = sc(0,:,:,:)
sc0 = sc(1,:,:,:) - sc(0,:,:,:)

sd0 = sd(0,:,:,:)
sd0 = sd(1,:,:,:) - sd(0,:,:,:)

sf0 = sf(0,:,:,:)
sf0 = sf(1,:,:,:) - sf(0,:,:,:)

Mth = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
mth = Mth
pfx = pf0
scx = sc0
sdx = sd0
sfx = sf0

printVarSummary(pfx)
do i = 0,6
    mth(i) = Mth(i+5)
    pfx(i,:,:) = pf0(i+5,:,:)
    scx(i,:,:) = sc0(i+5,:,:)
    sdx(i,:,:) = sd0(i+5,:,:)
    sfx(i,:,:) = sf0(i+5,:,:)
end do
do i = 0,4
    mth(i+7) = Mth(i)
    pfx(i+7,:,:) = pf0(i,:,:)
    scx(i+7,:,:) = sc0(i,:,:)
    sdx(i+7,:,:) = sd0(i,:,:)
    sfx(i+7,:,:) = sf0(i,:,:)
end do
printVarSummary(pfx)
print(mth)

rfx = pfx
rfx = pfx - sfx
;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)

;=========================================================
; plot 
;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/a20200207/"
figname = "Clim_diff_precip_snowfall_snowcover_snowdepth"

wks = gsn_open_wks("pdf", figpath+figname)
plot = new(6*4, graphic)
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
;res@lbLabelBarOn          = False

res@mpMinLatF            = latstr      ; range to zoom in on
res@mpMaxLatF            = latend
res@mpMinLonF            = lonstr
res@mpMaxLonF            = lonend
; res@mpCenterLonF         = 180.
res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025

res@gsnStringFontHeightF = 0.035
;res@gsnRightString       = "K"

res@pmLabelBarOrthogonalPosF = 0.12
res@pmLabelBarWidthF = 0.7
res@lbLabelFontHeightF = 0.02
res@cnLevelSelectionMode = "ManualLevels"

do i = 0, 5
    res@gsnRightString = mth(i+4)

    ;res@gsnLeftString  = "Precipitation"
    res@gsnLeftString = "Rainfall"
    res@gsnCenterString = pfx@units
    res@cnFillPalette = "precip_diff_12lev"
    
    res@cnMinLevelValF = -5.
    res@cnMaxLevelValF = 5.
    res@cnLevelSpacingF = 0.5
    plot(4*i) = gsn_csm_contour_map(wks, rfx(i+4,:,:), res)
    plot_topo(4*i) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(4*i), plot_topo(4*i))

    res@gsnLeftString = "Snowfall"
    res@gsnCenterString = sfx@units
    res@cnFillPalette = "precip_diff_12lev"
    res@cnMinLevelValF = -0.3
    res@cnMaxLevelValF = 0.3
    res@cnLevelSpacingF = 0.05
    plot(4*i+1) = gsn_csm_contour_map(wks, sfx(i+4,:,:), res)
    plot_topo(4*i+1) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(4*i+1), plot_topo(4*i+1))

    res@gsnLeftString = "Snow Cover"
    res@gsnCenterString = scx@units
    res@cnFillPalette = "MPL_RdBu"
    res@cnMinLevelValF = -0.1
    res@cnMaxLevelValF = 0.1
    res@cnLevelSpacingF = 0.01
    plot(4*i+2) = gsn_csm_contour_map(wks, scx(i+4,:,:), res)
    plot_topo(4*i+2) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(4*i+2), plot_topo(4*i+2))

    res@gsnLeftString = "Snow Depth"
    res@gsnCenterString = sdx@units
    res@cnFillPalette = "MPL_RdBu"
    res@cnMinLevelValF = -20
    res@cnMaxLevelValF = 20
    res@cnLevelSpacingF = 2
    plot(4*i+3) = gsn_csm_contour_map(wks, sdx(i+4,:,:), res)
    plot_topo(4*i+3) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(4*i+3), plot_topo(4*i+3))
end do

pres = True
pres@gsnMaximize = True

gsn_panel(wks, plot, (/6,4/), pres)







; do i = 0, 11
;     res@gsnRightString = mth(i)

;     ;res@gsnLeftString  = "Precipitation"
;     res@gsnLeftString = "Rainfall"
;     res@gsnCenterString = pfx@units
;     res@cnFillPalette = "precip3_16lev"
    
;     res@cnMinLevelValF = 1.
;     res@cnMaxLevelValF = 8.
;     res@cnLevelSpacingF = 0.5
;     plot(4*i) = gsn_csm_contour_map(wks, rfx(i,:,:), res)
;     plot_topo(4*i) = gsn_csm_contour(wks, topo, res_topo)
;     overlay(plot(4*i), plot_topo(4*i))

;     res@gsnLeftString = "Snowfall"
;     res@gsnCenterString = sfx@units
;     res@cnFillPalette = "precip3_16lev"
;     res@cnMinLevelValF = 0.1
;     res@cnMaxLevelValF = 5.
;     res@cnLevelSpacingF = 0.5
;     plot(4*i+1) = gsn_csm_contour_map(wks, sfx(i,:,:), res)
;     plot_topo(4*i+1) = gsn_csm_contour(wks, topo, res_topo)
;     overlay(plot(4*i+1), plot_topo(4*i+1))

;     res@gsnLeftString = "Snow Cover"
;     res@gsnCenterString = scx@units
;     res@cnFillPalette = "WhiteBlue"
;     res@cnMinLevelValF = 0.1
;     res@cnMaxLevelValF = 1.
;     res@cnLevelSpacingF = 0.1
;     plot(4*i+2) = gsn_csm_contour_map(wks, scx(i,:,:), res)
;     plot_topo(4*i+2) = gsn_csm_contour(wks, topo, res_topo)
;     overlay(plot(4*i+2), plot_topo(4*i+2))

;     res@gsnLeftString = "Snow Depth"
;     res@gsnCenterString = sdx@units
;     res@cnFillPalette = "precip3_16lev"
;     res@cnMinLevelValF = 30
;     res@cnMaxLevelValF = 300
;     res@cnLevelSpacingF = 30
;     plot(4*i+3) = gsn_csm_contour_map(wks, sdx(i,:,:), res)
;     plot_topo(4*i+3) = gsn_csm_contour(wks, topo, res_topo)
;     overlay(plot(4*i+3), plot_topo(4*i+3))
; end do

; pres = True
; pres@gsnMaximize = True

; gsn_panel(wks, plot, (/12,4/), pres)


end