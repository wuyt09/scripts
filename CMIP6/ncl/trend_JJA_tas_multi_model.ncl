begin

lonstr = -20
lonend = 160
latstr = 0
latend = 70

yrStrt = 1979
yrLast = 2014
year = ispan(yrStrt, yrLast, 1)
season = "JJA"

;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)
printVarSummary(topo)
;=========================================================================
diri = "/home/yangsong3/data-observation/ERA5-monthly/surface/"
fili = "t2m_mon_1x1_1979_2019.nc"
f = addfile(diri+fili, "r")
lat = f->latitude(::-1)
lon = f->longitude
nlat = dimsizes(lat)
nlon = dimsizes(lon)
;=========================================================================
model = (/"ACCESS-CM2","ACCESS-ESM1-5","BCC-CSM2-MR","CAS-ESM2-0","CESM2-WACCM","CIESM","CMCC-CM2-SR5","CMCC-ESM2","CanESM5",\
         "E3SM-1-1","FGOALS-f3-L","FGOALS-g3","GFDL-ESM4","IPSL-CM6A-LR","KACE-1-0-G","MIROC6","MPI-ESM1-2-HR","MPI-ESM1-2-LR",\
         "MRI-ESM2-0","NorESM2-LM","NorESM2-MM","TaiESM1"/)
nmodel = dimsizes(model)

; x = new((/nmodel,nyear,nlat,nlon/), float)
rc = new((/nmodel,nlat,nlon/), float)
rc!0 = "model"
rc&model = model
prob = rc
do i = 0, 2 ;nmodel-1
    diri = "/home/yangsong3/CMIP6/"+model(i)+"/historical/"
    fili = systemfunc("ls "+diri+"tas_Amon_"+model(i)+"_historical_r1i1p1f1_g*.nc")
    print(fili)
    nfili = dimsizes(fili)

    if (nfili.ge.2) then
        f1 = addfiles(fili, "r")
        TIME = f1[:]->time
        YYYY = cd_calendar(TIME, -1)/100
        iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
        X = f1[:]->tas(iYYYY,:,:)
    else
        f1 = addfile(fili, "r")
        TIME = f1->time
        YYYY = cd_calendar(TIME, -1)/100
        iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
        X = f1->tas(iYYYY,:,:)
    end if
    printVarSummary(X)
    xlat = X&lat
    xlon = X&lon

    x_med = linint2_Wrap(xlon, xlat, X, True, lon, lat, 0)
    printVarSummary(x_med)
    x_season = month_to_season(x_med, season)

    rc_med = regCoef_n(year, x_season, 0, 0)
    rc(i,:,:) = rc_med*10
    copy_VarCoords(x_season(0,:,:), rc(i,:,:))

    tval = onedtond(rc_med@tval, dimsizes(rc_med))
    df   = onedtond(rc_med@nptxy, dimsizes(rc_med))
    b    = tval
    b    = 0.5
    prob(i,:,:) = betainc(df/(df+tval^2), df/2.0, b)
    copy_VarCoords(x_season(0,:,:), prob(i,:,:))
    delete([/X,xlat,xlon,x_med,x_season/])
end do
printVarSummary(rc)
printVarSummary(prob)

sig = 0.05

;================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/CMIP6/ncl/"
figname = "trend_JJA_tas_multi_model"
wks = gsn_open_wks("pdf", figpath+figname)

nplot = nmodel+1
plot = new(nplot, graphic)
plot_topo = plot
plot_prob = plot
dum1 = plot
dum2 = plot
dum3 = plot
dum4 = plot

 res_topo            = True
 res_topo@cnFillOn   = False
 res_topo@cnLinesOn  = True
 res_topo@gsnFrame   = False
 res_topo@gsnDraw    = False
 res_topo@cnLineLabelsOn         = False
 res_topo@cnLineThicknessF       = 5
 res_topo@cnLineColor            = "forestgreen"
 res_topo@cnInfoLabelOn          = False
 res_topo@cnLevelSelectionMode   ="ExplicitLevels"
 res_topo@cnLevels               = (/2000/)
 res_topo@gsnLeftString = ""
 res_topo@gsnRightString = ""


;----------------------------------------------
; draw prob
;----------------------------------------------
 res3  = True
 res3@gsnDraw   = False
 res3@gsnFrame  = False

 ;res3@cnFillOn  = False 
 res3@cnLinesOn = False
 res3@cnLineLabelsOn = False
 res3@cnInfoLabelOn  = False
 ;res3@lbLabelBarOn   = False

 res3@gsnRightString = ""

 res3@cnLevelSelectionMode = "ManualLevels"
 res3@cnMinLevelValF = 0.00
 res3@cnMaxLevelValF = 0.10
 res3@cnLevelSpacingF = 0.01

 ; res3@cnFillPattern = 17 ;(/17/)
 res3@cnFillColor  =  "black"
 res3@gsnLeftString = ""
 res3@cnFillDotSizeF = 0.003
 res3@cnFillScaleF  = 1.2
 ;----------------------------------------------
 ;----------------------------------------------

res = True
res@gsnDraw = False
res@gsnFrame = False
res@gsnMaximize = True

res@mpOutlineOn = True
res@mpFillOn = True
res@mpGeophysicalLineThicknessF = 5
res@mpGeophysicalLineColor = "grey25" ;"red" 
; res@mpFillDrawOrder = "PostDraw"
; res@mpCenterLonF = 180
res@mpMinLatF = latstr
res@mpMaxLatF = latend
res@mpMinLonF = lonstr
res@mpMaxLonF = lonend


res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False
; res@cnFillMode = "CellFill"

res@gsnStringFontHeightF = 0.025
res@tmXBLabelFontHeightF = 0.02
res@tmYLLabelFontHeightF = 0.02

res@cnFillPalette = "CBR_coldhot"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -1
res@cnMaxLevelValF = 1
res@cnLevelSpacingF = 0.2
res@lbLabelBarOn = False

res@pmLabelBarOrthogonalPosF = 0.12
res@lbLabelFontHeightF = 0.015
; res@tiMainString = "surface temperature"

;-----------------------------------
gres = True
gres@gsLineColor = "blue"
gres@gsLineThicknessF = 5.0
gres@gsLineDashPattern = 16
 glon1 = (/120,120,85,85,120/)
 glat1 = (/45,60,60,45,45/) ; north
 glon2 = (/140,140,105,105,140/)
 glat2 = (/30,45,45,30,30/) ; east
 ; glon3 = (/100,100,50,50,100/)
 ; glat3 = (/10,25,25,10,10/) ; south
 glon3 = (/90,90,70,70,90/)
 glat3 = (/5,25,25,5,5/) ; south
 glon4 = (/60,60,-10,-10,60/)
 glat4 = (/30,55,55,30,30/) ; west
;-----------------------------------

do i = 0, nmodel-1
    res@gsnLeftString = model(i)
    plot(i) = gsn_csm_contour_map(wks, rc(i,:,:), res)
    plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(i), plot_topo(i))

    plot_prob(i) = gsn_csm_contour(wks,prob(i,:,:),res3)
    plot_prob(i) = ShadeLtContour(plot_prob(i),sig,17)
    overlay(plot(i), plot_prob(i))

    dum1(i) = gsn_add_polyline(wks, plot(i), glon1, glat1, gres)
    dum2(i) = gsn_add_polyline(wks, plot(i), glon2, glat2, gres)
    dum3(i) = gsn_add_polyline(wks, plot(i), glon3, glat3, gres)
    dum4(i) = gsn_add_polyline(wks, plot(i), glon4, glat4, gres)
end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
pres@pmLabelBarOrthogonalPosF = 0.2
pres@gsnPanelMainString = season+" trend ("+yrStrt+"-"+yrLast+")"
gsn_panel(wks, plot, (/4,6/), pres)


end