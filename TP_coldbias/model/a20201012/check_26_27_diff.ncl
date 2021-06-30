
begin

  lonstr = 60 ;10
  lonend = 110 ;150
  latstr = 20 ;-10
  latend = 50 ;70

yrStrt = 198102
yrLast = 200601

nvname = (/"rad_conv","lwu","lwd","swu","swd","lwu_toa","lwd_toa","swu_toa","swd_toa"/)  ; ,"error"
nvars = dimsizes(nvname)

nvname2 = (/"rad_conv","sw","lw","lwu","lwd","swu","swd"/)  ; ,"error"
nvars2 = dimsizes(nvname2)

;================read metadata=============================================

casename = (/"Ctrl_FAMIP","Hist_FAMIP"/) ;,"Hist_IO","Hist_PO","Hist_AO","Hist_PIO","Hist_PAO","Hist_IAO"

ncase = dimsizes(casename)
season = "DJF"

ssrd = new((/ncase,96,144/), float)

do i = 0, ncase-1
    diri = "/home/yangsong3/data-model/wuyt/TPbias/TPbias_"+casename(i)+"/a20191206/"
    fili = casename(i)+".cam.h0.FLDS.1979-2005.nc"
    f = addfile(diri+fili, "r")

    TIME = f->time
    YYYY = cd_calendar(TIME, -1)
    iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

    SSRD = f->FLDS(iYYYY,:,:)
    ssrdx = month_to_season(SSRD, season)
    ssrd(i,:,:) = dim_avg_n_Wrap(ssrdx, 0)
end do
printVarSummary(ssrd)

lon = tofloat(ssrd&lon)
lat = tofloat(ssrd&lat)

nlat = dimsizes(lat)
nlon = dimsizes(lon)
delete(f)
;============================================================================
diri2 = "/home/yangsong3/wuyt/sysu/scripts/TP_coldbias/model/a20200910/cfram20200917/" ;cfram20200805
fili2 = diri2+"baseline_radsfc_ranc_1.grd"
fili3 = diri2+"warm_radsfc_ranc_1.grd"

fili4 = diri2+"baseline_radranc_1.grd"
fili5 = diri2+"warm_radranc_1.grd"

UNDEF = -999.
x = new((/nvars,nlat,nlon/),float,UNDEF)
x!0 = "var"
x!1 = "lat"
x!2 = "lon"
x&var = nvname
x&lat = lat
x&lon = lon
x&lat@units = "degrees_north"
x&lon@units = "degrees_east"
x@long_name = "forcing"
x@units     = "W/m2"

y = x

X = new((/nvars,27,nlat,nlon/),float,UNDEF)
X!0 = "var"
X!1 = "lev"
X!2 = "lat"
X!3 = "lon"
X&var = nvname
X&lat = lat
X&lon = lon
X&lat@units = "degrees_north"
X&lon@units = "degrees_east"
X@long_name = "forcing"
X@units     = "W/m2"
Y = X

do nv = 0,nvars-1
    x(nv,:,:) = fbindirread(fili2,nv,(/nlat,nlon/),"float")
    y(nv,:,:) = fbindirread(fili3,nv,(/nlat,nlon/),"float")
end do

do nv = 0,nvars2-1
    do nl = 0, 26
        X(nv,nl,:,:) = fbindirread(fili4,nv*27+nl,(/nlat,nlon/),"float")
        Y(nv,nl,:,:) = fbindirread(fili5,nv*27+nl,(/nlat,nlon/),"float")
    end do
end do

ilev = 24
diff1 = y(2,:,:) - x(2,:,:)
copy_VarCoords(x(2,:,:), diff1)
diff2 = Y(4,ilev,:,:) - X(4,ilev,:,:)
copy_VarCoords(X(4,ilev,:,:), diff2)
printVarSummary(diff1)

diff = diff2 - diff1
copy_VarCoords(diff1, diff)
printVarSummary(diff)
;=======================================
; this is for outline of TP
;=======================================

ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo ;(::-1,:)
;==========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/a20201012/" ;diri2
figname = "check_26_27_diff"
wks  = gsn_open_wks("pdf",figpath+figname)
gsn_define_colormap(wks,"temp_diff_18lev")


plot = new(3, graphic)
plot_topo = plot

;----------------------
res_topo            = True
res_topo@cnFillOn   = False
res_topo@cnLinesOn  = True
res_topo@gsnFrame   = False
res_topo@gsnDraw    = False
res_topo@cnLineLabelsOn         = False
res_topo@cnLineThicknessF       = 3
res_topo@cnLineColor            = "grey40"
;res_topo@cnLineDashPattern      = "1"
res_topo@cnInfoLabelOn          = False
res_topo@cnLevelSelectionMode   ="ExplicitLevels"
res_topo@cnLevels               = (/1500/)
res_topo@gsnLeftString = ""
res_topo@gsnRightString = ""
;==========================================================

res = True
res@gsnDraw   = False
res@gsnFrame  = False
;res@gsnSpreadColors = True
res@gsnAddCyclic = False

res@mpFillOn  = True
;res@mpMinLonF = 180.
;res@mpMaxLonF = 360.
res@mpMinLatF = latstr
res@mpMaxLatF = latend
res@mpMinLonF = lonstr
res@mpMaxLonF = lonend
;res@mpCenterLonF = 85.
res@mpGeophysicalLineThicknessF = 2

res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False

res@cnLevelSelectionMode = "ManualLevels"
res@cnMaxLevelValF = 10
res@cnMinLevelValF = -10
res@cnLevelSpacingF = 1

; do i = 0, dimsizes(title)-1
;     res@gsnLeftString = title(i)
;     plot(i) = gsn_csm_contour_map(wks, x_diff(i,:,:), res)
;     plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
;     overlay(plot(i), plot_topo(i))
; end do

plot(0) = gsn_csm_contour_map(wks, diff1, res)
plot_topo(0) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(0), plot_topo(0))

plot(1) = gsn_csm_contour_map(wks, diff2, res)
plot_topo(1) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(1), plot_topo(1))

plot(2) = gsn_csm_contour_map(wks, diff, res)
plot_topo(2) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(2), plot_topo(2))

pres = True
pres@gsnMaximize = True
pres@gsnFrame    = True
pres@gsnDraw     = True
pres@gsnPanelLabelBar = True

pres@lbLabelFontHeightF = 0.015
pres@pmLabelBarWidthF = 0.8    ;plot manager 
pres@pmLabelBarHeightF = 0.06

gsn_panel(wks,plot,(/1,3/),pres)

end