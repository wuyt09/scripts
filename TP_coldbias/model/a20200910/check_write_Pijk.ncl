
begin

lonstr = 75
lonend = 105
latstr = 25
latend = 40
yrStrt = 198102
yrLast = 200601


casename = (/"Ctrl_FAMIP","Hist_FAMIP"/) ;,"Hist_IO","Hist_PO","Hist_AO","Hist_PIO","Hist_PAO","Hist_IAO"

ncase = dimsizes(casename)
season = "DJF"

var = (/"hyam","hybm","PS"/) ;"G0=netS-netR-LH-SH" "DownL"
nvar = dimsizes(var)
nlev = 26
nlat = 96
nlon = 144
P0 = 1.e+5
P0@units = "Pa"

i = 0
P = new((/300,nlev,nlat,nlon/), float) ;, "1e+36"
; do i = 0, ncase-1
    diri = "/home/yangsong3/data-model/wuyt/TPbias/TPbias_"+casename(i)+"/a20191206/"

    fili1 = casename(i)+".cam.h0.hyam.1979-2005.nc"
    fili2 = casename(i)+".cam.h0.hybm.1979-2005.nc"
    fili3 = casename(i)+".cam.h0.PS.1979-2005.nc"

    f1 = addfile(diri+fili1, "r")
    f2 = addfile(diri+fili2, "r")
    f3 = addfile(diri+fili3, "r")
    TIME = f1->time
    YYYY = cd_calendar(TIME, -1)
    iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

    hyam = tofloat(f1->hyam(iYYYY,:))
    hybm = tofloat(f2->hybm(iYYYY,:))
    ps = f3->PS(iYYYY,:,:)
    printVarSummary(hyam)
    printVarSummary(hybm)
    printVarSummary(ps)

    do ilev = 0, nlev-1
        do ilat = 0, nlat-1
            do ilon = 0, nlon-1
                P(:,ilev,ilat,ilon) = hyam(:,ilev)*P0+hybm(:,ilev)*ps(:,ilat,ilon)
            end do
        end do
    end do
    copy_VarCoords(ps, P(:,0,:,:))
    copy_VarCoords(hyam(0,:), P(0,:,0,0))
    printVarSummary(P)

    P_season = month_to_season(P, season)
    p_season = dim_avg_n_Wrap(P_season, 0)
    p_season = p_season/100.
    printVarSummary(p_season)
; end do
outpath = "/home/yangsong3/wuyt/sysu/scripts/TP_coldbias/model/a20200910/cfram20200917/data/"
outname = "P_3D_Ctrl.dat"
system("rm -f "+outpath+outname)
fbindirwrite(outpath+outname, p_season)

;=========================================================
; outline of tp
;=========================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)

;=========================================================
; plot
;=========================================================

; p_fig = p_season/100.
; copy_VarCoords(p_season, p_fig)
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/a20200910/"
figname = "check_write_Pijk"


wks = gsn_open_wks("pdf", figpath+figname)
plot = new(nlev, graphic)
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
res@cnFillPalette         = "CBR_coldhot"
;res@lbLabelBarOn          = False

res@mpMinLatF            = latstr      ; range to zoom in on
res@mpMaxLatF            = latend
res@mpMinLonF            = lonstr
res@mpMaxLonF            = lonend
; res@mpCenterLonF         = 180.
res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025

res@gsnStringFontHeightF = 0.035


; res@gsnRightString = "Diff"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = 50
res@cnMaxLevelValF = 600
res@cnLevelSpacingF = 50
do i = 0,nlev-1
    ; res@gsnLeftString = var(i)
    plot(i) = gsn_csm_contour_map(wks, p_season(i,:,:), res)
    plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(i), plot_topo(i))
end do

pres = True
pres@gsnMaximize = True
;pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = season+" mean pressure at diff levels" ;(500hPa to 200hPa)
; pres@pmLabelBarOrthogonalPosF = -0.02
; pres@lbLabelFontHeightF = 0.015

gsn_panel(wks, plot, (/3,9/), pres)


end 





