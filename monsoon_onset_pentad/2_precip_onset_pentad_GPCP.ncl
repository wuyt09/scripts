begin

yrStrt = 1979
yrLast = 2018
mnStrt = 1
mnLast = 12

year = ispan(yrStrt, yrLast, 1)
month = ispan(mnStrt, mnLast, 1)
penta = ispan(1, 73, 1)
nyear = dimsizes(year)
nmonth = dimsizes(month)
npenta = dimsizes(penta)

bob_lat1 = 5
bob_lat2 = 15
bob_lon1 = 90
bob_lon2 = 100

scs_lat1 = 10
scs_lat2 = 20
scs_lon1 = 110
scs_lon2 = 120

ism_lat1 = 10
ism_lat2 = 20
ism_lon1 = 60
ism_lon2 = 85


lat = fspan(-88.75, 88.75, 72)
lon = fspan(1.25, 358.75, 144)
lat@units = "degree_north"
lon@units = "degree_east"
nlat = dimsizes(lat)
nlon = dimsizes(lon)

diri = "/home/yangsong3/data-observation/GPCP/pentad/"
prec = new((/nyear,npenta,nlat,nlon/), float)
prec@_FillValue = -999.0
prec!0 = "year"
prec&year = year
prec!1 = "penta"
prec&penta = penta
prec!2 = "lat"
prec&lat = lat
prec!3 = "lon"
prec&lon = lon
do i = 0, nyear-1
    fili = "gpcp_pen_v2.2.lnx."+year(i)
    do p = 0, npenta-1
        prec(i,p,:,:) = fbindirread(diri+fili, p, (/nlat,nlon/), "float")
    end do
end do
printVarSummary(prec)

; print(prec(0,20,{10:20},{100:120}))


;===================================================================
; land sea mask
;===================================================================
lat_new = fspan(-90, 90, 181)
lon_new = fspan(0, 359, 360)
lat_new@units = "degree_north"
lon_new@units = "degree_east"

a      = addfile ("/home/yangsong3/wuyt/data/lsm_1279l4_0.1x0.1.grb_v4_unpack.nc","r")
; lsdata = doubletofloat(a->lsm(0,:,:))
lsdata = a->lsm(0,::-1,:)
printVarSummary(lsdata)
lsm = linint2_Wrap(lsdata&longitude, lsdata&latitude, lsdata, True, lon_new, lat_new, 0)
printVarSummary(lsm)

prec_newx = linint2_Wrap(lon, lat, prec, True, lon_new, lat_new, 0)
; prec_new = linint2_Wrap(lon, lat, prec, True, lsm&lon, lsm&lat, 0)


prec_newx!2 = "lat"
prec_newx!3 = "lon"

pens = 18
prec_new = prec_newx(:,pens:npenta-1,:,:)
nlat = dimsizes(prec_new&lat)
nlon = dimsizes(prec_new&lon)
printVarSummary(prec_new)


prec_bob = wgt_areaave_Wrap(prec_new(:,:,{bob_lat1:bob_lat2},{bob_lon1:bob_lon2}), 1., 1., 0)
prec_scs = wgt_areaave_Wrap(prec_new(:,:,{scs_lat1:scs_lat2},{scs_lon1:scs_lon2}), 1., 1., 0)
prec_ism = wgt_areaave_Wrap(prec_new(:,:,{ism_lat1:ism_lat2},{ism_lon1:ism_lon2}), 1., 1., 0)

peak = prec_new(:,0,:,:)
BOB = new(nyear, integer)
SCS = BOB
ISM = BOB

do yy = 0, nyear-1
    do ilat = 0, nlat-1
        do ilon = 0, nlon-1
            ; a = ind(prec(yy,:,ilat,ilon).gt.0 .and. diff(yy,:,ilat,ilon).lt.0 .and. last(yy,:,ilat,ilon).gt.0)
            if (lsm(ilat,ilon).eq.0) then
                a := ind(prec_new(yy,:,ilat,ilon).gt.5.)
            else
                a := ind(prec_new(yy,:,ilat,ilon).gt.3.)
            end if
            ; printVarSummary(a)
            peak(yy,ilat,ilon) = (/a(0)/) + pens
            delete(a)
        end do
    end do

    a = ind(prec_bob(yy,:).ge.5.)
    BOB(yy) =(/a(0)/) + pens
    delete(a)

    a = ind(prec_scs(yy,:).ge.5.)
    SCS(yy) =(/a(0)/) + pens
    delete(a)

    a = ind(prec_ism(yy,:).ge.5.)
    ISM(yy) =(/a(0)/) + pens
    delete(a)
end do

print(BOB+" "+SCS+" "+ISM)

write_table("prec_GPCP_onset_BOB_SCS_ISM.txt", "w", [/"YEAR","BOB","SCS","ISM"/], "%10s %10s %10s %10s")
write_table("prec_GPCP_onset_BOB_SCS_ISM.txt", "a", [/year,BOB,SCS,ISM/], "%0.4i %0.2i %0.2i %0.2i")


; bob = wgt_areaave_Wrap(peak(:,{bob_lat1:bob_lat2},{bob_lon1:bob_lon2}), 1., 1., 0)
; scs = wgt_areaave_Wrap(peak(:,{scs_lat1:scs_lat2},{scs_lon1:scs_lon2}), 1., 1., 0)
; ism = wgt_areaave_Wrap(peak(:,{ism_lat1:ism_lat2},{ism_lon1:ism_lon2}), 1., 1., 0)
; print(bob+" "+scs+" "+ism)

peak@units = "penta"
peak@long_name = "the onset pentad define by prec"

outpath = "/home/yangsong3/wuyt/sysu/scripts/monsoon_onset_pentad/"
outname = "2_precip_onset_pentad_GPCP_"+yrStrt+"_"+yrLast+".nc"

system("rm -rf "+outpath+outname)

fout = addfile(outpath+outname, "c")
setfileoption(fout, "DefineMode", True)

fout->prec_peak = peak


;================================================================
; plot
;================================================================
lonstr = 50
lonend = 150
latstr = 0
latend = 30

peak_clm = dim_avg_n_Wrap(peak, 0)

figpath = "/home/yangsong3/wuyt/sysu/figures/monsoon_onset_pentad/"
figname = "2_precip_GPCP_onset"
wks = gsn_open_wks("pdf", figpath+figname)
plot = new(1, graphic)


res = True
res@gsnDraw = False
res@gsnFrame = False
res@gsnMaximize = True

res@mpOutlineOn = True
res@mpFillOn = True
res@mpGeophysicalLineThicknessF = 2
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

res@gsnStringFontHeightF = 0.025
res@tmXBLabelFontHeightF = 0.02
res@tmYLLabelFontHeightF = 0.02

res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = 20
res@cnMaxLevelValF = 40
res@cnLevelSpacingF = 2

plot(0) = gsn_csm_contour_map(wks, peak_clm, res)
draw(plot)

end