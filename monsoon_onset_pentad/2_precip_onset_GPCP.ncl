begin

yrStrt = 1979
yrLast = 2017
mnStrt = 1
mnLast = 12

year = ispan(yrStrt, yrLast, 1)
month = ispan(mnStrt, mnLast, 1)
penta = ispan(1, 73, 1)
nyear = dimsizes(year)
nmonth = dimsizes(month)
npenta = dimsizes(penta)

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
do i = 0, 1;nyear-1
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

prec_new = linint2_Wrap(lon, lat, prec, True, lon_new, lat_new, 0)
; prec_new = linint2_Wrap(lon, lat, prec, True, lsm&lon, lsm&lat, 0)
printVarSummary(prec_new)

prec_new!2 = "lat"
prec_new!3 = "lon"
nlat = dimsizes(prec_new&lat)
nlon = dimsizes(prec_new&lon)


; diff = prec
; last = prec
; do n = 1, npenta-2
;     diff(:,n,:,:) = prec(:,n-1,:,:)*prec(:,n,:,:)    
;     last(:,n,:,:) = prec(:,n+1,:,:)*prec(:,n,:,:)
; end do
; printVarSummary(diff)

peak = prec_new(:,0,:,:)
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
            peak(yy,ilat,ilon) = (/a(0)/)
            delete(a)
        end do
    end do
end do

bob_lat1 = 5
bob_lat2 = 15
bob_lon1 = 85
bob_lon2 = 100
print(lsm({bob_lat1:bob_lat2},{bob_lon1:bob_lon2})+"    "+peak(0,{bob_lat1:bob_lat2},{bob_lon1:bob_lon2}))
printVarSummary(peak)


peak@units = "penta"
peak@long_name = "the onset pentad define by prec"

outpath = "/home/yangsong3/wuyt/sysu/scripts/monsoon_onset_pentad/"
outname = "2_precip_onset_pentad_GPCP.nc"

system("rm -rf "+outpath+outname)

fout = addfile(outpath+outname, "c")
setfileoption(fout, "DefineMode", True)

fout->prec_peak = peak

end