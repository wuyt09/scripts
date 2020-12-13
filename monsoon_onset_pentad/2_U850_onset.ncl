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

nlat = 181
nlon = 360
nlev = 19


f = addfile("1_U850_pentad_ERA5_1979_2017.nc", "r")
u = f->u850
printVarSummary(u)

diff = u
last = u
do n = 1, npenta-2
    diff(:,n,:,:) = u(:,n-1,:,:)*u(:,n,:,:)    
    last(:,n,:,:) = u(:,n+1,:,:)*u(:,n,:,:)
end do
printVarSummary(diff)

peak = u(:,0,:,:)
do yy = 0, nyear-1
    do ilat = 0, nlat-1
        do ilon = 0, nlon-1
            a = ind(u(yy,:,ilat,ilon).gt.0 .and. diff(yy,:,ilat,ilon).lt.0 .and. last(yy,:,ilat,ilon).gt.0)
            peak(yy,ilat,ilon) = (/a(0)/)
            delete(a)
        end do
    end do
end do
printVarSummary(peak)

bob_lat1 = 5
bob_lat2 = 15
bob_lon1 = 85
bob_lon2 = 100
print(peak(0,{bob_lat1:bob_lat2},{bob_lon1:bob_lon2}))
printVarSummary(peak)


peak@units = "penta"
peak@long_name = "the onset pentad define by U850 reverse"

outpath = "/home/yangsong3/wuyt/sysu/scripts/monsoon_onset_pentad/"
outname = "2_U850_onset_pentad.nc"

system("rm -rf "+outpath+outname)

fout = addfile(outpath+outname, "c")
setfileoption(fout, "DefineMode", True)

fout->u850_peak = peak

end