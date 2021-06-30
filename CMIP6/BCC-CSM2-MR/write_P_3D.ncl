
begin

case = "BCC-CSM2-MR"

yrStrt = 1979
yrLast = 2014
year = ispan(yrStrt, yrLast, 1)
season = (/"JJA"/)


diri = "/home/yangsong3/CMIP6/"+case+"/historical/"

filip = systemfunc("ls "+diri+"ps_Amon_"+case+"_historical_r1i1p1f1_g*.nc")
print(filip)
nfili = dimsizes(filip)
if (nfili.ge.2) then
    f1 = addfiles(filip, "r")
    TIME = f1[:]->time
    YYYY = cd_calendar(TIME, -1)/100
    iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
    ps = f1[:]->ps(iYYYY,:,:)
else
    f1 = addfile(filip, "r")
    TIME = f1->time
    YYYY = cd_calendar(TIME, -1)/100
    iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
    ps = f1->ps(iYYYY,:,:)
end if
printVarSummary(ps)
delete([/f1,TIME,YYYY,iYYYY/])

filic = systemfunc("ls "+diri+"cl_Amon_"+case+"_historical_r1i1p1f1_g*.nc")
print(filic)
nfilic = dimsizes(filic)
if (nfilic.ge.2) then
    f1 = addfiles(filic, "r")
    hyam = tofloat(f1[0]->a)
    hybm = tofloat(f1[0]->b)
    level = tofloat(f1[0]->lev)
else
    f1 = addfile(filic, "r")
    hyam = tofloat(f1->a)
    hybm = tofloat(f1->b)
    level = tofloat(f1->lev)
end if
printVarSummary(hyam)

ntime = dimsizes(ps&time)
nlev = dimsizes(level)
nlat = dimsizes(ps&lat)
nlon = dimsizes(ps&lon)

P0 = 1.e+5
P0@units = "Pa"
P = new((/ntime,nlev,nlat,nlon/), float) ;, "1e+36"

do ilev = 0, nlev-1
    do ilat = 0, nlat-1
        do ilon = 0, nlon-1
            P(:,ilev,ilat,ilon) = hyam(ilev)*P0+hybm(ilev)*ps(:,ilat,ilon)
        end do
    end do
end do
copy_VarCoords(ps, P(:,0,:,:))
copy_VarCoords(hyam(0,:), P(0,:,0,0))
printVarSummary(P)

end