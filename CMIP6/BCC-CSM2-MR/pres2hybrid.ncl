
begin

case = "BCC-CSM2-MR"

yrStrt = 1979
yrLast = 2014
year = ispan(yrStrt, yrLast, 1)
season = (/"JJA"/)

varname = (/"ta","hus"/)
nvar = dimsizes(varname)

outname = (/"ta","hus"/)
outpath = "./hist_79-19/"

diri = "/home/yangsong3/CMIP6/BCC-CSM2-MR/historical/"

base_s = ispan(1979, 2004, 1)
base_e = ispan(1988, 2013, 1)
warm_s = ispan(1980, 2005, 1)
warm_e = ispan(1989, 2014, 1)
n = dimsizes(base_s)        
print("n = "+n)


filic = systemfunc("ls "+diri+"cl_Amon_"+case+"_historical_r1i1p1f1_g*.nc")
nfilic = dimsizes(filic)

if (nfilic.ge.2) then
    fc = addfiles(filic, "r")
    TIMEc = fc[:]->time
    YYYYc = cd_calendar(TIMEc, -1)/100
    iYYYYc = ind(YYYYc.ge.yrStrt .and. YYYYc.le.yrLast)
    CL = short2flt(fc[:]->cl(iYYYYc,:,:,:))
    ps = short2flt(fc[:]->ps(iYYYYc,:,:))
    hyao = tofloat(fc[:]->a)
    hybo = tofloat(fc[:]->b)
else
    fc = addfile(filic, "r")
    TIMEc = fc->time
    YYYYc = cd_calendar(TIMEc, -1)/100
    iYYYYc = ind(YYYYc.ge.yrStrt .and. YYYYc.le.yrLast)
    CL = short2flt(fc->cl(iYYYY,:,:,:))
    ps = short2flt(fc->ps(iYYYY,:,:))
    hyao = tofloat(fc->a)
    hybo = tofloat(fc->b)
end if
printVarSummary(CL)
p0 = 100000 ;Pa


do i = 0, 0 ;nvar-1
    fili = systemfunc("ls "+diri+varname(i)+"_Amon_"+case+"_historical_r1i1p1f1_g*.nc")
    nfili = dimsizes(fili)

    if (nfili.ge.2) then
        f1 = addfiles(fili, "r")
        TIME = f1[:]->time
        YYYY = cd_calendar(TIME, -1)/100
        iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
        X0 = short2flt(f1[:]->$varname(i)$(iYYYY,:,:,:))
    else
        f1 = addfile(fili, "r")
        TIME = f1->time
        YYYY = cd_calendar(TIME, -1)/100
        iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
        X0 = short2flt(f1->$varname(i)$(iYYYY,:,:,:))
    end if
    printVarSummary(X0)

    plev = X0&plev

    X = pres2hybrid(plev, ps, p0, X0, hyao, hybo, 0)
    printVarSummary(X)

end do


end