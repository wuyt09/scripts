
begin

case = "BCC-CSM2-MR"

yrStrt = 1979
yrLast = 2014
year = ispan(yrStrt, yrLast, 1)
season = (/"JJA"/)

varname1 = (/"cl","clw","cli"/)
nvar1 = dimsizes(varname1)

varname2 = (/"ta","hus"/)
nvar2 = dimsizes(varname2)

outname1 = (/"cc","clwc","ciwc"/)
outname2 = (/"ta","hus"/)
outpath = "./hist_79-14/"

diri = "/home/yangsong3/CMIP6/"+case+"/historical/"

base_s = ispan(1979, 2004, 1)
base_e = ispan(1988, 2013, 1)
warm_s = ispan(1980, 2005, 1)
warm_e = ispan(1989, 2014, 1)
n = dimsizes(base_s)        
print("n = "+n)

do i = 0, nvar1-1
    basename = outname1(i)+"_base.dat"
    warmname = outname1(i)+"_warm.dat"

    filic = systemfunc("ls "+diri+varname1(i)+"_Amon_"+case+"_historical_r1i1p1f1_g*.nc")
    nfilic = dimsizes(filic)
    if (nfilic.ge.2) then
        fc = addfiles(filic, "r")
        TIMEc = fc[:]->time
        YYYYc = cd_calendar(TIMEc, -1)/100
        iYYYYc = ind(YYYYc.ge.yrStrt .and. YYYYc.le.yrLast)
        CL = short2flt(fc[:]->$varname1(i)$(iYYYYc,:,:,:))
        ps = short2flt(fc[:]->ps(iYYYYc,:,:))
        hyao = tofloat(fc[0]->a)
        hybo = tofloat(fc[0]->b)
    else
        fc = addfile(filic, "r")
        TIMEc = fc->time
        YYYYc = cd_calendar(TIMEc, -1)/100
        iYYYYc = ind(YYYYc.ge.yrStrt .and. YYYYc.le.yrLast)
        CL = short2flt(fc->$varname1(i)$(iYYYY,:,:,:))
        ps = short2flt(fc->ps(iYYYY,:,:))
        hyao = tofloat(fc->a)
        hybo = tofloat(fc->b)
    end if
    printVarSummary(CL)

    x = month_to_season(CL, season)
    printVarSummary(x)
    x!0 = "year"
    x&year = year

    nlat = dimsizes(x&lat)
    nlon = dimsizes(x&lon)
    nlev = dimsizes(x&lev)

    X_b = new((/n,10,nlev,nlat,nlon/), float)
    X_w = X_b

    do ll = 0, n-1 
        X_b(ll,:,:,:,:) = x({year|base_s(ll):base_e(ll)},lev|:,lat|:,lon|:)
        X_w(ll,:,:,:,:) = x({year|warm_s(ll):warm_e(ll)},lev|:,lat|:,lon|:)
    end do

    xb = dim_avg_n_Wrap(X_b, 1)
    xw = dim_avg_n_Wrap(X_w, 1)
    copy_VarCoords(x(0:n-1,:,:,:), xb)
    copy_VarCoords(x(0:n-1,:,:,:), xw)
    printVarSummary(xb)

    system("rm -f "+outpath+basename)
    fbindirwrite(outpath+basename,xb)
    system("rm -f "+outpath+warmname)
    fbindirwrite(outpath+warmname,xw)
end do

p0 = 100000 ;Pa

;============================================================
; pres 2 hybrid
;============================================================
do i = 0, nvar2-1
    basename = outname2(i)+"_base.dat"
    warmname = outname2(i)+"_warm.dat"

    fili = systemfunc("ls "+diri+varname2(i)+"_Amon_"+case+"_historical_r1i1p1f1_g*.nc")
    nfili = dimsizes(fili)

    if (nfili.ge.2) then
        f1 = addfiles(fili, "r")
        TIME = f1[:]->time
        YYYY = cd_calendar(TIME, -1)/100
        iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
        X0 = short2flt(f1[:]->$varname2(i)$(iYYYY,:,:,:))
    else
        f1 = addfile(fili, "r")
        TIME = f1->time
        YYYY = cd_calendar(TIME, -1)/100
        iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
        X0 = short2flt(f1->$varname2(i)$(iYYYY,:,:,:))
    end if
    printVarSummary(X0)

    plev = X0&plev

    X = pres2hybrid(plev, ps, p0, X0, hyao, hybo, 0)
    copy_VarCoords(CL, X)
    printVarSummary(X)

    x = month_to_season(X, season)
    printVarSummary(x)
    x!0 = "year"
    x&year = year

    nlat = dimsizes(x&lat)
    nlon = dimsizes(x&lon)
    nlev = dimsizes(x&lev)

    X_b = new((/n,10,nlev,nlat,nlon/), float)
    X_w = X_b

    do ll = 0, n-1 
        X_b(ll,:,:,:,:) = x({year|base_s(ll):base_e(ll)},lev|:,lat|:,lon|:)
        X_w(ll,:,:,:,:) = x({year|warm_s(ll):warm_e(ll)},lev|:,lat|:,lon|:)
    end do

    xb = dim_avg_n_Wrap(X_b, 1)
    xw = dim_avg_n_Wrap(X_w, 1)
    copy_VarCoords(x(0:n-1,:,:,:), xb)
    copy_VarCoords(x(0:n-1,:,:,:), xw)
    printVarSummary(xb)

    system("rm -f "+outpath+basename)
    fbindirwrite(outpath+basename,xb)
    system("rm -f "+outpath+warmname)
    fbindirwrite(outpath+warmname,xw)
end do


end