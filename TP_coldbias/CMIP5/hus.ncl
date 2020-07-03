
begin

models = (/"ACCESS1-0","bcc-csm1-1","CanESM2","CCSM4","CESM1-BGC","CESM1-CAM5",\ ;"ACCESS1-3",
        "CMCC-CM","CMCC-CMS","CNRM-CM5","CSIRO-Mk3-6-0","GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M",\
        "GISS-E2-H","GISS-E2-R","inmcm4","IPSL-CM5A-LR","IPSL-CM5A-MR",\
        "IPSL-CM5B-LR","MIROC5","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-M","NorESM1-ME"/) ;"HadGEM2-CC","HadGEM2-ES",

; models = (/"bcc-csm1-1","CanESM2","CCSM4","CESM1-BGC","CESM1-CAM5",\
;         "CMCC-CM","CMCC-CMS","CNRM-CM5","CSIRO-Mk3-6-0","GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M",\
;         "GISS-E2-H","GISS-E2-R","inmcm4","IPSL-CM5A-LR","IPSL-CM5A-MR",\
;         "IPSL-CM5B-LR","MIROC5","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-M","NorESM1-ME"/)

yrStrt = 1951
yrLast = 2000

nm = dimsizes(models)
print(nm)

diriera = "/home/yangsong3/data-observation/wuyt/ERAI/monthly/1x1/"
filiera = "t2m.mon.mean.7901.1804.nc"
fera = addfile(diriera+filiera, "r")
TIMEera = fera->time
yyyyera = cd_calendar(TIMEera, -1)/100
iyyyyera = ind(yyyyera.ge.1979 .and. yyyyera.le.2017)
TS_era = short2flt(fera->t2m(iyyyyera,::-1,:))
ts_era = clmMonTLL(TS_era)
ts_era!1 = "lat"
ts_era!2 = "lon"
;printVarSummary(ts_era)
lat = ts_era&lat;itude
lon = ts_era&lon;gitude

LEV = (/1000, 925, 850, 700, 600, 500, 400, 300/)
LEV  = LEV*100
diri = "/home/yangsong3/wuyt/sysu/scripts/TP_coldbias/CMIP5/"
do i = 0, nm-1
    fili = systemfunc("ls "+diri+"hus_Amon_"+models(i)+"_historical_r1i1p1_18*-20*.nc")
    print(i+"   "+fili)
    f = addfile(fili, "r")
    TIME = f->time
    YYYY = cd_calendar(TIME, -1)/100
    ;print(YYYY)
    iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
    TS = f->hus(iYYYY,{LEV},:,:)
    printVarSummary(TS)
    ; time = cd_calendar(TS&time, -1) ;/100
    ; print(i+"   "+time)
    ; exit
    ;Ts = clmMonTLL(TS)
    ;printVarSummary(Ts)

    ts = linint2_Wrap(TS&lon, TS&lat, TS, True, lon, lat, 0)
    printVarSummary(ts)

    filiout = "./post-pro/hus_mon_"+models(i)+"_historical_r1i1p1_195101-200012.nc"
    system("rm -f "+filiout)
    fout = addfile(filiout, "c")
    fout->hus = ts

    delete([/fili, f, TIME, YYYY, iYYYY, TS, ts, fout/])

    print("end "+models(i))
end do


end
