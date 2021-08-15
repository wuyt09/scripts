begin

lonstr = 68
lonend = 107
latstr = 7
latend = 43

yrStrt = 1951
yrLast = 2014
year = ispan(yrStrt, yrLast, 1)
nyear = dimsizes(year)

LEV = 200
;=======================================================
diriobs = "/home/yangsong3/data-observation/GPCP/"
filiobs = "precip.mon.mean.nc"
fobs = addfile(diriobs+filiobs, "r")
lat = fobs->lat({latstr:latend})
lon = fobs->lon({lonstr:lonend})
nlat = dimsizes(lat)
nlon = dimsizes(lon)
;=======================================================================
; dirigcm = "/home/yangsong3/wuyt/sysu/scripts/model_precip/CMIP6/"


; models := (/"CAMS-CSM1-0","CanESM5","CESM2","CESM2-WACCM","CNRM-CM6-1","CNRM-ESM2-1",\
;            "E3SM-1-0","EC-Earth3","EC-Earth3-Veg","FGOALS-f3-L","FGOALS-g3","GISS-E2-1-G-CC","GISS-E2-1-G",\
;            "HadGEM3-GC31-LL","IPSL-CM6A-LR","MPI-ESM1-2-HR","NESM3","SAM0-UNICON","UKESM1-0-LL"/)
models = (/"BCC-CSM2-MR","BCC-ESM1","CAMS-CSM1-0","CanESM5","CESM2-FV2","CESM2","CESM2-WACCM-FV2",\
          "CESM2-WACCM","CIESM","E3SM-1-0","E3SM-1-1-ECA","E3SM-1-1","EC-Earth3","EC-Earth3-Veg",\
          "FGOALS-f3-L","FGOALS-g3","GISS-E2-1-G-CC","GISS-E2-1-G","GISS-E2-1-H","IPSL-CM6A-LR",\
          "KACE-1-0-G","MCM-UA-1-0","MIROC6","MPI-ESM1-2-HR","MRI-ESM2-0","NESM3","SAM0-UNICON"/)

nm = dimsizes(models)

do i = 0, nm-1
    dirigcm = "/home/yangsong3/CMIP6/"+models(i)+"/historical/"
    filiu := systemfunc("ls "+dirigcm+"ua_Amon_"+models(i)+"_historical_r1i1p1f1_g*.nc")
    filiv := systemfunc("ls "+dirigcm+"va_Amon_"+models(i)+"_historical_r1i1p1f1_g*.nc")
    numu = dimsizes(filiu)
    print(numu)

    if (numu.lt.2) then
        fu := addfile(filiu, "r")
        fv := addfile(filiv, "r")

        TIME := fu->time
        YYYY := cd_calendar(TIME, -1)/100
        iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
        U := fu->ua(iYYYY,{LEV},{latstr:latend},{lonstr:lonend})
        V := fv->va(iYYYY,{LEV},{latstr:latend},{lonstr:lonend})

        u := (U(6::12,:,:)+U(7::12,:,:))/2
        v := (V(6::12,:,:)+V(7::12,:,:))/2
        copy_VarCoords(U(0,:,:), u(0,:,:))
        copy_VarCoords(V(0,:,:), v(0,:,:))

    else
        fu := addfiles(filiu, "r")
        fv := addfiles(filiv, "r")

        TIME := fu[:]->time
        YYYY := cd_calendar(TIME, -1)/100
        iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
        U := fu[:]->ua(iYYYY,{LEV},{latstr:latend},{lonstr:lonend})
        V := fv[:]->va(iYYYY,{LEV},{latstr:latend},{lonstr:lonend})

        u := (U(6::12,:,:)+U(7::12,:,:))/2
        v := (V(6::12,:,:)+V(7::12,:,:))/2
        copy_VarCoords(U(0,:,:), u(0,:,:))
        copy_VarCoords(V(0,:,:), v(0,:,:))
    end if
    ; printVarSummary(u)

    





    
    print("end of "+models(i))
end do




end