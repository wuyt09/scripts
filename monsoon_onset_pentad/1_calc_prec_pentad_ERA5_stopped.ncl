begin
; this pentad data are quite different with GPCP pentad data, careful


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

mdiri = "/home/yangsong3/data-observation/ERA5-daily/Single-level/prect/"

prec = new((/nyear,npenta,nlat,nlon/), float)

do i = 0, 1 ;nyear-1
    filip = systemfunc("ls "+mdiri+"prect."+year(i)+"*.daily.nc")
    fp = addfiles(filip, "r")    
    PRECX = short2flt(fp[:]->tp(:,::-1,:))
    PRECX = PRECX * 1000. ; convert units from m/day to mm/day

    if (dimsizes(PRECX&time).eq.366) then
        PREC = PRECX(0:364,:,:)
        PREC(0:58,:,:) = PRECX(0:58,:,:)
        PREC(59:364,:,:) = PRECX(60:365,:,:)
    else
        PREC = PRECX
    end if
    ; printVarSummary(U)
    delete([/PRECX/])

    do p = 0, 72
        prec(i,p,:,:) = dim_sum_n_Wrap(PREC(5*p:5*p+4,:,:), 0)/5.0
    end do
    copy_VarCoords(PREC(0,:,:), prec(0,0,:,:))
end do
prec!0 = "year"
prec&year = year
prec!1 = "pentad"
prec&pentad = penta
printVarSummary(prec)


latn = fspan(-88.75, 88.75, 72)
lonn = fspan(1.25, 358.75, 144)
latn@units = "degree_north"
lonn@units = "degree_east"
prec_new = linint2_Wrap(prec&longitude, prec&latitude, prec, True, lonn, latn, 0)
print(prec_new(0,20,{10:20},{100:120}))

; prec@units = "mm/day"
; prec@long_name = "zonal wind at 850 hPa derived from ERA5 daily data"

; outpath = "/home/yangsong3/wuyt/sysu/scripts/monsoon_onset_pentad/"
; outname = "1_U850_pentad_ERA5_1979_2017.nc"

; system("rm -rf "+outpath+outname)

; fout = addfile(outpath+outname, "c")
; setfileoption(fout, "DefineMode", True)

; fout->u850 = u


end
