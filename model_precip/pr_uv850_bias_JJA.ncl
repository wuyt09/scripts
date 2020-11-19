
begin


season = "JJA"

;============================================
; GPCP
;============================================

yrStrt = 1979
yrLast = 2014

diri = "/home/yangsong3/data-observation/GPCP/"
fili = "precip.mon.mean.nc"
f = addfile(diri+fili, "r")

TIME = f->time
YYYY = cd_calendar(TIME, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

precip = f->precip(iYYYY,:,:)

prec = month_to_season(precip, season)
prec_ave = dim_avg_n_Wrap(prec, 0)
prec_var = dim_variance_n_Wrap(prec, 0)
printVarSummary(prec)




;===================================================================
; yrStrt = 1979
; yrLast = 2014
; year = ispan(yrStrt, yrLast, 1)

dirigcm = "/home/yangsong3/wuyt/sysu/scripts/model_precip/CMIP6/"

models = (/"CAMS-CSM1-0","CanESM5","CESM2","CESM2-WACCM","CNRM-CM6-1","CNRM-ESM2-1",\
           "E3SM-1-0","EC-Earth3","EC-Earth3-Veg","FGOALS-f3-L","FGOALS-g3","GISS-E2-1-G-CC","GISS-E2-1-G",\
           "HadGEM3-GC31-LL","IPSL-CM6A-LR","MPI-ESM1-2-HR","NESM3","SAM0-UNICON","UKESM1-0-LL"/)

nm := dimsizes(models)
do i = 0, nm-1
    filigcm = systemfunc("ls "+dirigcm+"pr_Amon_"+models(i)+"_historical_r1i1p1f*_g*_*.nc")
    fgcm = addfile(filigcm, "r")
    TIME := fgcm->time
    YYYY := cd_calendar(TIME, -1)/100
    iYYYY := ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
    PR_GCMs = fgcm->pr(iYYYY,:,:)
    PR_GCMs = PR_GCMs * 86400
    PR_GCMx = linint2_Wrap(PR_GCMs&lon, PR_GCMs&lat, PR_GCMs, True, plon, plat, 0)
    pr_gcm = month_to_season(PR_GCMx, season) ;PR_GCMx
    pr_gcm_ave = dim_avg_n_Wrap(pr_gcm, 0)
    pr_gcm_var = dim_variance_n_Wrap(pr_gcm, 0)

    prob_pr(i,:,:) = ttest(pr_obs_ave, pr_obs_var, nyro, pr_gcm_ave, pr_gcm_var, nyr, False, False)
end do

