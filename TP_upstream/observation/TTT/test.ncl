
begin

yrStrt = 1979
yrLast = 2016
year = ispan(yrStrt, yrLast, 1)
ny = dimsizes(year)
Season = (/"MAM"/)
ns = dimsizes(Season)

;---------------------------------------------
path = "./"
name = "TTT_TP_dtrend_standardize_NCEP2.txt"

x = asciiread(path+name, (/ny,ns/), "float")
x@_FillValue = -32767
x!0 = "year"
x&year = year
x!1 = "season"
x&season = Season
printVarSummary(x)

index = x(:,0) ;MAM
indlow = ind(index .le. -0.8)
indhig = ind(index .ge. 0.8)
numhig = dimsizes(indhig)
numlow = dimsizes(indlow)
yearhig = year(indhig)
yearlow = year(indlow)
print(yearhig)
print(yearlow)
;---------------------------------------------

;===================Composite========================
diri = "/home/yangsong3/data-observation/GPCP/"
fili = "precip.mon.mean.nc" 
f = addfile(diri+fili, "r")

TIME = f->time
YYYY = cd_calendar(TIME, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

p = f->precip(iYYYY,:,:)
printVarSummary(p)

;----------------------------------------------
pclm = clmMonTLL(p)
pano = calcMonAnomTLL(p, pclm)
psea = month_to_season(pano, "MAM")
printVarSummary(psea)
psea!0 = "year"
psea&year = year
printVarSummary(psea)
print(psea(4,{25:40},{75:105}))

;-----------------------------------------------
pSea = month_to_season(p, "MAM")
Pclm = dim_avg_n_Wrap(pSea, 0)
PSea = conform_dims(dimsizes(pSea), Pclm, (/1,2/))
Pano = pSea
Pano = pSea - PSea
print(Pano(4,{25:40},{75:105}))

end