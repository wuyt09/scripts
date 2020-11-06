
begin

lonstr = 0
lonend = 360
latstr = -30
latend = 70


yrStrt = 198102
yrLast = 200601


casename = (/"Ctrl_FAMIP","Hist_FAMIP"/) ;,"Hist_FAMIP_CMIP6","Hist_PO""Hist_PO","Hist_IO","Hist_PO","Hist_AO","Hist_PIO","Hist_PAO","Hist_IAO"
ncase = dimsizes(casename)

season = "DJF"
nlat = 96
nlon = 144

LEV = 300


nyear = 25
t = new((/ncase,nyear,96,144/), float)
SF = new((/nyear*12,96,144/), float)
VP = SF

do i = 0, ncase-1
    diri = "/home/yangsong3/data-model/wuyt/TPbias/TPbias_"+casename(i)+"/a20191206/"

    fili1 = casename(i)+".cam.h0.T.1979-2005.nc" 
    fili2 = casename(i)+".cam.h0.V.1979-2005.nc"    
    fili3 = casename(i)+".cam.h0.U.1979-2005.nc"
    f1 = addfile(diri+fili1, "r")
    f2 = addfile(diri+fili2, "r")
    f3 = addfile(diri+fili3, "r")

    TIME = f2->time
    YYYY = cd_calendar(TIME, -1)
    iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

    T = f1->T(iYYYY,{LEV},:,:)
    V = f2->V(iYYYY,{LEV},:,:)
    U = f3->U(iYYYY,{LEV},:,:)

    uv2sfvpg(U, V, SF, VP)
    printVarSummary(SF)
    copy_VarCoords(U, SF)
    copy_VarCoords(U, VP)

    t(i,:,:,:) = month_to_season(SF, season)
    copy_VarCoords(U(0,:,:), t(i,0,:,:))
end do

printVarSummary(t)

