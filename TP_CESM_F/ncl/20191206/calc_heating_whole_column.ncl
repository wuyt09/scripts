
begin

diri = "/home/yangsong3/data-model/wuyt/TP_CTRL_F2000_F19G16_CAM4/pro/daily/"

filips = "F2000_CTRL.PS.daily.06-30.nc"
fps = addfile(diri+filips, "r")
PS = fps->PS
printVarSummary(PS)

q1 = new((/25,365,96,144/), float)
Cpair = 1004.07
G = 9.8 

do i = 0, 24
    yearx = sprintf("%2i", (i+6))
    filia = "F2000_CTRL.plevel.DTCOND."+yearx+".daily.nc"
    filib = "F2000_CTRL.plevel.DTV."+yearx+".daily.nc"
    filic = "F2000_CTRL.plevel.QRS."+yearx+".daily.nc"
    filid = "F2000_CTRL.plevel.QRL."+yearx+".daily.nc"
    fa = addfile(diri+filia, "r")
    fb = addfile(diri+filib, "r")
    fc = addfile(diri+filic, "r")
    fd = addfile(diri+filid, "r")
    DTCOND = fa->DTCOND
    DTV = fb->DTV
    QRS = fc->QRS
    QRL = fd->QRL
    Q1 = DTCOND
    ;Q1 = DTCOND+DTV+QRS+QRL
    Q1 = DTCOND+DTV
    ;===========================================
    ; whole column 
    ;===========================================
    Q1 = Q1*Cpair/G
    plev = Q1&lev
    plev = plev*100
    plev@units = "Pa"
    psfc = PS(i,:,:,:)
    ptop = min(plev)
    dp = dpres_plevel_Wrap(plev, psfc, ptop, 0)

    Qdp = Q1
    Qdp = Q1*dp

    q1(i,:,:,:) = dim_sum_n_Wrap(Qdp, 1)
    copy_VarCoords(Q1(:,0,:,:), q1(i,:,:,:))

delete([/yearx,Q1,plev,psfc,ptop,dp,Qdp/])
end do

year = ispan(6, 30, 1)
q1!0 = "year"
q1&year = year
q1@units = "W/m^2"
q1@long_name = "Daily Q1 of whole column (sum of 4 terms): derived by CTRL run"
printVarSummary(q1)
g_nlat = dimsizes(q1&lat)
g_nlon = dimsizes(q1&lon)
g_nyr  = dimsizes(year)

filiout = "./heating.daily.4terms.6-30.nc"
system("rm -rf "+filiout)
fout = addfile(filiout, "c")
setfileoption(fout, "DefineMode", True)

;Set All field attribution
fileAtt     = True
fileAtt@creation_date=systemfunc("date")
fileattdef(fout,fileAtt)

;Define Coordinate
dimNames    =(/"year","time","lat","lon"/)
dimSizes    =(/g_nyr,365,g_nlat,g_nlon/)
dimUnlim    =(/False,False,False,False/)
filedimdef(fout,dimNames,dimSizes,dimUnlim)

;Define var, type and dim
filevardef(fout,"year",typeof(q1&year),getvardims(q1&year))
filevardef(fout,"time",typeof(q1&time),getvardims(q1&time))
filevardef(fout,"lat",typeof(q1&lat),getvardims(q1&lat))
filevardef(fout,"lon",typeof(q1&lon),getvardims(q1&lon))
filevardef(fout,"q1",typeof(q1),getvardims(q1))

;Define Attribute
filevarattdef(fout,"q1",q1)
filevarattdef(fout,"year",q1&year)
filevarattdef(fout,"time",q1&time)
filevarattdef(fout,"lat",q1&lat)
filevarattdef(fout,"lon",q1&lon)

fout->year=(/q1&year/)
fout->time=(/q1&time/)
fout->lat=(/q1&lat/)
fout->lon=(/q1&lon/)
fout->q1=(/q1/)

end