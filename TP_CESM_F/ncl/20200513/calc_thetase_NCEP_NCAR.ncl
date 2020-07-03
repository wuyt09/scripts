
begin

casename = (/"CTRL","NOQ1","NOTP"/) ;
ncase = dimsizes(casename)

nlat = 96
nlon = 144

LEV = (/1000, 925, 850, 700, 600, 500, 400, 300/)


    diri = "/home/yangsong3/data-observation/NCEP_NCAR/ltmmean/"
    tname = "air.day.1981-2010.ltm.nc"    
    qname = "shum.day.1981-2010.ltm.nc"
    print(tname)

    ft = addfile(diri+tname, "r")
    fq = addfile(diri+qname, "r")

    t_clm = ft->air(:,{LEV},:,:)
    q_clm = fq->shum

    nlvl = dimsizes(t_clm&level)
    nlat = dimsizes(t_clm&lat)
    nlon = dimsizes(t_clm&lon)
    
    lev = t_clm&level
    ; levx = conform_dims(dimsizes(t_clm), lev, (/1/))
    ; w_clm = mixhum_convert(q_clm, "q", (/0,0/))
    ; tlcl = tlcl_mixr_bolton(t_clm,w_clm,levx,(/1,0,0,1/))
    ; printVarSummary(tlcl)
    ; thetase = pot_temp_equiv_tlcl(levx,t_clm,tlcl,w_clm,(/0,1,0,1/))
    ; delete([/tname,qname,ft,fq,t_clm,q_clm,lev,levx,w_clm,tlcl/])


    thetase = pot_temp_equiv(lev, t_clm, q_clm, 1, "q")
    ;-------------------------------------------------------------

    diriout = "/home/yangsong3/wuyt/sysu/scripts/TP_CESM_F/ncl/20200513/"
    filiout = "NCEP_NCAR.ThetaSe.daily.clm.nc" ;.v2
    system("rm "+diriout+filiout)
    fout = addfile(diriout+filiout, "c")

        ;Set All field attribution
        fileAtt     = True
        fileAtt@creation_date=systemfunc("date")
        fileattdef(fout,fileAtt)

    ;Define Coordinate
        dimNames    =(/"time","level","lat","lon"/)
        dimSizes    =(/-1,nlvl,nlat,nlon/)
        dimUnlim    =(/True,False,False,False/)
        filedimdef(fout,dimNames,dimSizes,dimUnlim)

    ;Define var, type and dim
        filevardef(fout,"time",typeof(thetase&time),getvardims(thetase&time))
        filevardef(fout,"level",typeof(thetase&level),getvardims(thetase&level))
        filevardef(fout,"lat",typeof(thetase&lat),getvardims(thetase&lat))
        filevardef(fout,"lon",typeof(thetase&lon),getvardims(thetase&lon))
        filevardef(fout,"thetase",typeof(thetase),getvardims(thetase))

    ;Define Attribute
        filevarattdef(fout,"thetase",thetase)
        filevarattdef(fout,"time",thetase&time)
        filevarattdef(fout,"level",thetase&level)
        filevarattdef(fout,"lat",thetase&lat)
        filevarattdef(fout,"lon",thetase&lon)
        
        fout->time=(/thetase&time/)
        fout->level=(/thetase&level/)
        fout->lat=(/thetase&lat/)
        fout->lon=(/thetase&lon/)
        fout->thetase=(/thetase/)

end