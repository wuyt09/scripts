
begin

casename = (/"CTRL","NOQ1","NOTP"/) ;
ncase = dimsizes(casename)

nlat = 96
nlon = 144

do i = 1, ncase-1
    diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(i)+"_F2000_F19G16_CAM4/pro/daily/"
    tname = systemfunc("ls "+diri+"F2000_"+casename(i)+".plevel.T.*.daily.nc")    
    qname = systemfunc("ls "+diri+"F2000_"+casename(i)+".plevel.Q.*.daily.nc")
    print(tname)

    ft = addfiles(tname, "r")
    fq = addfiles(qname, "r")

    T = ft[:]->T
    t = reshape(T, (/25,365,37,nlat,nlon/))
    t_clm = dim_avg_n_Wrap(t, 0)
    copy_VarCoords(T(0:364,:,:,:), t_clm)

    Q = fq[:]->Q
    q = reshape(Q, (/25,365,37,nlat,nlon/))
    q_clm = dim_avg_n_Wrap(q, 0)
    copy_VarCoords(Q(0:364,:,:,:), q_clm)

    nlvl = dimsizes(T&lev)
    nlat = dimsizes(T&lat)
    nlon = dimsizes(T&lon)
    
    lev = T&lev
    levx = conform_dims(dimsizes(t_clm), lev, (/1/))
    w_clm = mixhum_convert(q_clm, "q", (/0,0/))
    tlcl = tlcl_mixr_bolton(t_clm,w_clm,levx,(/1,0,0,1/))
    printVarSummary(tlcl)
    ; thetase = pot_temp_equiv(lev, t_clm, q_clm, 1, "q")
    thetase = pot_temp_equiv_tlcl(levx,t_clm,tlcl,w_clm,(/0,1,0,1/))
    delete([/tname,qname,ft,fq,T,t,t_clm,Q,q,q_clm,lev,levx,w_clm,tlcl/])
    ;-------------------------------------------------------------

    diriout = "/home/yangsong3/wuyt/sysu/scripts/TP_CESM_F/ncl/20200513/"
    filiout = "F2000_"+casename(i)+".plevel.ThetaSe.daily.clm.v2.nc"
    system("rm "+diriout+filiout)
    fout = addfile(diriout+filiout, "c")

        ;Set All field attribution
        fileAtt     = True
        fileAtt@creation_date=systemfunc("date")
        fileattdef(fout,fileAtt)

    ;Define Coordinate
        dimNames    =(/"time","lev","lat","lon"/)
        dimSizes    =(/-1,nlvl,nlat,nlon/)
        dimUnlim    =(/True,False,False,False/)
        filedimdef(fout,dimNames,dimSizes,dimUnlim)

    ;Define var, type and dim
        filevardef(fout,"time",typeof(thetase&time),getvardims(thetase&time))
        filevardef(fout,"lev",typeof(thetase&lev),getvardims(thetase&lev))
        filevardef(fout,"lat",typeof(thetase&lat),getvardims(thetase&lat))
        filevardef(fout,"lon",typeof(thetase&lon),getvardims(thetase&lon))
        filevardef(fout,"thetase",typeof(thetase),getvardims(thetase))

    ;Define Attribute
        filevarattdef(fout,"thetase",thetase)
        filevarattdef(fout,"time",thetase&time)
        filevarattdef(fout,"lev",thetase&lev)
        filevarattdef(fout,"lat",thetase&lat)
        filevarattdef(fout,"lon",thetase&lon)
        
        fout->time=(/thetase&time/)
        fout->lev=(/thetase&lev/)
        fout->lat=(/thetase&lat/)
        fout->lon=(/thetase&lon/)
        fout->thetase=(/thetase/)
end do

end