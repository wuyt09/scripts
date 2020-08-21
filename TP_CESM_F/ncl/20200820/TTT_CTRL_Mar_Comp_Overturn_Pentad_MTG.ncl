
begin

   latS = 0.
   latN = 30.
   lonL = 50.
   lonR = 150.

  latS1  = 10
  latN1  = 20
  lonL1  = 60
  lonR1  = 85.0

  latS2  = 5
  latN2  = 15
  lonL2  = 90
  lonR2  = 100.0

  latS3  = 10
  latN3  = 20
  lonL3  = 110
  lonR3  = 120.0

lonstr = 40
lonend = 150
latstr = -15
latend = 50


casename = (/"CTRL"/)

ncase = dimsizes(casename)
nlev = 37
nlat = 34
nlon = 45

levs = 600
levt = 300
u = new((/ncase,12,nlev,nlat,nlon/), float)
v = u

i = 0
; do i = 0, ncase-1
    diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(i)+"_F2000_F19G16_CAM4/pro/"
    filit = "F2000_"+casename(i)+".plevel.T.mon.0601-3012.nc"
    filip = "F2000_"+casename(i)+".PS.mon.0601-3012.nc"

    ft = addfile(diri+filit, "r")
    fp = addfile(diri+filip, "r")
    T = ft->T ;(:,:,{latstr:latend},{lonstr:lonend})
    PS = fp->PS ;(:,:,{latstr:latend},{lonstr:lonend})
    ; printVarSummary(u1)
    plev = T&lev
    plev = plev*100
    dp = dpres_plevel_Wrap(plev, PS, min(plev), 0)
    Tdp = T*dp
    copy_VarCoords(T, Tdp)
    t = dim_sum_n_Wrap(Tdp(:,{levs:levt},:,:), 1)/dim_sum_n_Wrap(dp(:,{levs:levt},:,:), 1)
    copy_VarCoords(T(:,0,:,:), t)


; end do
; u!0 = "case"
; u&case = casename
; copy_VarCoords(u, v)
; printVarSummary(u)




tpmaskpath = "/home/yangsong3/wuyt/sysu/scripts/TP_CESM_F/model_need_input_data/"
; tpmaskname = "TPmask_1000m_C20190410.nc"
tpmaskname = "TPmask_1500m_C20200713.nc"
ftp = addfile(tpmaskpath+tpmaskname, "r")

tpmask = ftp->tp
printVarSummary(tpmask)

t =  mask(t, tpmask.ge.0.5, True)

t_TP = wgt_areaave_Wrap(t(:,{25:40},{75:105}), 1, 1, 0)
t_mar = t_TP(2::12)
printVarSummary(t_mar)

index = dim_standardize_Wrap(t_mar, opt)
end