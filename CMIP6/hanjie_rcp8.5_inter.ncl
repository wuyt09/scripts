load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
; ==============================================================
begin
; ==============================================================
; dir & in
; ==============================================================
  diri = "/homes/hfan7/data/bcc-csm1-1/bcc-csm1-1_Amon_rcp85_r1i1p1/"   ; input directory
  fils_cl = systemfunc ("ls "+diri+"cl_Amon_bcc-csm1-1_rcp85_r1i1p1_*")    ; file paths
  fin_cl = addfiles(fils_cl, "r")

  fils_cli = systemfunc ("ls "+diri+"cli_Amon_bcc-csm1-1_rcp85_r1i1p1_*")    ; file paths
  fin_cli = addfiles(fils_cli, "r")

  fils_clw = systemfunc ("ls "+diri+"clw_Amon_bcc-csm1-1_rcp85_r1i1p1_*")    ; file paths
  fin_clw = addfiles(fils_clw, "r")

  fils_hus = systemfunc ("ls "+diri+"hus_Amon_bcc-csm1-1_rcp85_r1i1p1_*")    ; file paths
  fin_hus = addfiles(fils_hus, "r")

  fils_ta = systemfunc ("ls "+diri+"ta_Amon_bcc-csm1-1_rcp85_r1i1p1_*")    ; file paths
  fin_ta = addfiles(fils_ta, "r")

  ; fils_tro3 = systemfunc ("ls "+diri+"tro3_Amon_bcc-csm1-1_rcp85_r1i1p1_*")    ; file paths
  ; fin_tro3 = addfiles(fils_tro3, "r")

  ; fils_zg = systemfunc ("ls "+diri+"zg_Amon_bcc-csm1-1_rcp85_r1i1p1_*")    ; file paths
  ; fin_zg = addfiles(fils_zg, "r")

  fils_hfls = systemfunc ("ls "+diri+"hfls_Amon_bcc-csm1-1_rcp85_r1i1p1_*")    ; file paths
  fin_hfls = addfiles(fils_hfls, "r")

  fils_hfss = systemfunc ("ls "+diri+"hfss_Amon_bcc-csm1-1_rcp85_r1i1p1_*")    ; file paths
  fin_hfss = addfiles(fils_hfss, "r")

  fils_huss = systemfunc ("ls "+diri+"huss_Amon_bcc-csm1-1_rcp85_r1i1p1_*")    ; file paths
  fin_huss = addfiles(fils_huss, "r")

  fils_ps = systemfunc ("ls "+diri+"ps_Amon_bcc-csm1-1_rcp85_r1i1p1_*")    ; file paths
  fin_ps = addfiles(fils_ps, "r")

  fils_rsds = systemfunc ("ls "+diri+"rsds_Amon_bcc-csm1-1_rcp85_r1i1p1_*")    ; file paths
  fin_rsds = addfiles(fils_rsds, "r")

  fils_rsdt = systemfunc ("ls "+diri+"rsdt_Amon_bcc-csm1-1_rcp85_r1i1p1_*")    ; file paths
  fin_rsdt = addfiles(fils_rsdt, "r")

  fils_rsus = systemfunc ("ls "+diri+"rsus_Amon_bcc-csm1-1_rcp85_r1i1p1_*")    ; file paths
  fin_rsus = addfiles(fils_rsus, "r")

  fils_ts = systemfunc ("ls "+diri+"ts_Amon_bcc-csm1-1_rcp85_r1i1p1_*")    ; file paths
  fin_ts = addfiles(fils_ts, "r")

; ==============================================================
; basic setting 
; ==============================================================
  time   = fin_cl[:]->time 
  XXXX     = cd_calendar(time,-1 )
  t_index_start = ind(XXXX.eq.205101)
  t_index_end   = ind(XXXX.eq.210012)

  plev = (/10.,20.,30.,50.,70.,100.,150.,\
          200.,250.,300.,400.,500.,\
          600.,700.,850.,925.,1000./)

  plev!0        = "plev"                  ; variable/dim name 
  plev&plev     =  plev                   ; create coordinate variable
  plev@long_name = "pressure"               ; attach some attributes
  plev@units     = "hPa"

  pin = fin_ta[0]->plev(::-1)

  a = fin_cl[0]->a(::-1)
  b = fin_cl[0]->b(::-1)
  psfc = fin_cl[:]->ps(t_index_start:t_index_end,:,:)
  psfc_annual = month_to_annual(psfc, 1)
  psfc_annual_comp = dim_avg_n_Wrap(psfc_annual, 0)
  p0 = 0.01*fin_cl[0]->p0

; ==============================================================
; read data
; ==============================================================
  hfls_org = fin_hfls[:]->hfls(t_index_start:t_index_end,:,:)
  hfls = -hfls_org
  copy_VarAtts(hfls_org, hfls)
  copy_VarCoords(hfls_org, hfls)
  copy_VarMeta(hfls_org, hfls)
  printMinMax(hfls, False)

  hfss_org = fin_hfss[:]->hfss(t_index_start:t_index_end,:,:)
  hfss = -hfss_org
  copy_VarAtts(hfss_org, hfss)
  copy_VarCoords(hfss_org, hfss)
  copy_VarMeta(hfss_org, hfss)
  printMinMax(hfss, False)

  huss = fin_huss[:]->huss(t_index_start:t_index_end,:,:)
  printMinMax(huss, False)

  ps = fin_ps[:]->ps(t_index_start:t_index_end,:,:)
  printMinMax(ps, False)

  rsds = fin_rsds[:]->rsds(t_index_start:t_index_end,:,:)
  printMinMax(rsds, False)

  rsdt = fin_rsdt[:]->rsdt(t_index_start:t_index_end,:,:)
  printMinMax(rsdt, False)

  rsus = fin_rsus[:]->rsus(t_index_start:t_index_end,:,:)
  printMinMax(rsus, False)

  ts = fin_ts[:]->ts(t_index_start:t_index_end,:,:)
  printMinMax(ts, False)

  cl_org = fin_cl[:]->cl(t_index_start:t_index_end,::-1,:,:)
  cl = cl_org * 0.01 
  copy_VarAtts(cl_org, cl)
  copy_VarCoords(cl_org, cl)
  copy_VarMeta(cl_org, cl)
  printMinMax(cl, False)

  cli = fin_cli[:]->cli(t_index_start:t_index_end,::-1,:,:)
  printMinMax(cli, False)

  clw = fin_clw[:]->clw(t_index_start:t_index_end,::-1,:,:)
  printMinMax(clw, False)

  hus = fin_hus[:]->hus(t_index_start:t_index_end,::-1,:,:)
  printMinMax(hus, False)

  ta = fin_ta[:]->ta(t_index_start:t_index_end,::-1,:,:)
  printMinMax(ta, False)

  ; tro3_org = fin_tro3[:]->tro3(t_index_start:t_index_end,::-1,:,:)
  ; tro3 = tro3_org*1e-9
  ; copy_VarAtts(tro3_org, tro3)
  ; copy_VarCoords(tro3_org, tro3)
  ; copy_VarMeta(tro3_org, tro3)
  ; printMinMax(tro3, False)

  ; zg = fin_zg[:]->zg(t_index_start:t_index_end,::-1,:,:)
  ; printMinMax(zg, False)

; ==============================================================
; annual composite 
; ==============================================================
  hfls_annual = month_to_annual(hfls, 1)
  hfls_annual_comp = dim_avg_n_Wrap(hfls_annual, 0)
  printVarSummary(hfls_annual_comp)
  printMinMax(hfls_annual_comp, False)

; ==============================================================
  hfss_annual = month_to_annual(hfss, 1)
  hfss_annual_comp = dim_avg_n_Wrap(hfss_annual, 0)
  printVarSummary(hfss_annual_comp)
  printMinMax(hfss_annual_comp, False)

;==============================================================
  huss_annual = month_to_annual(huss, 1)
  huss_annual_comp = dim_avg_n_Wrap(huss_annual, 0)
  printVarSummary(huss_annual_comp)
  printMinMax(huss_annual_comp, False)

; ==============================================================
  ps_annual = month_to_annual(ps, 1)
  ps_annual_comp = dim_avg_n_Wrap(ps_annual, 0)
  printVarSummary(ps_annual_comp)
  printMinMax(ps_annual_comp, False)

; ==============================================================
  rsds_annual = month_to_annual(rsds, 1)
  rsds_annual_comp = dim_avg_n_Wrap(rsds_annual, 0)
  printVarSummary(rsds_annual_comp)
  printMinMax(rsds_annual_comp, False)

; ==============================================================
  rsdt_annual = month_to_annual(rsdt, 1)
  rsdt_annual_comp = dim_avg_n_Wrap(rsdt_annual, 0)
  printVarSummary(rsdt_annual_comp)
  printMinMax(rsdt_annual_comp, False)

; ==============================================================
  rsus_annual = month_to_annual(rsus, 1)
  rsus_annual_comp = dim_avg_n_Wrap(rsus_annual, 0)
  printVarSummary(rsus_annual_comp)
  printMinMax(rsus_annual_comp, False)

; ==============================================================
  ts_annual = month_to_annual(ts, 1)
  ts_annual_comp = dim_avg_n_Wrap(ts_annual, 0)
  printVarSummary(ts_annual_comp)
  printMinMax(ts_annual_comp, False)

; ============================================================== 
  cl_annual = month_to_annual(cl, 1)
  cl_annual_comp = dim_avg_n_Wrap(cl_annual, 0)
  printVarSummary(cl_annual_comp)
  printMinMax(cl_annual_comp, False)

; ==============================================================
  cli_annual = month_to_annual(cli, 1)
  cli_annual_comp = dim_avg_n_Wrap(cli_annual, 0)
  printVarSummary(cli_annual_comp)
  printMinMax(cli_annual_comp, False)

; ==============================================================
  clw_annual = month_to_annual(clw, 1)
  clw_annual_comp = dim_avg_n_Wrap(clw_annual, 0)
  printVarSummary(clw_annual_comp)
  printMinMax(clw_annual_comp, False)

; ==============================================================
  hus_annual = month_to_annual(hus, 1)
  hus_annual_comp = dim_avg_n_Wrap(hus_annual, 0)
  printVarSummary(hus_annual_comp)
  printMinMax(hus_annual_comp, False)

; ==============================================================
  ta_annual = month_to_annual(ta, 1)
  ta_annual_comp = dim_avg_n_Wrap(ta_annual, 0)
  printVarSummary(ta_annual_comp)
  printMinMax(ta_annual_comp, False)

; ==============================================================
  ; tro3_annual = month_to_annual(tro3, 1)
  ; tro3_annual_comp = dim_avg_n_Wrap(tro3_annual, 0)
  ; printVarSummary(tro3_annual_comp)
  ; printMinMax(tro3_annual_comp, False)

; ==============================================================
; interpolation to p grid
; ==============================================================
  phis = new(dimsizes(ta_annual_comp({1000},:,:)), "float")
  tbot = ta_annual_comp({1000},:,:)

  cl_p = vinth2p_ecmwf(cl_annual_comp, a, b, plev, psfc_annual_comp, 2, p0, 1, True, 0, tbot, phis)
  printVarSummary(cl_p)
  printMinMax(cl_p, False)

  cli_p = vinth2p_ecmwf(cli_annual_comp, a, b, plev, psfc_annual_comp, 2, p0, 1, True, 0, tbot, phis)
  printVarSummary(cli_p)
  printMinMax(cli_p, False)

  clw_p = vinth2p_ecmwf(clw_annual_comp, a, b, plev, psfc_annual_comp, 2, p0, 1, True, 0, tbot, phis)
  printVarSummary(clw_p)
  printMinMax(clw_p, False)

  hus_p = int2p_n_Wrap(pin, hus_annual_comp, pin, -2, 0)
  printVarSummary(hus_p)
  printMinMax(hus_p, False)

  ta_p = int2p_n_Wrap(pin, ta_annual_comp, pin, -2, 0)
  printVarSummary(ta_p)
  printMinMax(ta_p, False)
  
  ; tro3_p = int2p_n_Wrap(pin, tro3_annual_comp, pin, -2, 0)
  ; printVarSummary(tro3_p)
  ; printMinMax(tro3_p, False)

; ==============================================================
; output data as dat
; ==============================================================
  fbindirwrite("slhf_warm.dat", hfls_annual_comp)

  fbindirwrite("sshf_warm.dat", hfss_annual_comp)

  fbindirwrite("huss_warm.dat", huss_annual_comp)

  fbindirwrite("sp_warm.dat", ps_annual_comp)

  fbindirwrite("ssrd_warm.dat", rsds_annual_comp)

  fbindirwrite("solarin_warm.dat", rsdt_annual_comp)

  fbindirwrite("ssru_warm.dat", rsus_annual_comp)

  fbindirwrite("skt_warm.dat", ts_annual_comp)

  fbindirwrite("cc_warm.dat", cl_p)

  fbindirwrite("ciwc_warm.dat", cli_p)

  fbindirwrite("clwc_warm.dat", clw_p)

  fbindirwrite("hus_warm.dat", hus_p)

  fbindirwrite("t_warm.dat", ta_p)

  ; fbindirwrite("o3_warm.dat", tro3_p)

end