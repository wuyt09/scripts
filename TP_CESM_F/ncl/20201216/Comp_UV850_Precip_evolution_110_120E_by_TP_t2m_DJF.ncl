
begin

yrStrt = 1979
yrLast = 2018
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


fx = "time_series_TP_t2m_DJF_MAM_OBS.txt"
x = asciiread(fx, (/nyear,3/), "float")
inx = x(0:nyear-2,1)

hig_ind = ind(inx.ge.1.)
low_ind = ind(inx.le.-1.)
print(year(hig_ind))
print(year(low_ind))

hig_ind = hig_ind+1
low_ind = low_ind+1

;-------------------------------------------------------------------
;-------------------------------------------------------------------
lonW = 110
lonE = 120
LEV = 850

U = new((/365,nlat,nlon/), float)
V = U
PREC = U

u = new((/73,nlat,nlon/), float)
v = u
prec = u
u_band = new((/nyear,73,nlat/), float)
v_band = u_band
prec_band = u_band

diriu = "/home/yangsong3/data-observation/ERA5-daily/Multi-levels/uwind/"
diriv = "/home/yangsong3/data-observation/ERA5-daily/Multi-levels/vwind/"
dirip = "/home/yangsong3/data-observation/ERA5-daily/Single-level/prect/"

do i = 0, nyear-1
  irecS = 0
  do j = 0, nmonth-1
      mth = sprinti("%0.2i", month(j))
      filiu = "uwind."+year(i)+"-"+mth+".daily.nc"
      filiv = "vwind."+year(i)+"-"+mth+".daily.nc"
      filip = "prect."+year(i)+"-"+mth+".daily.nc"

      fu = addfile(diriu+filiu, "r")    
      UX := short2flt(fu->u(:,{LEV},::-1,:))
      fv = addfile(diriv+filiv, "r")    
      VX := short2flt(fv->v(:,{LEV},::-1,:))
      fp = addfile(dirip+filip, "r")    
      PX := short2flt(fp->tp(:,::-1,:))
      PX = PX*24*1000 
      PX@units = "mm/day"

      ndays := dimsizes(UX&time)
      if (ndays.eq.29) then
          irecL = irecS+ndays-2
          U(irecS:irecL,:,:) = UX(0:27,:,:)
          V(irecS:irecL,:,:) = VX(0:27,:,:)
          PREC(irecS:irecL,:,:) = PX(0:27,:,:)
          irecS = irecS+ndays-1
      else
          irecL = irecS+ndays-1
          U(irecS:irecL,:,:) = UX 
          V(irecS:irecL,:,:) = VX 
          PREC(irecS:irecL,:,:) = PX 
          irecS = irecS+ndays
      end if
      ; print(irecL)
  end do
  ; print(irecL)
  do p = 0, 72
      u(p,:,:) = dim_sum_n_Wrap(U(5*p:5*p+4,:,:), 0)/5.0
      v(p,:,:) = dim_sum_n_Wrap(V(5*p:5*p+4,:,:), 0)/5.0
      prec(p,:,:) = dim_sum_n_Wrap(PREC(5*p:5*p+4,:,:), 0)/5.0
  end do
  copy_VarCoords(U(0,:,:), u(0,:,:))
  copy_VarCoords(V(0,:,:), v(0,:,:))
  copy_VarCoords(PREC(0,:,:), prec(0,:,:))
  u_band(i,:,:) = dim_avg_n_Wrap(u(:,:,{lonW:lonE}), 2)
  v_band(i,:,:) = dim_avg_n_Wrap(v(:,:,{lonW:lonE}), 2)
  prec_band(i,:,:) = dim_avg_n_Wrap(prec(:,:,{lonW:lonE}), 2)
end do

u_band!0 = "year"
u_band&year = year
u_band!1 = "pentad"
u_band&pentad = penta
copy_VarCoords(u_band, v_band)
copy_VarCoords(u_band, prec_band)

u_hig = dim_avg_n_Wrap(u_band(hig_ind,:,:), 0)
u_low = dim_avg_n_Wrap(u_band(low_ind,:,:), 0)
v_hig = dim_avg_n_Wrap(v_band(hig_ind,:,:), 0)
v_low = dim_avg_n_Wrap(v_band(low_ind,:,:), 0)
prec_hig = dim_avg_n_Wrap(prec_band(hig_ind,:,:), 0)
prec_low = dim_avg_n_Wrap(prec_band(low_ind,:,:), 0)

ud_hig = u_hig(latitude|:,pentad|:)
vd_hig = v_hig(latitude|:,pentad|:)
precd_hig = prec_hig(latitude|:,pentad|:)
ud_low = u_low(latitude|:,pentad|:)
vd_low = v_low(latitude|:,pentad|:)
precd_low = prec_low(latitude|:,pentad|:)
printVarSummary(ud_hig)

;-------------------------------------------------------------------
;-------------------------------------------------------------------

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20201216/"
figname = "Comp_UV850_Precip_evolution_110_120E_by_TP_t2m_DJF"

wks = gsn_open_wks("pdf", figpath+figname)

pens = 19
pene = 65
latstr = 10
latend = 50

nplot = 2
plot = new(nplot, graphic)
plotd = plot


res          = True                ; plot mods desired
res@gsnDraw  = False
res@gsnFrame = False
res@vpHeightF = 0.6
res@vpWidthF  = 0.8


; res@gsnScalarContour = True
res@cnFillOn         = True                ; color on
res@cnFillPalette    = "precip2_17lev"    ; set color map
res@lbLabelBarOn = False
res@cnInfoLabelOn = False
res@cnLinesOn  = False
res@cnLineLabelsOn = False

res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = 2
res@cnMaxLevelValF = 20
res@cnLevelSpacingF  = 1                  ; contour spacing


res@tmYLMinorOn = False
res@tmYLMode = "Explicit"
res@tmYLValues = ispan(latstr, latend, 5)
res@tmYLLabels = (/"10N","15N","20N","25N","30N","35N","40N","45N","50N"/)
res@tmYLLabelFontHeightF = 0.02


res@tmXBPrecision = 2
res@tiXAxisString    = "Pentad"              ; y axis title
res@gsnStringFontHeightF = 0.03

vcres = True
vcres@gsnDraw  = False
vcres@gsnFrame = False
vcres@vcRefLengthF            = 0.04            ;the length of the reference magnitude
vcres@vcRefMagnitudeF         = 5.0
vcres@vcMinDistanceF          = 0.018            ;the minimum distance of neighboring vectors
vcres@vcGlyphStyle            = "CurlyVector"   ;the style of vector
vcres@vcLineArrowThicknessF   = 3.0

vcres@vcRefAnnoOrthogonalPosF = -0.4          ;move ref vector down

vcres@vcRefAnnoString2On = False
vcres@vcRefAnnoFontHeightF = 0.02
vcres@vcLineArrowColor = "black" ;grey50
; res@vcVectorDrawOrder = "PostDraw"


  res@gsnLeftString = "High DJF TP t2m"
  plot(0) = gsn_csm_contour(wks, precd_hig({latstr:latend},pens-1:pene-1), res)
  plotd(0) = gsn_csm_vector(wks, ud_hig({latstr:latend},pens-1:pene-1), vd_hig({latstr:latend},pens-1:pene-1), vcres)
  overlay(plot(0), plotd(0))

  res@gsnLeftString = "Low DJF TP t2m"
  plot(1) = gsn_csm_contour(wks, precd_low({latstr:latend},pens-1:pene-1), res)
  plotd(1) = gsn_csm_vector(wks, ud_low({latstr:latend},pens-1:pene-1), vd_low({latstr:latend},pens-1:pene-1), vcres)
  overlay(plot(1), plotd(1))



pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = "UV"+LEV+" & Precip ("+lonW+"-"+lonE+"E)"

gsn_panel(wks, plot, (/1,2/), pres)


end