
begin

lonstr = 70
lonend = 170
latstr = -5
latend = 60

casename = (/"CTRL","NOQ1","NOTP_NOQ1","NOTP"/)

do i = 0, dimsizes(casename)-1
diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(i)+"_F2000_F19G16_CAM4/pro/"
filiu = "F2000_"+casename(i)+".plevel.U.mon.0601-3012.nc"
filiv = "F2000_"+casename(i)+".plevel.V.mon.0601-3012.nc"
filiz = "F2000_"+casename(i)+".plevel.Z3.mon.0601-3012.nc"

fu = addfile(diri+filiu, "r")
fv = addfile(diri+filiv, "r")
fz = addfile(diri+filiz, "r")
u1 = fu->U(:,:,{latstr:latend},{lonstr:lonend})
v1 = fv->V(:,:,{latstr:latend},{lonstr:lonend})
z1 = fz->Z3(:,:,{latstr:latend},{lonstr:lonend})
printVarSummary(u1)

u = clmMonTLLL(u1)
v = clmMonTLLL(v1)
z = clmMonTLLL(z1)
printVarSummary(u)

;===============================================================================
; Chen et al. (2000)
; averaged V 10m over (10-25N 110-130E) and (25-40N, 120-140E) 
;===============================================================================
v_chen = v(:,36,:,:) ; 1000hPa
index_chen1 = wgt_areaave_Wrap(v_chen(:,{10:25},{110:130}), 1, 1, 0)
index_chen2 = wgt_areaave_Wrap(v_chen(:,{25:40},{120:140}), 1, 1, 0)
index_chen = (index_chen1 + index_chen2)/2
index_chen = (-1.) * index_chen
copy_VarCoords(index_chen1, index_chen)
printVarSummary(index_chen)

;===============================================================================
; Yang et al. (2002)
; V850 (20N-40N, 100E-140E))
;===============================================================================
u_yang = v(:,{850},:,:)
index_yang = wgt_areaave_Wrap(u_yang(:,{20:40},{100:140}), 1, 1, 0)
index_yang = (-1.) * index_yang
copy_VarCoords(u_yang(:,0,0), index_yang)
printVarSummary(index_yang)

;===============================================================================
; Li and Yang (2010)
; {[U200(30N-35N, 90E-160E) - U200(50N-60N, 70E-170E)] + 
; [U200(30N-35N, 90E-160E) - U200(5S-10N, 90E-160E)]}/2
;===============================================================================
u_lyang = u(:,{200},:,:)
u_lyang1 = wgt_areaave_Wrap(u_lyang(:,{30:35},{90:160}), 1, 1, 0)
u_lyang2 = wgt_areaave_Wrap(u_lyang(:,{50:60},{70:170}), 1, 1, 0)
u_lyang3 = wgt_areaave_Wrap(u_lyang(:,{-5:10},{90:160}), 1, 1, 0)
index_lyang = (u_lyang1 - u_lyang2 + u_lyang1 - u_lyang3)/2
copy_VarCoords(u_lyang(:,0,0), index_lyang)
printVarSummary(index_lyang)

;===============================================================================
; Jhun and Lee (2004)
; U300 [(27.5N-37.5N, 110E-170E) - (50N-60N, 80E-140E)]
;===============================================================================
u_jhun = u(:,{300},:,:)
u1_jhun = wgt_areaave_Wrap(u_jhun(:,{27.5:37.5},{110:170}), 1, 1, 0)
u2_jhun = wgt_areaave_Wrap(u_jhun(:,{50:60},{80:140}), 1, 1, 0)
index_jhun = u1_jhun - u2_jhun
copy_VarCoords(u1_jhun, index_jhun)
printVarSummary(index_jhun)

;===============================================================================
; Sun and Li (1997)
; Z500 (30-45N, 125-145E)
;===============================================================================
z_sun = z(:,{500},:,:)
index_sun = wgt_areaave_Wrap(z_sun(:,{30:45},{125:145}), 1, 1, 0)
index_sun = (-1.) * index_sun
copy_VarCoords(z_sun(:,0,0), index_sun)
printVarSummary(index_sun)

;===============================================================================
; wtite data to txt file
;===============================================================================
filipath = "./"
filiname = "calc_monsoon_index_12month_"+casename(i)+".txt"

write_table(filipath+filiname, "w", [/"Chen","Yang","LYang","Jhun","Sun"/], "%10s %10s %10s %10s %10s")
write_table(filipath+filiname, "a", [/index_chen,index_yang,index_lyang,index_jhun,index_sun/], "%10.2f %10.2f %10.2f %10.2f %10.2f")
;print("end of "+casename(i))
end do

end
