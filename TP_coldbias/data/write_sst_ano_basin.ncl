
begin
;=========================================================
; Note: the changes has been add as annual cycle, 
;       each year has the same anomalies.
;       but I still confused by how the anomaly data derived.
;                            2019-9-27 Wu Yuting
;=========================================================


; yrStrt = 1850
; yrLast = 2013
; year = ispan(yrStrt, yrLast, 1)

fili_ctrl = "sst_HadOIBl_bc_1.9x2.5_1850_2012_c130411.nc"
fili_hist = "Hist_bias_sst_cmip5_bc_1.9x2.5_c20191204.nc"

f_ctrl = addfile(fili_ctrl, "r") ; the land has value!!
f_hist = addfile(fili_hist, "r") ; the model doesnot read data over land!

sst1 = f_ctrl->SST_cpl
sst2 = f_hist->SST_cpl
printVarSummary(sst1)
printVarSummary(sst2)
lat = sst1&lat
lon = sst1&lon
; print(lat)
; print(lon)
; exit
hist_io = sst1
hist_po = sst1
hist_ao = sst1
hist_pio = sst1
hist_pao = sst1
hist_iao = sst1

;==============================================================
; Indian Ocean
;==============================================================
hist_io(:,{-30:30},{30:105}) = sst2(:,{-30:30},{30:105})
hist_io(:,{-30:-10},{105:130}) = sst2(:,{-30:-10},{105:130})

hist_pio(:,{-30:30},{30:105}) = sst2(:,{-30:30},{30:105})
hist_pio(:,{-30:-10},{105:130}) = sst2(:,{-30:-10},{105:130})

hist_iao(:,{-30:30},{30:105}) = sst2(:,{-30:30},{30:105})
hist_iao(:,{-30:-10},{105:130}) = sst2(:,{-30:-10},{105:130})

;==============================================================
; Pacific Ocean
;==============================================================
fili_landsea = "/home/yangsong3/soft/ncl650/lib/ncarg/data/cdf/landsea.nc"
flandsea = addfile(fili_landsea, "r")
LANDSEA = flandsea->LSMASK
landsea = linint2_Wrap(LANDSEA&lon, LANDSEA&lat, LANDSEA, True, sst1&lon, sst1&lat, 0)
printVarSummary(landsea)
;print(landsea({-30:30},{240:290}))
; x1 = hist_po({-30:30},{130:290})
; x2 = sst2({-30:30},{130:290})
; xlandsea = landsea({-30:30},{130:290})
; poindex = ind(ndtooned(xlandsea).eq.0)

; print(landsea({20},{240:270}))
; print(landsea({16},{240:270}))
; print(landsea({12},{240:270}))
; exit
hist_po(:,{-10:30},{105:130}) = sst2(:,{-10:30},{105:130})
hist_po(:,{-30:30},{130:240}) = sst2(:,{-30:30},{130:240})

hist_pio(:,{-10:30},{105:130}) = sst2(:,{-10:30},{105:130})
hist_pio(:,{-30:30},{130:240}) = sst2(:,{-30:30},{130:240})

hist_pao(:,{-10:30},{105:130}) = sst2(:,{-10:30},{105:130})
hist_pao(:,{-30:30},{130:240}) = sst2(:,{-30:30},{130:240})

do j = 32, 63 ; -29.36 to 29.36 56, 57;
    count = 0
    do i = 96, 116 ; 240 to 290
        if (landsea(j,i).lt.1 .and. count.eq.0) then
            hist_po(:,j,i) = sst2(:,j,i)
            hist_pio(:,j,i) = sst2(:,j,i)
            hist_pao(:,j,i) = sst2(:,j,i)
            print("lat "+lat(j)+"   lon "+(360-lon(i))) ;lon(i)
        else
            count = count + 1
        end if
    end do
end do
print("=====================================================")

;==============================================================
; Atlantic Ocean
;==============================================================
hist_ao(:,{-30:30},{0:30}) = sst2(:,{-30:30},{0:30})
hist_pao(:,{-30:30},{0:30}) = sst2(:,{-30:30},{0:30})
hist_iao(:,{-30:30},{0:30}) = sst2(:,{-30:30},{0:30})
do j = 32, 63 ; -29.36 to 29.36
    count = 0
    do i = 104, 143 ; 260 to 357.5
        if (landsea(j,i).eq.1) then
            count = count + 1
        end if
        if (landsea(j,i).eq.0 .and. count.ge.1) then ;
            hist_ao(:,j,i) = sst2(:,j,i)
            hist_pao(:,j,i) = sst2(:,j,i)
            hist_iao(:,j,i) = sst2(:,j,i)
            ;print("lat "+lat(j)+"   lon "+(360-lon(i))) ;lon(i)
        end if
    end do
end do
print("=====================================================")

SST1 = dim_avg_n_Wrap(sst1, 0)
SST2 = dim_avg_n_Wrap(sst2, 0)
io = dim_avg_n_Wrap(hist_io, 0)
po = dim_avg_n_Wrap(hist_po, 0)
ao = dim_avg_n_Wrap(hist_ao, 0)
pio = dim_avg_n_Wrap(hist_pio, 0)
pao = dim_avg_n_Wrap(hist_pao, 0)

;=======================================================
; write data
;=======================================================
date = systemfunc("date +%Y%m%d")
; Hist_IO = "IO_sst_cmip5_bc_1.9x2.5_c"+date+".nc"
; Hist_PO = "PO_sst_cmip5_bc_1.9x2.5_c"+date+".nc"
; Hist_AO = "AO_sst_cmip5_bc_1.9x2.5_c"+date+".nc"
Hist_PIO = "PIO_sst_cmip5_bc_1.9x2.5_c"+date+".nc"
Hist_PAO = "PAO_sst_cmip5_bc_1.9x2.5_c"+date+".nc"
; Hist_IAO = "IAO_sst_cmip5_bc_1.9x2.5_c"+date+".nc"
; system("rm "+Hist_IO)
; system("rm "+Hist_PO)
; system("rm "+Hist_AO)
system("rm "+Hist_PIO)
system("rm "+Hist_PAO)
; system("rm "+Hist_IAO)
; system("cp "+fili_ctrl+" "+Hist_IO)
; system("cp "+fili_ctrl+" "+Hist_PO)
; system("cp "+fili_ctrl+" "+Hist_AO)
system("cp "+fili_ctrl+" "+Hist_PIO)
system("cp "+fili_ctrl+" "+Hist_PAO)
; system("cp "+fili_ctrl+" "+Hist_IAO)
; fili_hist_io = addfile(Hist_IO, "w")
; fili_hist_po = addfile(Hist_PO, "w")
; fili_hist_ao = addfile(Hist_AO, "w")
fili_hist_pio = addfile(Hist_PIO, "w")
fili_hist_pao = addfile(Hist_PAO, "w")
; fili_hist_iao = addfile(Hist_IAO, "w")
; fili_hist_io->SST_cpl = hist_io
; fili_hist_po->SST_cpl = hist_po
; fili_hist_ao->SST_cpl = hist_ao
fili_hist_pio->SST_cpl = hist_pio
fili_hist_pao->SST_cpl = hist_pao
; fili_hist_iao->SST_cpl = hist_iao

;=======================================================
; plot 
;=======================================================
wks = gsn_open_wks("pdf", "write_sst_ano_basin")
plot = new(8, graphic)

res = True
res@gsnDraw  = False
res@gsnFrame = False
res@gsnMaximize = True

res@mpOutlineOn = True
res@mpFillOn = True
res@mpGeophysicalLineThicknessF = 2
res@mpGeophysicalLineColor = "grey25" ;"red" 
res@mpFillDrawOrder = "PostDraw"

res@mpCenterLonF = 180
res@mpMinLatF = -30
res@mpMaxLatF = 30

res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False

res@gsnStringFontHeightF = 0.015
res@tmXBLabelFontHeightF = 0.015
res@tmYLLabelFontHeightF = 0.015

res@lbLabelBarOn = False
res@pmLabelBarOrthogonalPosF = 0.25
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 0.5
res@cnMaxLevelValF = 30
res@cnMinLevelValF = 20
res@cnFillPalette = "CBR_coldhot"
res@gsnRightString = "K"
res@gsnStringFontHeightF = 0.025

res@gsnLeftString = "Ctrl"
plot(0) = gsn_csm_contour_map(wks, SST1, res)

res@gsnLeftString = "Hist"
plot(1) = gsn_csm_contour_map(wks, SST2, res)

res@gsnLeftString = "IndianO"
plot(2) = gsn_csm_contour_map(wks, io, res)

res@gsnLeftString = "PacificO"
plot(3) = gsn_csm_contour_map(wks, po, res)

res@gsnLeftString = "AtlanticO"
plot(4) = gsn_csm_contour_map(wks, ao, res)

res@gsnLeftString = "PacificO-InidanO"
plot(5) = gsn_csm_contour_map(wks, pio, res)

res@gsnLeftString = "PacificO-AtlanticO"
plot(6) = gsn_csm_contour_map(wks, pao, res)

res@gsnLeftString = "PacificO-AtlanticO"
plot(7) = gsn_csm_contour_map(wks, pio, res)

pres = True
pres@gsnMaximize = True
pres@gsnDraw = True
pres@gsnFrame = True
pres@gsnPanelLabelBar = True

gsn_panel(wks, plot, (/4,2/), pres)


end