
begin

;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20191116/"
figname = "monsoon_index_area"

lonstr = 30
lonend = 160
latstr = -10
latend = 60

wks = gsn_open_wks("pdf", figpath+figname)
plot = new(1, graphic)


res = True
res@gsnDraw = False
res@gsnFrame = False
res@gsnMaximize = True

res@mpOutlineOn = True
res@mpFillOn = True
res@mpGeophysicalLineThicknessF = 2
res@mpGeophysicalLineColor = "grey45" ;"grey25"

res@mpMinLatF = latstr
res@mpMaxLatF = latend
res@mpMinLonF = lonstr
res@mpMaxLonF = lonend

res@gsnStringFontHeightF = 0.03
res@tmXBLabelFontHeightF = 0.02
res@tmYLLabelFontHeightF = 0.02


plot(0) = gsn_csm_map_ce(wks, res)


gres = True
gres@gsLineThicknessF = 3.0
;--------------------------------------------------------------
latWY1 = 0
latWY2 = 20
lonWY1 = 40
lonWY2 = 110
gres@gsLineColor = "red"
glonWY = (/lonWY2,lonWY2,lonWY1,lonWY1,lonWY2/)
glatWY = (/latWY1,latWY2,latWY2,latWY1,latWY1/)
dumWY = gsn_add_polyline(wks, plot(0), glonWY, glatWY, gres)
;--------------------------------------------------------------
;--------------------------------------------------------------
latSEA11 = 5
latSEA12 = 15
lonSEA11 = 90
lonSEA12 = 130

latSEA21 = 22.5
latSEA22 = 32.5
lonSEA21 = 110
lonSEA22 = 140
gres@gsLineColor = "blue"
glonSEA1 = (/lonSEA12,lonSEA12,lonSEA11,lonSEA11,lonSEA12/)
glatSEA1 = (/latSEA11,latSEA12,latSEA12,latSEA11,latSEA11/)

glonSEA2 = (/lonSEA22,lonSEA22,lonSEA21,lonSEA21,lonSEA22/)
glatSEA2 = (/latSEA21,latSEA22,latSEA22,latSEA21,latSEA21/)
dumSEA1 = gsn_add_polyline(wks, plot(0), glonSEA1, glatSEA1, gres)
dumSEA2 = gsn_add_polyline(wks, plot(0), glonSEA2, glatSEA2, gres)
;--------------------------------------------------------------
;--------------------------------------------------------------
latSA1 = 10
latSA2 = 30
lonSA1 = 70
lonSA2 = 110
gres@gsLineColor = "green4"
glonSA = (/lonSA2,lonSA2,lonSA1,lonSA1,lonSA2/)
glatSA = (/latSA1,latSA2,latSA2,latSA1,latSA1/)
dumSA = gsn_add_polyline(wks, plot(0), glonSA, glatSA, gres)
;--------------------------------------------------------------
;--------------------------------------------------------------
latEA11 = 25
latEA12 = 35
lonEA11 = 110
lonEA12 = 150

latEA21 = 40
latEA22 = 50
lonEA21 = 110
lonEA22 = 150
gres@gsLineColor = "black"
glonEA1 = (/lonEA12,lonEA12,lonEA11,lonEA11,lonEA12/)
glatEA1 = (/latEA11,latEA12,latEA12,latEA11,latEA11/)

glonEA2 = (/lonEA22,lonEA22,lonEA21,lonEA21,lonEA22/)
glatEA2 = (/latEA21,latEA22,latEA22,latEA21,latEA21/)
dumEA1 = gsn_add_polyline(wks, plot(0), glonEA1, glatEA1, gres)
dumEA2 = gsn_add_polyline(wks, plot(0), glonEA2, glatEA2, gres)
;--------------------------------------------------------------

draw(plot)

end