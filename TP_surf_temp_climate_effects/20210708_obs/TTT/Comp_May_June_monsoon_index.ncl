begin

yrStrt = 1979
yrLast = 2020
year = ispan(yrStrt, yrLast, 1)
nyear = dimsizes(year)

monsoon_name = (/"WY","SA-Goswami","SEA-Wang","EA-LKY"/)

fili1 = "../May_monsoon_index_"+yrStrt+"-"+yrLast+"_ERA5.txt"
fili2 = "../June_monsoon_index_"+yrStrt+"-"+yrLast+"_ERA5.txt"

x1 = asciiread(fili1, (/nyear,4/), "float")
x2 = asciiread(fili2, (/nyear,4/), "float")
print(x1(nyear-1,:))
print(x2(nyear-1,:))

x1_clm = dim_avg_n_Wrap(x1, 0)
x2_clm = dim_avg_n_Wrap(x2, 0)
;-----------------------------------------------------------
indexfili = "./time_series_TP_TT_Mar_Apr_May_ERA5.txt"
inx = asciiread(indexfili, (/nyear,4/), "float")
inx!0 = "year"
inx&year = year

index = inx(:,2)
printVarSummary(index)

indhig := ind(index.ge.1.)
indlow := ind(index.le.-1.)
yearhig := year(indhig)
yearlow := year(indlow)
print(yearlow)
print(yearhig)
numhig = dimsizes(yearhig)
numlow = dimsizes(yearlow)

May_hig = dim_avg_n_Wrap(x1(indhig,:), 0)
May_low = dim_avg_n_Wrap(x1(indlow,:), 0)
June_hig = dim_avg_n_Wrap(x2(indhig,:), 0)
June_low = dim_avg_n_Wrap(x2(indlow,:), 0)

May_diff = May_hig - May_low
June_diff = June_hig - June_low

; print(monsoon_name+"  clm  "+sprintf("%4.2f", x1_clm))
; print(monsoon_name+"  clm  "+sprintf("%4.2f", x2_clm))
print(monsoon_name+"    "+sprintf("%4.2f", May_diff/x1_clm))
print(monsoon_name+"    "+sprintf("%4.2f", June_diff/x2_clm))

; print(monsoon_name+"  May High  "+sprintf("%4.2f", May_hig))
; print(monsoon_name+"  May Low  "+sprintf("%4.2f", May_low))
; print(monsoon_name+"  June High  "+sprintf("%4.2f", June_hig))
; print(monsoon_name+"  June Low  "+sprintf("%4.2f", June_low))

end