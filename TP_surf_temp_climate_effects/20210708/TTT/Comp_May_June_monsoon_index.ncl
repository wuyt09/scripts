begin

yrStrt = 1979
yrLast = 2020
year = ispan(yrStrt, yrLast, 1)
nyear = dimsizes(year)

monsoon_name = (/"WY","SA-Goswami","SEA-Wang","EA-LKY"/)

fili1 = "../May_monsoon_index_1979_2020_ERA5.txt"
fili2 = "../June_monsoon_index_1979_2020_ERA5.txt"

x1 = asciiread(fili1, (/nyear,4/), "float")
x2 = asciiread(fili2, (/nyear,4/), "float")
print(x1(nyear-1,:))
print(x2(nyear-1,:))



end