
begin

yrStrt = 1981
yrLast = 2005

f = addfile("./ghg_hist_1765-2005_c091218.nc", "r")

TIME = f->time
YYYY = cd_calendar(TIME, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

CO2 = f->CO2(iYYYY)
CH4 = f->CH4(iYYYY)
N2O = f->N2O(iYYYY)

co2 = avg(CO2)
ch4 = avg(CH4)
n2o = avg(N2O)
print("CO2 "+co2)
print("CH4 "+ch4)
print("N2O "+n2o)

end 