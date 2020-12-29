
begin

x = fspan(1, 36, 36)
x!0 = "time"
x&time = ispan(1, 36, 1)
x_djf = month_to_season(x, "DJF")
x_jja = month_to_season(x, "JJA")

print(x_jja)
print(x_djf)

end