
varname=("albedo" "drdt"  "co2" "o3" "solar" "ta" "wv")
#varname=("albedo" "co2" "drdt" "o3" "solar")
size=${#varname[*]}
for ((i=0; i<${size}; i++))
do
    nohup ./1-${varname[i]}.out > 1-${varname[i]}.txt &
done





