
varname=("albedo" "cloud"  "ts" "o3" "ta" "wv" "warm")
#varname=("albedo" "co2" "drdt" "o3" "solar")
size=${#varname[*]}
for ((i=0; i<${size}; i++))
do
    nohup ./1-${varname[i]}.out > 1-${varname[i]}.txt &
done





