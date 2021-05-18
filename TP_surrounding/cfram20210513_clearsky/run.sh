
#varname=("albedo" "warm" "co2" "drdt" "o3" "solar" "ta" "wv")
varname=("albedo" "co2" "ta" "wv" "o3" "solar" "warm")
size=${#varname[*]}
for ((i=0; i<${size}; i++))
do
    nohup ./1-${varname[i]}.out > 1-${varname[i]}.txt &
done





