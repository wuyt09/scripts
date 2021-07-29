
#varname=("albedo" "cloud" "co2" "drdt" "o3" "solar" "ta" "ts" "wv")
varname=("albedo" "cloud" "warm" "ts" "o3" "ta" "wv")
#varname=("cloud")
size=${#varname[*]}
for ((i=0; i<${size}; i++))
do
	ifort 1-${varname[i]}.f -o 1-${varname[i]}.out
done





