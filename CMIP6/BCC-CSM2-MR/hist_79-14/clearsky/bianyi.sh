
varname=("albedo" "co2" "drdt" "o3" "solar" "t" "wv" "warm")
#varname=("cloud")
size=${#varname[*]}
for ((i=0; i<${size}; i++))
do
	ifort 1-${varname[i]}.f -o 1-${varname[i]}.out
done





