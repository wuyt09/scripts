
varname=("warm" "albedo" "cloud" "co2" "drdt" "o3" "solar" "t" "wv")
size=${#varname[*]}
for ((i=0; i<${size}; i++))
do
	#ifort ${varname[i]}.f -o ${varname[i]}.out
    nohup ./${varname[i]}.out > ${varname[i]}.txt &
done





