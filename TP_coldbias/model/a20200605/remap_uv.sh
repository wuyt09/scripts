module load cdo-bin
casename=("Ctrl_FAMIP" "Hist_FAMIP")
size=${#casename[*]}
for ((i=0; i<${size}; i++))
do
	cdo
	#cdo -remapnn,r360x181 output/CanESM2_r1_r5_AR1_u850hist.nc output/CanESM2_remap_r1_r5_AR1_u850hist.nc
	#cdo -remapnn,r360x181 output/CanESM2_r1_r5_AR2_u850hist.nc output/CanESM2_remap_r1_r5_AR2_u850hist.nc
	cdo -remapbil,r360x181 /home/yangsong3/data-model/wuyt/TPbias/TPbias_${casename[i]}/a20191206/${casename[i]}.cam.h0.T.1979-2005.nc ./${casename[i]}.cam.h0.T.1979-2005.1x1.nc
	cdo -remapbil,r360x181 /home/yangsong3/data-model/wuyt/TPbias/TPbias_${casename[i]}/a20191206/${casename[i]}.cam.h0.PS.1979-2005.nc ./${casename[i]}.cam.h0.PS.1979-2005.1x1.nc
	echo "${casename[i]} is done!"
done
