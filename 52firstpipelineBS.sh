#!/bin/bash
# To run simply write ./pipelineBS.sh "sourcename" "night" dispCut "period"<-- So, the first two arguments (and the last one) have to be in quotes (strings in quotes, numbers without quotes)

Source=$1
Night=$2
DispCut=$3
BinWidth=$4
Period=$5

PathIn="/data/magic/users-ifae/edosouto/MAGIC/BigSearch/Melibea/${Period}"
PathOut="/data/magic/users-ifae/edosouto/MAGIC/BigSearch/Code/MARS3version/${Period}"
CodePath="/data/magic/users-ifae/edosouto/MAGIC/BigSearch/Code/MARS3version"

# Security check in case we're recreating a LightCurves file
rm -rf ${PathOut}/05to35/LightCurves_${Night}_${Source}_${DispCut}.root
rm -rf ${PathOut}/35to50/LightCurves_${Night}_${Source}_${DispCut}.root

#Create a counter to fill the pid array correctly
declare -i ct=0

echo "${Source} ${Night} with DispCut ${DispCut}deg and bin width ${BinWidth}s."
echo "Running EventExtractorExcl"

#run processes and store pids in array
#for File in $PathIn/05to35/$Source/$Night/20*.root
##TODO add the next zenith range
#do
#    echo $File
#    root -l -b -q $CodePath/EventExtractorLowZd.cpp+\(\"$File\",$DispCut\) &
#    pids[${ct}]=$!
#    ct+=1
#done
#
## wait for all pids
#for pid in ${pids[*]}
#do
#    wait $pid
#done
## The intermediate files EL_night_source_run_dispcut.root have been created

echo "Number of runs processed = ${#pids[@]}"
echo "Running CellCreator for the low Zd range"
echo "${PathOut}/05to35/EL_${Night}_${Source}*${DispCut}.root"

if [ `ls -1 ${PathOut}/05to35/EL_${Night}_${Source}*${DispCut}.root 2>/dev/null | wc -l ` -gt 0 ];
then
    echo "ok"
    for f in ${PathOut}/05to35/EL_${Night}_${Source}*${DispCut}.root
    do
        echo ${f}
        root -l -b -q ${CodePath}/CellCreator.cpp+\(\"${f}\",\"${Night}\",\"${Source}\",${DispCut}\)
    done
else
    echo "No EL files in the low range. Exiting."
fi

