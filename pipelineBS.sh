#!/bin/bash
# To run simply write ./pipelineBS.sh "sourcename" "night" dispCut "period"<-- So, the first two arguments (and the last one) have to be in quotes (strings in quotes, numbers without quotes)

Source=$1
Night=$2
DispCut=$3
BinWidth=$4
Period=$5
Create=$6

PathIn="/data/magic/users-ifae/edosouto/MAGIC/BigSearch/Melibea/${Period}"
PathOut="/data/magic/users-ifae/edosouto/MAGIC/BigSearch/Code/MARS3version/${Period}"
CodePath="/data/magic/users-ifae/edosouto/MAGIC/BigSearch/Code/MARS3version"

mkdir -p $PathOut/05to35/Results
mkdir -p $PathOut/35to50/Results

# I removed the acceptance maps part from here because it leads to error when sending many jobs

if [ "$Create" = true ] ; then

    #TODO fix this absurdity, change the fucking name
    # Security check in case we're recreating a LightCurves file
    rm -rf ${PathOut}/05to35/LightCurves_${Night}_${Source}_${DispCut}.root
    rm -rf ${PathOut}/35to50/LightCurves_${Night}_${Source}_${DispCut}.root
    
    #Create a counter to fill the pid array correctly
    declare -i ct=0
    
    echo "${Source} ${Night} with DispCut ${DispCut}deg and bin width ${BinWidth}s."
    echo "Running EventExtractorExcl"
    
    #run processes and store pids in array
    # The ST.03.03 period is special, it has a different format because the RFs and
    # all that were made for the full 05 to 50 Zd range, so there is no low and mid separation
    if [ "$Period" == "ST.03.03" ]; then
        for File in $PathIn/$Source/$Night/20*.root
        do
            echo $File
            # Yes, we're running everything twice, and if it doesn't have X range data, 
            # no lightcurve file will be created
            root -l -b -q $CodePath/EventExtractorLowZd.cpp+\(\"$File\",$DispCut,\"$Period\"\) &
            pids[${ct}]=$!
            ct+=1
        done
    else
        for File in $PathIn/05to35/$Source/$Night/20*.root
        do
            echo $File
            root -l -b -q $CodePath/EventExtractorLowZd.cpp+\(\"$File\",$DispCut,\"$Period\"\) &
            pids[${ct}]=$!
            ct+=1
        done
    fi
    
    # wait for all pids
    for pid in ${pids[*]}
    do
        wait $pid
    done
    # The intermediate files EL_night_source_run_dispcut.root have been created
    
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
    
        # Remove middle files 
        echo "Removing middle files"
        echo "rm -rf ${PathOut}/05to35/EL_${Night}_${Source}*${DispCut}.root"
        rm -rf ${PathOut}/05to35/EL_${Night}_${Source}*${DispCut}.root
    
        echo "Running FullCameraLC on low Zd range"
        
        if [ -f "${PathOut}/05to35/LightCurves_${Night}_${Source}_${DispCut}.root" ]; then
        
            echo "LightCurves file ${PathOut}/05to35/LightCurves_${Night}_${Source}_${DispCut}.root exists" 
    
            # Add the TNtuple with the events from the full camera light curve
            if [ "$Period" == "ST.03.03" ]; then
                root -l -b -q ${CodePath}/FullCameraLCCoordsTeVLZd.cpp+\(\"${PathIn}/${Source}/${Night}/20*root\",\"${PathOut}/05to35/LightCurves_${Night}_${Source}_${DispCut}.root\",\"$Source\",$DispCut,\"$Period\"\)
                if [ $? -ne 0 ]; then
                    echo "FullCameraLCCoords for LowZd, ST03.03 failed, aborting."
                    exit 1
                fi
            else
                root -l -b -q ${CodePath}/FullCameraLCCoordsTeVLZd.cpp+\(\"${PathIn}/05to35/${Source}/${Night}/20*root\",\"${PathOut}/05to35/LightCurves_${Night}_${Source}_${DispCut}.root\",\"$Source\",$DispCut,\"$Period\"\)
                if [ $? -ne 0 ]; then
                    echo "FullCameraLCCoords for LowZd, ${Period} failed, aborting."
                    exit 1
                fi
            fi
        else
            echo "Low Zd range LightCurves file does not exist. Exiting"
        fi
        
    else
        echo "No EL files in the low range. Exiting."
    fi

    ###################################################
    ##
    ##            Medium Zd range
    ##
    ##################################################
    # reset the counter for th epid array
    ct=0
    #I have to run these loops sequentially, bc the pids don't work like this
    # and CellCreator executes before some of them are finished
    if [ "$Period" == "ST.03.03" ]; then
        for File in $PathIn/$Source/$Night/20*.root
        do
            echo $File
            root -l -b -q $CodePath/EventExtractorMidZd.cpp+\(\"$File\",$DispCut,\"$Period\"\) &
            pidsM[${ct}]=$!
            ct+=1
        done
    else

        for File in $PathIn/35to50/$Source/$Night/20*.root
        do
            echo $File
            root -l -b -q $CodePath/EventExtractorMidZd.cpp+\(\"$File\",$DispCut,\"$Period\"\) &
            pidsM[${ct}]=$!
            ct+=1
        done
    fi

    # wait for all pids
    for pid in ${pidsM[*]}
    do
        wait $pid
    done
    
    echo "Number of runs processed = ${#pidsM[@]}"
    echo "Running CellCreator for the medium Zd range"
    echo "${PathOut}/35to50/EL_${Night}_${Source}*${DispCut}.root"
    
    if [ `ls -1 ${PathOut}/35to50/EL_${Night}_${Source}*${DispCut}.root 2>/dev/null | wc -l ` -gt 0 ];
    then
        echo "ok"
        for f in ${PathOut}/35to50/EL_${Night}_${Source}*${DispCut}.root
        do
            echo ${f}
            root -l -b -q ${CodePath}/CellCreator.cpp+\(\"${f}\",\"${Night}\",\"${Source}\",${DispCut}\)
        done
        # The LightCurves_night_source_dispcut.root have been created 
    
        echo "Removing middle files"
        echo "rm -rf ${PathOut}/35to50/EL_${Night}_${Source}*${DispCut}.root"
        rm -rf ${PathOut}/35to50/EL_${Night}_${Source}*${DispCut}.root
        
        echo "Running FullCameraLC on medium Zd range"
    
        if [ -f "${PathOut}/35to50/LightCurves_${Night}_${Source}_${DispCut}.root" ]; then
        
            echo "LightCurves file ${PathOut}/35to50/LightCurves_${Night}_${Source}_${DispCut}.root exists."
    
            # Add the TNtuple with the events from the full camera light curve
            if [ "$Period" == "ST.03.03" ]; then
                root -l -b -q ${CodePath}/FullCameraLCCoordsTeVMZd.cpp+\(\"${PathIn}/${Source}/${Night}/20*root\",\"${PathOut}/35to50/LightCurves_${Night}_${Source}_${DispCut}.root\",\"$Source\",$DispCut,\"$Period\"\)
                if [ $? -ne 0 ]; then
                    echo "FullCameraLCCoords for MidZd, ${Period} failed, aborting."
                    exit 1
                fi
            else
                root -l -b -q ${CodePath}/FullCameraLCCoordsTeVMZd.cpp+\(\"${PathIn}/35to50/${Source}/${Night}/20*root\",\"${PathOut}/35to50/LightCurves_${Night}_${Source}_${DispCut}.root\",\"$Source\",$DispCut,\"$Period\"\)
                if [ $? -ne 0 ]; then
                    echo "FullCameraLCCoords for MidZd, ${Period} failed, aborting."
                    exit 1
                fi
            fi 
        else
        
            echo "Medium Zd range LightCurves file does not exist. Exiting"
        
        fi
    
    else
        echo "No EL files in the medium range. Exiting."
    fi
else
    echo "Only analysis part selected"
fi

# That part runs only if we want to create the Lightcurve files from the start.
# We may want to re-run the analysis part more than once, so we put it here
# so that we have only one pipeline file

if [ -f "${PathOut}/05to35/LightCurves_${Night}_${Source}_${DispCut}.root" ]; then
        
    echo "LightCurves file ${PathOut}/05to35/LightCurves_${Night}_${Source}_${DispCut}.root exists" 
    echo "Running UltraAnalysis Low Zd range"
    
    root -l -b -q ${CodePath}/AnalysisLZdNFE.cpp+\(\"${PathOut}/05to35/LightCurves_${Night}_${Source}_${DispCut}.root\",$DispCut,${BinWidth},\"$Period\"\)
    if [ $? -ne 0 ]; then
        echo "AnalysisNFE for LowZd, ${Period} failed, aborting."
        exit 1
    fi
    
else
    echo "Low Zd range LightCurves file does not exist. Exiting"
fi

if [ -f "${PathOut}/35to50/LightCurves_${Night}_${Source}_${DispCut}.root" ]; then

    echo "LightCurves file ${PathOut}/35to50/LightCurves_${Night}_${Source}_${DispCut}.root exists."
    echo "Running UltraAnalysis Medium Zd range"

    root -l -b -q ${CodePath}/AnalysisMZdNFE.cpp+\(\"${PathOut}/35to50/LightCurves_${Night}_${Source}_${DispCut}.root\",$DispCut,${BinWidth},\"$Period\"\)
    if [ $? -ne 0 ]; then
        echo "AnalysisNFE for MidZd, ${Period} failed, aborting."
        exit 1
    fi

else
    echo "Medium Zd range LightCurves file does not exist. Exiting"
fi

