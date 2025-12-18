#!/bin/bash
# This script runs the acceptance map creator and the Nominal effective area calculation
# that is added to the acceptance map file
# This script has to be run for each period before anything else

Period=$1

PathIn="/data/magic/users-ifae/edosouto/MAGIC/BigSearch/Melibea/${Period}"
CodePath="/data/magic/users-ifae/edosouto/MAGIC/BigSearch/Code/MARS3version"


FILE="${PathIn}/AcceptanceMaps_${Period}.root"

if [ ! -f "$FILE" ]; then
    echo "$FILE does not exist. Running command..."
    echo "If you're running multiple jobs with this, it will fail"
    root -l -b -q ${CodePath}/AcceptanceMapsGen.cpp+\(\"${Period}\"\)
    if [ $? -ne 0 ]; then
        echo "Acceptance map generation failed, aborting."
        exit 1
    fi
    root -l -b -q ${CodePath}/NominalAeffComp.cpp\(\"${Period}\"\)
    status=$?
    if [ $status -ne 0 ]; then
        case $status in
            1) echo "Error: MC input file not found";;
            2) echo "Error: Can't open the txt file to write epsilon_ref";;
            *) echo "Unknown error, code $status";;
        esac
        exit $status
    fi
else
    echo "$FILE already exists. Skipping."
fi
