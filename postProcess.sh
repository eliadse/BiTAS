#!/bin/bash

period=$1
# Remove all periods
period_nodots=${period//./}

# Merge the results files of the period
hadd -f /data/magic/users-ifae/edosouto/MAGIC/BigSearch/Code/MARS3version/${period}/merged_${period_nodots}_100s.root /data/magic/users-ifae/edosouto/MAGIC/BigSearch/Code/MARS3version/${period}/*/Results/Results*100s.root

# Add the effective area to the merged file
root -l -b -q /data/magic/users-ifae/edosouto/MAGIC/BigSearch/Code/MARS3version/RingCollectionArea_LZd.cpp+\(\"${period}\"\)
root -l -b -q /data/magic/users-ifae/edosouto/MAGIC/BigSearch/Code/MARS3version/RingCollectionArea_MZd.cpp+\(\"${period}\"\)
