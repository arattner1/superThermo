#!/bin/sh
source /home/mnabil/OpenFOAM/OpenFOAM-v1612+/etc/bashrc

m4 system/blockMeshDict.m4 > system/blockMeshDict
blockMesh
checkMesh
cp -r orig0 0
topoSet
splitMeshRegions -cellZones -overwrite

# remove fluid fields from solid regions (important for post-processing)
for i in Pin1 Pin2
do
   rm -f 0*/$i/{nut,alphat,U,p_rgh}
done

for i in fluid Pin1 Pin2
do
   changeDictionary -region $i 
done

# remove fluid fields from solid regions (important for post-processing)
for i in Pin1 Pin2
do
   rm -f 0*/$i/{nut,alphat,U,p_rgh}
done

# Decompose
for i in fluid Pin1 Pin2
do
decomposePar -force -region $i 
done

# Run
mpirun -np 8 chtMultiRegionFoam_GAdaptDoubleCyc -parallel 

