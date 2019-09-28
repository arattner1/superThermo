#!/bin/sh

PROB_DIR=Case_1.1_300_H3_PF_Fine
MAP_DIR=Case_1.1_300_H3_PF_Coarse
MAP_TIME=0.0145

source /storage/home/mmn5264/OpenFOAM/SetupOpenFOAM4-Personal.sh

m4 system/blockMeshDict.m4 > system/blockMeshDict
blockMesh
#checkMesh
cp -r orig0 $MAP_TIME
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

# mapping fields only in fluid domain
mapFields ../$MAP_DIR -sourceRegion fluid -targetRegion fluid -mapMethod cellPointInterpolate

# mapping fields only in Pin1 domain
mapFields ../$MAP_DIR -sourceRegion Pin1 -targetRegion Pin1 -mapMethod cellPointInterpolate

# mapping fields only in Pin1 domain
mapFields ../$MAP_DIR -sourceRegion Pin2 -targetRegion Pin2 -mapMethod cellPointInterpolate

# remove fluid fields from solid regions (important for post-processing)
for i in Pin1 Pin2
do
   rm -f 0*/$i/{nut,alphat,U,p_rgh}
done

# Delete extra files 
rm -r $MAP_TIME/U*
rm -r $MAP_TIME/T*
rm -r $MAP_TIME/p*
rm -r $MAP_TIME/cell*
rm -r $MAP_TIME/alpha*
rm -r $MAP_TIME/nut*
rm -r $MAP_TIME/Pin1/p_rgh*
rm -r $MAP_TIME/Pin1/nut*
rm -r $MAP_TIME/Pin1/alpha*
rm -r $MAP_TIME/Pin1/U*
rm -r $MAP_TIME/Pin2/p_rgh*
rm -r $MAP_TIME/Pin2/nut*
rm -r $MAP_TIME/Pin2/alpha*
rm -r $MAP_TIME/Pin2/U*

# Decompose
for i in fluid Pin1 Pin2
do
decomposePar -force -region $i 
done
#------------------------------------------------------------------------------
