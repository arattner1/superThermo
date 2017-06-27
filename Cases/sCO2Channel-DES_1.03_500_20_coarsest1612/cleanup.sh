#!/bin/bash
source /storage/home/mmn5264/OpenFOAM/SetupOpenFOAM4-Personal.sh

#Cleanup script
rm constant/polyMesh/faces
rm constant/polyMesh/neighbour
rm constant/polyMesh/owner
rm constant/polyMesh/points
rm constant/polyMesh/boundary
rm -r processor*
rm -r dynamicCode
rm -r *.txt
rm -r *.dat
rm -r RunsCO2.e*
rm -r RunsCO2.o*
rm -r sCO2_*
rm -r postProcessing
