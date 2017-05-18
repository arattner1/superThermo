#!/bin/bash
source /storage/home/mmn5264/OpenFOAM/SetupOpenFOAM-Personal.sh
cd ${0%/*} || exit 1    # run from this directory

# Source tutorial clean functions
. $WM_PROJECT_DIR/bin/tools/CleanFunctions

#Cleanup script
cleanTimeDirectories
rm constant/polyMesh/faces
rm constant/polyMesh/neighbour
rm constant/polyMesh/owner
rm constant/polyMesh/points
rm constant/polyMesh/boundary
rm -r processor*
rm -r dynamicCode
