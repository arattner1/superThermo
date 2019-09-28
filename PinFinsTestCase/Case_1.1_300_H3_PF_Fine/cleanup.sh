#!/bin/bash
source /storage/home/mmn5264/OpenFOAM/SetupOpenFOAM4-Personal.sh

cd ${0%/*} || exit 1    # Run from this directory

# Source tutorial clean functions
. $WM_PROJECT_DIR/bin/tools/CleanFunctions

rm -rf VTK
rm -rf constant/cellToRegion*
rm -rf 0

rm -rf constant/cellTo*
rm constant/perturbUDict
rm -r cht*
rm -r processor*
rm -r dynamicCode
rm -r LES_*
rm -r Short*
rm -r log.*
rm -r Inlet*
rm -r Outlet*
rm -r Test*
rm -r CHT*
rm -r block*
rm -r post*
rm -r CHT*
rm -r buo*
rm -r F300*
rm -rf constant/fluid/polyMesh
rm -rf constant/Pin1/polyMesh
rm -rf constant/Pin2/polyMesh
rm -rf constant/polyMesh

