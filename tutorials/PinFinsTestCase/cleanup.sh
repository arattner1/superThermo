#!/bin/bash
source /home/mnabil/OpenFOAM/OpenFOAM-v1612+/etc/bashrc

rm -rf VTK
rm -rf constant/cellToRegion*
rm -rf 0

rm -rf constant/cellTo*
rm constant/perturbUDict
rm -rf constant/fluid/polyMesh
rm -rf constant/Pin1/polyMesh
rm -rf constant/Pin2/polyMesh
rm -rf constant/polyMesh
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


