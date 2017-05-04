#!/bin/bash
source /storage/work/asr20/OpenFOAM/SetupOpenFOAM.sh
m4 constant/polyMesh/blockMeshDict.m4 > constant/polyMesh/blockMeshDict
blockMesh
mirrorMesh -noFunctionObjects -overwrite
checkMesh
decomposePar -force
