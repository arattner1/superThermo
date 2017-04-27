#!/bin/bash
m4 constant/polyMesh/blockMeshDict.m4 > constant/polyMesh/blockMeshDict
blockMesh
mirrorMesh -noFunctionObjects -overwrite
checkMesh
decomposePar -force
mpirun -np 4 buoyantPimpleFoam -parallel
#buoyantPimpleFoam
