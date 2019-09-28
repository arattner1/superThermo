#!/bin/bash 

m4 system/blockMeshDict.m4 > system/blockMeshDict
blockMesh
#checkMesh
#cp -r orig0 1
#mapFields ../Case_1.1_500_40_05_NoHT_CoarseCyclic -mapMethod cellPointInterpolate
#m4 constant/perturbUDict1.m4 > constant/perturbUDict
#perturbUChannel
#rm constant/perturbUDict
#m4 constant/perturbUDict2.m4 > constant/perturbUDict
#perturbUChannel
decomposePar -force 
exit
