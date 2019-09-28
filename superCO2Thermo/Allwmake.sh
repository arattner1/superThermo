#!/bin/bash
 
#Build script for thermophysical libraries of supercritical carbon dioxide 
wmake libso superThermoFinal

#Build script for buoyantPimpleFoam_GAdaptCyc
wmake buoyantPimpleFoam_GAdaptCyc

#Build script for buoyantPimpleFoam_GAdaptDoubleCyc
wmake buoyantPimpleFoam_GAdaptDoubleCyc

#Build script for chtMultiRegionFoam_GAdaptCyc
wmake chtMultiRegionFoam_GAdaptCyc

#Build script for chtMultiRegionFoam_GAdaptDoubleCyc
wmake chtMultiRegionFoam_GAdaptDoubleCyc
