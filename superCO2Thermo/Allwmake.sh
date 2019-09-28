#!/bin/bash
 
#Build script for buoyantPimpleFoam_GAdaptCyc
wmake libso superThermoFinal
wmake buoyantPimpleFoam_GAdaptCyc

#Build script for buoyantPimpleFoam_GAdaptDoubleCyc
wmake buoyantPimpleFoam_GAdaptDoubleCyc

#Build script for chtMultiRegionFoam_GAdaptCyc
wmake chtMultiRegionFoam_GAdaptCyc

#Build script for chtMultiRegionFoam_GAdaptDoubleCyc
wmake chtMultiRegionFoam_GAdaptDoubleCyc
