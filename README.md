# superCO2Thermo
CFD simulation platform for high-heat-flux turbulent supercritical CO2 flow 

Copyright 2019: Mahdi Nabil, Alexander S Rattner

[Multiscale Thermal Fluids and Energy (MTFE) Laboratory](http://sites.psu.edu/mtfe/)

The Pennsylvania State University

## Overview
superCO2Thermo is an extensible open source volume-of-fluid (VOF) based simulation tool for studying convection heat transfer in high-heat-flux turbulent supercritical CO2 flow. The solver approach is based on the buoyantPimpleFoam and chtMultiRegionFoam codes developed by [OpenCFD Ltd.](http://openfoam.com/). Target applications for superCO2Thermo include:

* Channel/pipe sCO2 flow
* Enhanced heat exchangers working with sCO2, e.g., pin-fins enhanced geometry, etc.

If you use this solver in a project or scholarly work, we ask that you include a citation for [Nabil and Rattner (2019)](https://www.sciencedirect.com/science/article/pii/S0017931019308579). 

## Installation
The current version of the code uses the [OpenFOAM v1612+ libraries](https://www.openfoam.com/releases/openfoam-v1612+/).

**To Install:**
Navigate to a working folder in a shell terminal, clone the git code repository, and build.
```
$ git clone https://github.com/arattner1/superThermo
$ cd superThermo/superCO2Thermo
$ source Allwmake.sh
```

## Contribute
MTFE welcomes collaboration with other investigators studying supercritical flows. Please [contact us](mailto:Alex.Rattner@psu.edu) if you are interested in expanding the solver or find bugs to correct. We can also provide limited support (on a case-by-case basis) or consulting servies.

## Acknowledgements
This research was generously supported, in part, by the US National Science Foundation (NSF), and computing resources from NSF XSEDE. 

## References
* Nabil, M., Rattner, A.S., 2019. Large eddy simulations of high-heat-flux supercritical CO2 convection in microchannels: Mixed convection and non-uniform heating. International Journal of Heat and Mass Transfer 145: 118710. [DOI: 10.1016/j.ijheatmasstransfer.2019.118710](https://www.sciencedirect.com/science/article/pii/S0017931019308579).

