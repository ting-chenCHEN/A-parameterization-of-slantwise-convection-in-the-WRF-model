# A-parameterization-of-slantwise-convection-in-the-WRF-model

Authors: Ting-Chen Chen, Man-Kong (Peter) Yau and Daniel J. Kirshbaum

Department of Atmospheric and Oceanic Sciences, McGill University 

Montreal, Canada

Submitted to Journal of Atmospheric Sciences

This repository contains all the required modification of the source codes in the Weather Research and Forecasting (WRF) model (v4.1.5) for implementing a parameterization of slantwise convection (SC) developed by the present authors following a PhD work of L. Ma (2000). Physically, this SC scheme operates in locally-defined 2D cross-sections perpendicular to the deep-layer-averaged thermal winds. Its central concept is to eliminate positive slantwise convective available potential energy (SCAPE) by adjusting the momentum field. The environment is adjusted toward slantwise neutrality given a prescribed adjustment timescale (3-5 hours). Condensational heating and moisture loss associated with the upward motion are also parameterized. This scheme is a supplement to the existing cumulus parameterization scheme (CPS) for upright convection.

The main program of the SC scheme is named "module_cu_slantwise.F" under the folder "phys." 

To switch on the SC scheme in the simulation, add a variable "slcu_physics = 1" in the namelist.input 
within the &physics namelist record, e.g.,
 
 &physics   
 physics_suite                       = 'CONUS'      
 mp_physics                          = -1,    
 cu_physics                          =  1,   
 slcu_physics                        =  1,    


Note that this scheme has been coded successfully for OpenMP, and the test results show similar scaling performance as runs using only the Kain-Fritsch CPS. The implementation for MPI is not available at the moment.

To download and use this scheme, please cite 

https://doi.org/10.5281/zenodo.5602111

Chen, T.-C., M. K. Yau, and D. J. Kirshbaum, 2021: A parameterization of slantwise convection in the WRF model. Submitted to J. Atmos. Sci.


*The Advanced Research WRF model can be downloaded from https://www2.mmm.ucar.edu/wrf/users/download/get_sources.html (last access: April 2020).
