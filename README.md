**Package for MRI reconstruction. Supports:**
  - non-Cartesian regridding
  - iterative SENSE (Cartesian or non-Cartesian)
  - iterative SENSE for higher order models that may include: 
    - a B0 map
    - time varying spherical harmonics of phase accrual
    
**Before Usage**
  - run setPath.m to add all the directories to the Matlab path
    
**Includes a basic image viewer and ROI drawing tool called "bview"**
  - enter "bview" in the Matlab command prompt to use
    
**Tips**
  - step through the demos in the demos folder
    - recommended order is demo_regridding, demo_regSENSE, demo_highOrder

**Notes**
  - Matlab gpuArray functionality is used whenever possible, so a compatible GPU is strongly recommended

**Acknowledgement**
  - Please reference the following works for usage of the below functions:
    - nufftOp: Baron CA, Dwork N, Pauly JM, Nishimura DG. Rapid compressed sensing reconstruction of 3D non-Cartesian MRI. Magn. Reson. Med. 2018;79:2685–2692.

(c) 2020, Corey Baron
