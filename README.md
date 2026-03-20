## MatMRI

**MatMRI** is a GPU-enabled MATLAB package for advanced MRI reconstruction, modeling, and quantitative analysis. It is designed for flexible handling of both Cartesian and non-Cartesian acquisitions, with a strong emphasis on model-based reconstruction and field-monitoring–informed workflows.

**Package supports:**
- **Reconstruction**
  - non-Cartesian regridding
  - iterative SENSE (Cartesian or non-Cartesian)
  - iterative SENSE for higher-order signal models including:
    - B0 inhomogeneity correction
    - time-varying phase modeled with spherical harmonics (e.g., field probe–based encoding)
- **Diffusion MRI fitting**
  - spatially regularized diffusion kurtosis fitting with an axially symmetric model (`nii2kurt.m`)
  - free water–corrected kurtosis and micro-FA (`nii2uFA_fwe.m`)  
    - Cite: doi: 10.3389/fnins.2023.1074730

---

## Brainhack 2026 – UWO Project

This repository is being extended as part of **Brainhack 2026 (Western University)**.

### Project Goal
The objective is to significantly improve the existing MATLAB-based image viewer **`bview`** to better support modern fMRI analysis and trajectory validation workflows, particularly for non-Cartesian acquisitions (e.g., spiral imaging with field monitoring).

---

## Current Tool: bview

MatMRI includes a lightweight image viewer and ROI tool:

- Launch with:
  ```matlab
  bview
  ```

- Current capabilities:
  - Image visualization
  - ROI drawing
  - Basic interaction

---

## Planned Improvements

### fMRI Functionality

#### 1. Interactive Voxel Time Series
- Click on a voxel in the image → automatically display its time series
- Support **multi-dataset comparison**:
  - Simultaneous plotting of time courses from ≥2 datasets (e.g., spiral vs EPI, reconstruction variants)

#### 2. Time Series Processing
Add inline preprocessing options:
- Demeaning
- Normalization (e.g., z-score or max scaling)
- Percent signal change conversion

#### 3. Activation Map Overlay
- Overlay statistical maps (e.g., GLM t-maps) on anatomical/functional images
- Adjustable thresholding
- Colormap support with **dynamic colorbar display**
- Transparency control for overlay blending

#### 4. Automated SNR and tSNR Calculation
- Compute SNR (Signal-to-Noise Ratio) from a single volume: mean signal in a user-defined ROI divided by the standard deviation of a background/noise region
- Compute tSNR (Temporal SNR) across the time series: mean signal divided by its temporal standard deviation, computed voxel-wise
- Generate and display tSNR maps as an overlay on the reference image
- Report summary statistics (mean, median, min/max SNR/tSNR) across the whole brain or within a selected ROI

---

### Field Monitoring / Trajectory Validation

#### processSkope Integration
- `processSkope` outputs `.mat` files containing processed MRI field probe measurements (e.g., spherical harmonic coefficients over time)

#### New Feature: k-space Trajectory Visualization
- Load `processSkope` output directly in `bview`
- Reconstruct and display **measured k-space trajectories**
- Enable:
  - Visual validation of trajectory fidelity
  - Comparison against nominal trajectories
  - Identification of trajectory errors or drift

---

## Before Usage
- Run:
  ```matlab
  setPath.m
  ```
  to add all required directories to the MATLAB path.

---

## Tips
- Step through demos in `/demos`:
  - Recommended order:
    1. `demo_regridding`
    2. `demo_regSENSE`
    3. `demo_highOrder`
- Inspect function files for usage details and implementation specifics

---

## Notes
- Uses MATLAB `gpuArray` wherever possible  
  → A compatible GPU is strongly recommended for reconstruction workflows

---

## Acknowledgement

**Cite as:**
1. Varela-Mattatall G, Dubovan PI, Santini T, Gilbert KM, Menon RS, Baron CA.  
   *Single-shot spiral diffusion-weighted imaging at 7T using expanded encoding with compressed sensing.*  
   Magn Reson Med. 2023. doi: 10.1002/mrm.29666  

2. Baron CA (2021).  
   *MatMRI: A GPU enabled package for model based MRI image reconstruction.*  
   Zenodo. http://doi.org/10.5281/zenodo.4495476  

**Function-specific references:**
- **nii2kurt**: Hamilton et al., Imaging Neuroscience (2024)  
- **nii2uFA_fwe**: Arezza et al., Front Neurosci (2023)  
- **harmonicsFromRaw**: Dubovan et al., Magn Reson Med (2023)  
- **findDelAuto**: Dubovan et al., Magn Reson Med (2022)  
- **nufftOp**: Baron et al., Magn Reson Med (2018)  
- **sampHighOrder**:  
  - Wilm et al., IEEE TMI (2012)  
  - Baron et al., Magn Reson Med (2018)

---

(c) 2020, Corey Baron
