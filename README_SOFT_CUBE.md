# SOFT CUBE

**SOFT CUBE** (Synthesis of Flows and Temperatures based on a CFD model and urban building environments) is a computational framework designed to generate high-resolution 3D meteorological fields in urban areas.  
It reduces the computational cost of CFD simulations by reconstructing wind and air temperature fields using a pre-computed scenario database integrated with LDAPS forecasts.

---

## Overview

Urban flow and temperature distributions are influenced by:

- The vertical structure of background wind speed, wind direction, and air temperature entering the target area.
- Local surface characteristics (e.g., buildings, terrain) and surface heating.

SOFT CUBE separates the domain into two layers:

- **Roughness Dominant Layer (RL):** Flow and temperature changes are primarily governed by surface roughness.
- **Background Dominant Layer (BL):** Turbulent mixing dominates vertical distributions of momentum and energy.

---

## Computational Domain

- Size: **2000 m × 2000 m × 500 m** (x, y, z)
- Resolution: **10 m** (horizontal), **5 m** (vertical)
- Grid points: **200 × 200 × 100**

---

## Main Database

The main database consists of **6,400** CFD simulations combining:

- **8** inflow wind speed conditions
- **32** inflow wind directions
- **5** ratios of surface-to-air temperature
- **5** land-cover types (roof, wall, road, green, soil)

**Variables stored per grid point:**  
- Wind velocity components (u, v, w)  
- Air temperature change rates (rt) per land-cover type

---

## Auxiliary Databases

### 1. Vertical Variations in Wind and Temperature
To represent BL dynamics:
- **Vertical temperature gradients (γ):** –20 K/km to 20 K/km (10 K/km intervals) → 40 datasets.
- **Vertical wind speed gradients (σ):** –10 m/s/km to 30 m/s/km (10 m/s/km intervals).
- **Wind direction differences (λ):** 0° to 180° (45° intervals) → 200 datasets.

---

## Source Code

The SOFT CUBE main solver is implemented in **Fortran**.  
Compile using a standard Fortran compiler (e.g., `gfortran`):


## Software and Data Availability

- **Name of software:** SOFT-CUBE  
- **Developers:** Jang-Woon Wang and Jae-Jin Kim  
- **Contact:** jwwang@pukyong.ac.kr  
- **Date first available:** August 2025  
- **Program language:** Fortran  
- **Source code:** [GitHub Repository](https://github.com/jwwang-gif/SOFT-CUBE/tree/main)  
- **Data availability:** The full CFD database (~several TB) is hosted on an internal HPC server subject to institutional security policies.  
  - Public hosting of the full dataset is not possible.  
  - The software requires the complete dataset to operate; partial datasets are not supported.  
  - Source code and metadata are provided in the repository.  
  - Full access can be granted via a formal data-sharing agreement, subject to institutional approval.

---
