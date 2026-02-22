# Post-processing Turbulent Water Transport in a Hydraulically Smooth Pipeline

MATLAB-based post-processing tools for Computational Fluid Dynamics (CFD) analysis of turbulent water flow in a hydraulically smooth cylindrical pipeline.

## Overview

This repository contains MATLAB scripts for analyzing CFD simulation results of turbulent water transport through a pipeline. The tools process simulation data to extract turbulent flow parameters, velocity profiles, wall shear stress, Reynolds stresses, and turbulent viscosity distributions.

## Features

- **Turbulent Flow Characterization**
  - Reynolds number calculation for turbulent regime
  - Entrance length estimation for turbulent flow
  - Y+ (dimensionless wall distance) verification and monitoring
  - Fully developed turbulent flow verification

- **Advanced Flow Analysis**
  - Wall-normal velocity profiles at multiple axial locations
  - Wall shear stress computation and distribution
  - Reynolds stress analysis
  - Turbulent viscosity field extraction
  - Total shear stress decomposition (viscous + Reynolds)

- **Data Processing**
  - Import and process CFD simulation files (.phi format)
  - Extract velocity, pressure, turbulence kinetic energy, and dissipation fields
  - Handle both Cartesian and Body-Fitted Coordinate (BFC) grids
  - Dimensionless velocity (W+) and wall distance (Y+) calculations

- **Visualization**
  - Y+ distribution along the pipe length
  - Velocity profile comparisons at different axial locations  
  - Turbulent flow development visualization

## Files

- **CFD2_Processing.m**: Main post-processing script for turbulent pipeline flow analysis
  - Calculates Reynolds number (~75,000) and entrance length
  - Processes velocity, pressure, and turbulence fields
  - Computes wall shear stress and dimensionless parameters
  - Analyzes turbulent viscosity and Reynolds stresses
  - Verifies Y+ range (30 < Y+ < 130) for proper wall treatment
  
- **XYZ_reduced_19_20.m**: Utility script for importing CFD data files
  - Imports .phi files from CFD simulations
  - Handles grid coordinates and cell data
  - Supports both Cartesian and BFC grids
  - Allows selective variable processing

- **CFD_Test_Case_Documentation.pdf**: Test case specifications and documentation

## Physical Parameters

Default simulation parameters for turbulent water flow:
- **Pipe Diameter (D)**: 0.1 m
- **Bulk Velocity (Wb)**: 0.75 m/s
- **Water Density (ρ)**: 998.23 kg/m³
- **Kinematic Viscosity (ν)**: 1.0×10⁻⁶ m²/s
- **Pipe Length (L)**: 3.0 m
- **Reynolds Number (Re)**: ~75,000 (turbulent flow)

## Usage

### CFD Data Processing

1. Update the file paths in `XYZ_reduced_19_20.m`:
   ```matlab
   directoryLoad = 'path/to/your/simulation.phi';
   directorySave = 'path/to/save/output.mat';
   ```

2. Specify which variables to process (including turbulence quantities):
   ```matlab
   variables = {'P1', 'V1', 'W1', 'DWDY', 'ENUT', 'KE', 'EP', 'STRS', 'YPLS'};
   ```

3. Run the import script:
   ```matlab
   XYZ_reduced_19_20
   ```

### Turbulent Flow Analysis

1. Ensure the processed data (`main.mat`) is available with turbulence fields

2. Run the main processing script:
   ```matlab
   CFD2_Processing
   ```

3. The script will output:
   - Reynolds number and flow regime classification
   - Entrance length estimates for turbulent flow
   - Y+ values along the pipe (verification plots)
   - Velocity profiles at multiple locations (2m, 3m, 4m, 5m, 6m)
   - Wall shear stress distribution
   - Reynolds stresses and turbulent viscosity fields

## Output Variables

After processing, the following variables are available in the workspace:

### Grid Information
- `NY`, `NZ`: Number of cells in radial and axial directions
- `Y`, `Z`: Cell center coordinates
- `R`, `D`: Pipe radius and diameter

### Flow Fields
- `P`: Pressure field
- `W`: Axial velocity field
- `DWDY`: Axial velocity gradient in radial direction
- `DPDZ`: Pressure gradient in axial direction

### Turbulence Quantities
- `K`: Turbulence kinetic energy field
- `EPS`: Turbulence dissipation rate field
- `ENUT`: Eddy (turbulent) viscosity field
- `mu_turb`: Turbulent dynamic viscosity (ρ × ENUT)

### Wall Parameters
- `tau_wall`: Wall shear stress distribution
- `W_tau`: Friction velocity (u*)
- `YPLS1`: Y+ at first cell from wall
- `YPLS`: Full Y+ field
- `WPLS`: Dimensionless velocity (W+)

### Stress Components
- `tau_visc`: Viscous shear stress
- `tau_re`: Reynolds (turbulent) shear stress
- `tau_tot`: Total shear stress (viscous + Reynolds)

## Turbulence Modeling

The analysis supports standard turbulence models commonly used in CFD:
- k-ε (k-epsilon) model
- k-ω (k-omega) model
- Reynolds Stress Models (RSM)

Proper wall treatment is verified through Y+ monitoring (target range: 30-130 for wall functions).

## Requirements

- MATLAB (tested on R2019a or later)
- Sufficient memory for 3D turbulent field data processing

## Applications

This toolset is designed for:
- Turbulent flow research and analysis
- CFD simulation validation against experimental data
- Pipeline hydraulics and design
- Turbulence model performance assessment
- Educational purposes in advanced fluid mechanics and CFD

## Technical Notes

- **Y+ Verification**: The code checks that 30 < Y+ < 130 for proper wall function application
- **Fully Developed Flow**: Velocity profiles are extracted at multiple locations to verify flow development
- **Reynolds Stress**: Decomposition of total stress into viscous and turbulent components
- **Hydraulically Smooth**: Pipeline treated as hydraulically smooth (no roughness effects)

## License

This project is licensed under the GPL-3.0 License - see the [LICENSE](LICENSE) file for details.

## Author

Ziad Kandil

## Acknowledgments

Developed as part of MSc Mathematical Engineering coursework in Computational Fluid Dynamics, focusing on turbulent flow analysis and advanced post-processing techniques.
