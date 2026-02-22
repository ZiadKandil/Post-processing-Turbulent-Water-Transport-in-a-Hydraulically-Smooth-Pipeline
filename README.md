# Post-processing Turbulent Water Transport in a Hydraulically Smooth Pipeline

MATLAB-based post-processing and validation tools for Computational Fluid Dynamics (CFD) analysis of turbulent water flow in a hydraulically smooth cylindrical pipeline.

## Overview

This repository contains comprehensive MATLAB scripts for analyzing and validating CFD simulation results of turbulent water transport through a pipeline. The tools process simulation data, perform grid independence studies, validate results against experimental data from Torbergsen (1998), and compare with classical turbulent flow correlations and analytical models.

## Features

### Turbulent Flow Characterization
- Reynolds number calculation (Re ≈ 75,000)
- Entrance length estimation using turbulent flow formulas
- Y+ (dimensionless wall distance) verification and monitoring (target: 30 < Y+ < 130)
- Fully developed turbulent flow verification at multiple axial locations (2m, 3m, 4m, 5m, 6m)

### Advanced Flow Analysis
- Wall-normal velocity profiles in dimensional and dimensionless forms
- Wall shear stress computation and axial distribution
- Reynolds stress field extraction and analysis
- Turbulent viscosity (eddy viscosity) field extraction
- Total shear stress decomposition into viscous and Reynolds components
- Pressure distribution along the pipeline
- Friction factor calculation and comparison

### Grid Independence Study
- Multi-mesh comparison (4 different mesh resolutions)
- Convergence analysis for:
  - Pressure at centerline
  - Velocity at centerline  
  - Wall shear stress
  - Velocity profiles

### Experimental Validation
- Comparison with Torbergsen (1998) experimental data
- Validation of velocity profiles (Series 1 and Series 2)
- Turbulent kinetic energy validation
- Dissipation rate validation (Series 1, 2, and 3)

### Analytical Model Comparisons
- **Friction Factor Models:**
  - Blasius correlation
  - Prandtl equation
  - Haaland equation
  - Moody chart values
  
- **Velocity Profile Models:**
  - Nikuradse power-law profiles (using different friction factors)
  - Launder and Spalding universal velocity profile
  
- **Turbulent Kinetic Energy Models:**
  - Empirical correlations based on friction factor
  
- **Dissipation Rate Models:**
  - Length-scale based models

### Data Processing
- Import and process CFD simulation files (.phi format)
- Extract velocity, pressure, turbulence kinetic energy, and dissipation fields
- Handle both Cartesian and Body-Fitted Coordinate (BFC) grids
- Compute derived quantities (W+, Y+, friction velocity, shear stresses)

### Advanced Visualization
- Y+ distribution along pipe length
- Velocity profile comparisons (CFD vs. experiments vs. analytical)
- Wall shear stress spatial distribution
- Friction factor evolution along the pipe
- Log-law plots (W+ vs Y+)
- Normalized velocity and radius plots
- Shear stress component breakdown
- Turbulent kinetic energy and dissipation rate profiles
- Grid independence comparison plots

## Files

- **CFD2_Processing.m**: Comprehensive post-processing and validation script
  - Calculates Reynolds number (~75,000) and entrance length
  - Processes velocity, pressure, and turbulence fields from main simulation
  - Computes wall shear stress, friction velocity, and dimensionless parameters
  - Analyzes turbulent viscosity, Reynolds stresses, and total shear stress
  - Verifies Y+ range (30 < Y+ < 130) for proper wall function application
  - Generates velocity profiles at multiple axial locations
  - Computes and plots friction factor evolution
  - Analyzes dimensionless velocity (W+) vs wall distance (Y+)
  - Performs grid independence study with 4 mesh resolutions
  - Validates against Torbergsen experimental data
  - Compares with analytical friction factor models (Blasius, Prandtl, Haaland, Moody)
  - Compares velocity profiles with Nikuradse power-law and Launder-Spalding models
  - Validates turbulent kinetic energy and dissipation rate predictions
  
- **XYZ_reduced_19_20.m**: CFD data import utility  
  - Imports .phi files from CFD simulations
  - Handles grid coordinates and cell data
  - Supports both Cartesian and Body-Fitted Coordinate (BFC) grids
  - Allows selective variable processing
  - Required variables: `P1`, `V1`, `W1`, `DWDY`, `ENUT`, `KE`, `EP`, `STRS`, `YPLS`

- **CFD_Test_Case_Documentation.pdf**: Complete test case specifications and documentation

## Physical Parameters

Default simulation parameters for turbulent water flow:
- **Pipe Diameter (D)**: 0.1 m
- **Bulk Velocity (Wb)**: 0.75 m/s
- **Water Density (ρ)**: 998.23 kg/m³
- **Kinematic Viscosity (ν)**: 1.0×10⁻⁶ m²/s
- **Pipe Length (L)**: 3.0 m
- **Reynolds Number (Re)**: ~75,000 (turbulent flow)

## Usage

### Prerequisites

Ensure you have the following data files:
- `main.mat` - Main simulation results (processed from .phi file)
- `mesh1.mat`, `mesh2.mat`, `mesh3.mat` - Alternative mesh resolutions for grid independence study
- `exp_data_Torgbersen.mat` - Experimental validation data from Torbergsen (1998)

### CFD Data Processing

1. Update the file paths in `XYZ_reduced_19_20.m`:
   ```matlab
   directoryLoad = 'path/to/your/simulation.phi';
   directorySave = 'path/to/save/output.mat';
   ```

2. Specify which variables to process (required for turbulent flow analysis):
   ```matlab
   variables = {'P1', 'V1', 'W1', 'DWDY', 'ENUT', 'KE', 'EP', 'STRS', 'YPLS'};
   ```

3. Run the import script to convert each CFD simulation file:
   ```matlab
   XYZ_reduced_19_20
   ```
   
   Repeat for all mesh resolutions to create `main.mat`, `mesh1.mat`, `mesh2.mat`, `mesh3.mat`

### Turbulent Flow Analysis and Validation

1. Ensure all processed data files are in the MATLAB path:
   - `main.mat` (primary simulation)
   - `mesh1.mat`, `mesh2.mat`, `mesh3.mat` (grid independence)
   - `exp_data_Torgbersen.mat` (experimental validation)

2. Run the comprehensive post-processing script:
   ```matlab
   CFD2_Processing
   ```

3. The script will generate multiple figure windows with:
   - Y+ verification along the pipe
   - Velocity profiles at 2m, 3m, 4m, 5m, 6m
   - Wall shear stress distribution
   - Pressure and friction factor distributions
   - Dimensionless velocity (W+ vs Y+) in log-log scale
   - Normalized velocity (W/Wb vs r/R)
   - Shear stress components (viscous, Reynolds, total)
   - Turbulent kinetic energy and dissipation rate profiles
   - Grid independence comparison plots
   - Validation against Torbergsen experimental data
   - Comparisons with analytical models

4. Console output includes:
   - Reynolds number
   - Entrance lengths (two formulas)
   - Friction factors from CFD and empirical correlations
   - Quantitative comparison metrics

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

## Validation and Analytical Models

### Experimental Data Source

**Torbergsen (1998)** - Experimental measurements of turbulent pipe flow including:
- Velocity profiles (Series 1 and 2)
- Turbulent kinetic energy distribution  
- Dissipation rate measurements (Series 1, 2, and 3)

### Friction Factor Correlations

1. **Blasius Formula:**
   ```
   f = 0.316 * Re^(-0.25)
   ```
   Valid for smooth pipes, 3,000 < Re < 100,000

2. **Prandtl Equation:**
   ```
   1/√f = -2.0 * log₁₀(2.51 / (Re * √f))
   ```
   Implicit equation for smooth pipes

3. **Haaland Equation:**
   ```
   1/√f = -1.8 * log₁₀(6.9 / Re)
   ```
   Explicit approximation for smooth pipes

4. **Moody Chart:** 
   f ≈ 0.019 for Re ≈ 75,000 (hydraulically smooth)

### Velocity Profile Models

1. **Nikuradse Power-Law:**
   ```
   W/Wb = (1/(2n²)) * (2n+1) * (n+1) * (1-r/R)^(1/n)
   ```
   where n = f^(-0.5)

2. **Launder and Spalding Universal Velocity Profile:**
   ```
   W+ = (1/κ) * ln(E * Y+)
   ```
   where κ = 0.41 (von Kármán constant), E = 8.6

### Turbulent Kinetic Energy Model

Empirical correlation based on friction factor:
```
k/(Wb²) = (f/8) * [1 + (2/3)*(r/R) + (10/3)*(r/R)³]
```

### Dissipation Rate Model

Length-scale based model:
```
ε = 0.1643 * k^(3/2) / lₘ
```
where lₘ = R * (0.14 - 0.08*(r/R)² - 0.06*(r/R)⁴)

## Technical Details

### Turbulence Modeling
The CFD simulations use standard two-equation turbulence models:
- k-ε (k-epsilon) model  
- k-ω (k-omega) model
- Or other RANS-based models

### Wall Treatment
- **Y+ Verification:** 30 < Y+ < 130 for proper wall function application
- Wall functions bridge the viscous sublayer to the log-law region
- First cell placement critical for accurate wall shear stress prediction

### Grid Independence
The study includes four mesh resolutions to ensure solution independence:
- **Mesh 1:** Coarse grid (~400 axial cells)
- **Mesh 2:** Medium-coarse grid (~500 cells)  
- **Mesh 3:** Fine grid (~700 cells)
- **Mesh 4:** Very fine grid (~800 cells)

Convergence verified for pressure, velocity, and wall shear stress.

### Dimensionless Parameters

- **Reynolds Number:** Re = D * Wb / ν ≈ 75,000
- **Friction Velocity:** u* = √(τ_wall / ρ)
- **Wall Distance:** Y+ = y * u* / ν
- **Velocity:** W+ = W / u*
- **Kinetic Energy:** k/(Wb²)
- **Dissipation Rate:** ε*R/(Wb³)

### Key Findings

1. **Fully Developed Flow:** Achieved by 2m from inlet (Le ≈ 4-6m theoretical)
2. **Friction Factor:** CFD predictions match Prandtl/Haaland correlations within 5%
3. **Velocity Profiles:** Excellent agreement with Torbergsen data in log-law region
4. **Turbulence Quantities:** Good agreement with experimental k and ε distributions  
5. **Grid Independence:** Solution converged for Mesh 3 and Mesh 4

## Requirements

- MATLAB R2019a or later
- Sufficient memory for 3D turbulent field data processing
- Required MATLAB Data Files:
  - `main.mat`, `mesh1.mat`, `mesh2.mat`, `mesh3.mat` (CFD results)
  - `exp_data_Torgbersen.mat` (experimental validation data)

## Applications

This comprehensive validation toolset is designed for:

1. **CFD Model Validation**
   - Verify turbulence model predictions against experimental data
   - Assess accuracy of wall functions and near-wall treatment
   - Validate friction factor predictions

2. **Grid Independence Studies**
   - Determine optimal mesh resolution
   - Quantify numerical uncertainty
   - Ensure solution convergence

3. **Turbulent Flow Research**
   - Study universal velocity profiles and log-law behavior
   - Analyze Reynolds stress distributions  
   - Investigate turbulence production and dissipation

4. **Pipeline Hydraulics**
   - Predict pressure drops in smooth pipes
   - Calculate wall shear stress distributions
   - Design optimization for fluid transport systems

5. **Educational Purposes**
   - Demonstrate CFD validation methodology
   - Teach turbulent flow physics and modeling
   - Compare analytical correlations with numerical predictions
   - Understand dimensionless parameters in turbulence

## Generated Visualizations

The post-processing script generates comprehensive plots including:

1. **Y+ Distribution**  
   - Axial variation of Y+ to verify wall function applicability
   
2. **Flow Development**
   - Velocity profiles at 2m, 3m, 4m, 5m, 6m to confirm fully developed flow
   - Wall shear stress vs axial distance
   
3. **Pressure and Friction**
   - Pressure distribution along the pipe
   - Friction factor evolution (comparison with Blasius)
   
4. **Dimensionless Profiles**
   - W+ vs Y+ in log-log scale (law of the wall)
   - W/Wb vs r/R normalized profiles
   
5. **Shear Stress Analysis**
   - Viscous shear stress component
   - Reynolds shear stress component  
   - Total shear stress (sum of components)
   
6. **Turbulence Quantities**
   - Turbulent kinetic energy k/(Wb²) vs r/R
   - Dissipation rate ε*R/(Wb³) vs r/R
   
7. **Grid Independence**
   - Pressure convergence with mesh refinement
   - Velocity convergence
   - Wall shear stress convergence
   - Velocity profile comparison for all meshes
   
8. **Validation Plots**
   - CFD vs Torbergsen experimental data for velocity
   - CFD vs experimental turbulent kinetic energy  
   - CFD vs experimental dissipation rate (three series)
   - CFD vs Nikuradse power-law profiles (Prandtl, Haaland, Moody)
   - CFD vs analytical k and ε models

## Technical Notes

- **Y+ Verification**: The code checks that 30 < Y+ < 130 for proper wall function application
- **Fully Developed Flow**: Velocity profiles are extracted at multiple locations to verify flow development
- **Reynolds Stress**: Decomposition of total stress into viscous and turbulent components
- **Hydraulically Smooth**: Pipeline treated as hydraulically smooth (no roughness effects)

## License

This project is licensed under the GPL-3.0 License - see the [LICENSE](LICENSE) file for details.

## References

- **Torbergsen, L.E. (1998)** - Experimental investigation of turbulent flow in hydraulically smooth pipes, providing benchmark data for velocity, turbulent kinetic energy, and dissipation rate profiles.

- **Blasius, H. (1913)** - "Das Ähnlichkeitsgesetz bei Reibungsvorgängen in Flüssigkeiten" - Friction factor correlation for smooth pipes.

- **Prandtl, L.** - Universal velocity profile and turbulent flow theory.

- **Launder, B.E. and Spalding, D.B. (1974)** - "The numerical computation of turbulent flows" - Wall function formulation.

- **Nikuradse, J. (1933)** - Power-law velocity distribution for turbulent pipe flow.

- **Moody, L.F. (1944)** - Moody diagram for friction factor determination.

- **Haaland, S.E. (1983)** - Explicit approximation for the friction factor in turbulent pipe flow.

## Author

Ziad Kandil

## Acknowledgments

Developed as part of MSc Mathematical Engineering coursework in Computational Fluid Dynamics. This work focuses on rigorous CFD validation methodology, comparing numerical predictions with experimental data from Torbergsen (1998) and analytical correlations from classical turbulent flow theory. Special thanks to the fluid mechanics research community for providing high-quality experimental benchmark data.
