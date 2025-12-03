# Changelog

## [0.6.0] - 2025-12-03

### Major Features
- **YAML Configuration Support**: Complete migration from legacy `.conf` format to modern YAML
  - New `CYAMLReader` class for parsing YAML configuration files
  - ISO 8601 date format support: `"YYYY-MM-DDTHH:MM:SS"`
  - Hierarchical configuration structure with clear sections (run, geometry, boundary_conditions, numerics, transport, smoothing)
  - Automatic detection by file extension (`.yaml`/`.yml` vs `.conf`)
  - Full backward compatibility: both formats supported simultaneously
  - All test cases converted to YAML: `1022`, `reducedGRE`, `qGRE`, `reducedGRE_z0`, `x2`, `calibration`

- **UTM Coordinate Calculation**: Geographic positioning of river banks
  - New method `calculateRiverBankUTMCoordinates()` computes absolute UTM coordinates
  - Uses thalweg position, bank angles (radians), and lateral distances
  - Formula: `X_utm = X_thalweg + distance * cos(angle)`, `Y_utm = Y_thalweg + distance * sin(angle)`
  - Four new output variables: `xl_utm_x`, `xl_utm_y`, `xr_utm_x`, `xr_utm_y`
  - Integrated into NetCDF output system with proper metadata

### Critical Bug Fixes
- **Manning Friction Coefficient Bug (CRITICAL)**: Fixed zero friction throughout all simulations
  - **Issue**: `m_vManningNumberSquared` was calculated in `precomputeEstuaryData()` BEFORE Manning values were read from files, resulting in all zeros
  - **Impact**: All friction forces were zero, making Manning coefficient adjustments have NO effect on simulation results
  - **Fix**: Moved calculation to `calculateAlongEstuaryInitialConditions()` AFTER Manning values are loaded from `estuary[]` objects
  - **Result**: Friction now works correctly, tidal wave propagation speed now controllable via Manning coefficient

- **River Bank Position Initialization**: Fixed zero values at t=0
  - Added call to `calculateHydraulicParameters()` before first timestep output
  - Ensures `m_vCrossSectionLeftRBLocation` and `m_vCrossSectionRightRBLocation` are initialized
  - Bank positions now correctly interpolated in case 3 (interpolation branch) of `calculateHydraulicParameters()`

### New Test Case
- **calibration**: Real GRE estuary data for 2007-2011 period
  - Duration: 1,720 days (148,608,000 seconds)
  - Start date: 09/02/2007
  - Manning coefficient: 0.035 (tuned for tidal wave propagation)
  - Time-aligned forcing data: `tides.csv` (247,680 points), `hydro.csv` (1,721 daily values)
  - Time offset: 2,386,108,800 seconds from 1931-07-01 base date
  - YAML configuration with detailed comments

### Code Quality
- **Compiler Warning Cleanup**: Removed all warnings
  - Eliminated unused variable `dx` in `calculateSourceTerms()`
  - Eliminated unused variable `nCrossSectionsNumber` in `mergePredictorCorrector()`
  - Fixed member variable name mismatches in YAML reader integration

- **Dependencies**: Added yaml-cpp 0.8.0 via pkg-config
  - Integrated into CMake build system
  - Proper include directories and library linking

### Removed Features
- **SPONGE LAYER**: Removed selective boundary smoothing feature
  - Feature was documented in v0.5.0 but has been eliminated from codebase
  - Smoothing now controlled globally via configuration settings

### Technical Improvements
- Enhanced `simulation.h`: Added 4 new member variables for UTM coordinates (lines ~305-318)
- Enhanced `data_writer.cpp`: Added metadata definitions for UTM variables (lines ~124-147)
- Enhanced `data_reader.cpp`: Updated "full" variables list to include UTM coordinates (line ~362)
- New files: `include/yaml_reader.h` (66 lines), `src/yaml_reader.cpp` (323 lines)
- Updated `CMakeLists.txt`: Added yaml_reader.cpp to sources and yaml-cpp dependency

### Quality Assurance
- ✅ Code compiles without errors or warnings
- ✅ Both `.conf` and `.yaml` configuration formats validated
- ✅ Manning friction now correctly applied (major physics bug fixed)
- ✅ UTM coordinates calculation verified
- ✅ All test cases ready with YAML configurations
- ✅ Backward compatibility maintained

### Summary
This release includes a critical bug fix for Manning friction (was completely broken - all friction was zero), implements modern YAML configuration support, adds geographic UTM coordinate calculation for river banks, and creates a new calibration test case with real 2007-2011 data. The Manning bug fix is especially important as it means all previous simulations had incorrect physics. The codebase is now ready for proper model calibration with correctly functioning friction.

## [0.5.0] - 2025-12-02

### Documentation Improvements
- **Enhanced file-level documentation** across all major source files:
  - `data_writer.cpp`: Added detailed description of NetCDF file handling and output management
  - `data_reader.cpp`: Documented configuration parsing and input data validation
  - `main.cpp`: Clarified program entry point and exception handling flow
  - `cross_section.cpp`: Explained geometry storage and hydraulic parameter computation
  - `hydrograph.cpp`: Documented lateral inflow (tributary) handling and interpolation
  - `screen_presenter.cpp`: Described console output and progress reporting
- **Updated header file documentation** for `hydrograph.h`, `main.h`, `error_handling.h`, and `screen_presenter.h`

### Code Cleanup
- **Removed obsolete commented code** to improve readability:
  - Cleaned up large commented blocks in `updatePredictorBoundaries()` and `updateBoundaryConditions()`
  - Removed legacy code fragments from main simulation loop
  - Eliminated references to deprecated variables (e.g., `m_vCrossSectionBedSlopeDirection`)
- **Updated or removed TODO comments**:
  - TODO 003 → "Courant number accounting for density effects" (clarified intent)
  - TODO 014 → Removed (redundant with existing comments)
  - TODO 020/007 → "Time unit conversion (currently only seconds supported)" (documented limitation)
  - TODO 010 → "Store boundary condition value" (clarified purpose)
  - TODO 001 → Replaced generic placeholders with specific descriptions

### Technical Documentation
- **Added comprehensive function-level comments**:
  - `calculateBedSlope()`: Explained central, forward, and backward difference schemes for bed slope
  - `calculatePredictor()`/`calculateCorrector()`: Documented McCormack predictor-corrector method
  - `calculateBoundaryConditions()`: Described elevation-to-area interpolation for tidal boundaries
  - `calculate_GS_A_terms()`: Explained bed slope source terms and friction calculations
  - `calculateSourceTerms()`: Documented adaptive upwind scheme for steep bed slope changes
  - `mergePredictorCorrector()`: Described TVD flux limiters and monotonicity preservation

### Comment Quality Improvements
- **Enhanced key algorithm descriptions**:
  - Main time loop: "Main time-stepping loop: advance simulation until final time"
  - Boundary conditions: Clarified types (0=open, 1=reflective/flow, 2=elevation/tide)
  - SPONGE LAYER: Explained selective smoothing at boundaries to reduce discontinuities
  - Upwind scheme: Documented detection threshold for abrupt bed slope changes (>0.1%)
  - I1 term: Added note explaining why non-prismatic channel term is disabled (numerical instability)
- **Improved inline comments**:
  - Predictor-corrector steps clearly labeled as "forward differences" / "backward differences"
  - Flow direction logic: "Positive flow (downstream): use backward difference"
  - Boundary handling: "Boundaries handled separately" instead of vague references

### Quality Assurance
- ✅ Code compiles successfully without errors or warnings
- ✅ Test case 1022 validated (7-day simulation completes in ~4 seconds)
- ✅ No functional regressions introduced
- ✅ All changes preserve original numerical behavior

### Summary
This release focuses on code quality and maintainability improvements. The codebase is now significantly cleaner, better documented, and easier to understand for both current maintenance and future development. All changes are non-functional and preserve the original computational behavior of the model.