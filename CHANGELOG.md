# Changelog

## [0.5.0] - 2025-12-03

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