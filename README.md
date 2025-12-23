# BARRIER - Building Analysis capacity to improve the Resilience of the Guadalquivir river estuary

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Version](https://img.shields.io/badge/Version-0.6.0-blue.svg)](https://github.com/mcobosb/Barrier)
[![C++](https://img.shields.io/badge/C%2B%2B-20-blue.svg)](https://isocpp.org/)
[![CMake](https://img.shields.io/badge/CMake-3.16%2B-green.svg)](https://cmake.org/)
[![OpenMP](https://img.shields.io/badge/OpenMP-parallel-orange.svg)](https://www.openmp.org/)

## 📖 Overview

BARRIER is a high-performance 1D estuarine hydrodynamic model designed to simulate water flow, sediment transport, and salinity dynamics in the Guadalquivir river estuary. The project aims to **Build Analysis capacity to improve the Resilience** of this critical estuarine system through advanced computational modeling.

### Key Features

- **TVD McCormack Scheme**: Second-order accurate predictor-corrector method with flux limiters for stability
- **Saint-Venant Equations**: Solves 1D shallow water equations with source terms (friction, bed slope, density gradients)
- **Estuarine Processes**: Models tidal dynamics, salinity intrusion, sediment transport, and lateral inflows
- **Boundary Conditions**: Flexible upstream/downstream conditions (open flow, discharge, or water elevation)
- **Dry Bed Handling**: Robust treatment of intertidal zones and flood/dry transitions
- **High Performance**: OpenMP parallelization with optimized build modes (Release, ReleaseFast, PGO)
- **NetCDF Output**: CF-compliant scientific data format for analysis and visualization
- **Cross-platform**: Runs on Linux, macOS, and Windows systems

---

## 🏗️ Project Structure

```
Barrier/
├── src/                    # Source code implementation
│   ├── main.cpp            # Entry point and exception handling
│   ├── simulation.cpp      # Core TVD McCormack solver and time loop
│   ├── data_reader.cpp     # Configuration and input file parsing
│   ├── data_writer.cpp     # NetCDF output generation
│   ├── cross_section.cpp   # Geometry and hydraulic properties
│   ├── hydrograph.cpp      # Lateral inflow (tributary) handling
│   ├── screen_presenter.cpp # Console output and progress display
│   ├── utils.cpp           # Utility functions
│   └── error_handling.cpp  # Error management and logging
├── include/                # Header files
│   ├── simulation.h        # CSimulation class definition
│   ├── data_reader.h       # Input parsing interfaces
│   ├── data_writer.h       # Output interfaces
│   ├── cross_section.h     # CCrossSection class
│   ├── hydrograph.h        # CHydrograph class
│   ├── screen_presenter.h  # Console output
│   ├── utils.h             # Utility functions
│   ├── error_handling.h    # Error handling
│   └── main.h              # Global constants (G, DRY_AREA, etc.)
├── tests/                  # Test cases with input data
│   ├── 1022/               # Example: 251 sections, 7-day tidal simulation
│   │   ├── .conf           # Configuration file
│   │   ├── cross_sections.csv      # Cross-section geometry
│   │   ├── along_channel_data.csv  # Manning coefficients, bed elevation
│   │   ├── tides.csv       # Downstream boundary condition
│   │   └── *.nc            # NetCDF output files
│   └── [other test cases]
├── viewer/                 # Python visualization tools
│   ├── app.py              # Visualization application
│   └── environment.yml     # Python dependencies
├── docs/                   # Documentation
├── build/                  # Build directory (created by CMake)
├── CMakeLists.txt          # Build configuration
├── CHANGELOG.md            # Version history
├── LICENSE                 # GPL v3.0 license
└── README.md               # This file
```

---

## ⚙️ Configuration

The simulation is now configured through a **YAML file** (`config.yaml`). Each test case directory (e.g. `tests/1022/`) contains its own `config.yaml` file.

### Example YAML Configuration (`config.yaml`)

```yaml
run:
   name: "1022"
   start_date: "1999-01-01T00:00:00"
   duration: 604800  # seconds (7 days)
   timestep: 3600  # seconds (1 hour)
   output_variables: full
   log_level: 0  # no log file
   continue_simulation: true
   continue_netcdf_path: "1022_barrier_sim_20251215_1034_CS251_T7d_dt3600_BC12_SAL_TVD_DRY_CFL15.nc"

geometry:
   along_channel_file: along_channel_data
   cross_sections_file: cross_sections

initial_conditions:
   type: 0  # in calm

boundary_conditions:
   upstream:
      type: 1  # reflective
      file: ""
   downstream:
      type: 2  # elevation
      file: tides

forcing:
   tributaries_file: ""  # no tributaries

numerics:
   courant_number: 0.15
   tvd_limiter:
      enabled: true
      method: minmod  # 1 = MinMod
      psi_formula: tseng  # 2 = Tseng
      delta: 0.10
   surface_gradient: true
   source_term_balance: true
   beta_coefficient: true
   dry_bed: true
   murillo_condition: true

transport:
   salinity:
      enabled: true
      initial_file: ""
      upstream_condition: 0  # free
      downstream_condition: 2  # ocean
      beta: 0.00076
      dispersion_kh: 750.0  # m²/s
   sediment:
      enabled: false
      equation: 0  # van Rijn
      properties_file: ""
   density:
      enabled: false

smoothing:
   bathymetry: true
   solution: true
```

#### Main YAML Sections

| Section | Description |
|---------|-------------|
| `run` | General run settings: name, start date, duration, timestep, output, restart options |
| `geometry` | Input files for along-channel and cross-sections |
| `initial_conditions` | Type of initial state (0=in calm, 1=flow, 2=elevation) |
| `boundary_conditions` | Upstream and downstream boundary types and files |
| `forcing` | Tributary/river inflow files |
| `numerics` | Courant number, TVD limiter, numerical flags |
| `transport` | Salinity, sediment, and density transport options |
| `smoothing` | Bathymetry and solution smoothing flags |

See the example above for typical parameter names and values. All input files referenced (e.g. `along_channel_data`, `cross_sections`, `tides`) should be present in the same directory as the YAML file.

> **Note:** The old `.ini`, `.conf`, and `.config` formats are no longer supported. Use only `config.yaml` for new simulations.

## 📂 Input Data Files

All input files must be in **CSV format** with semicolon (`;`) delimiters and placed in the same directory as the `.conf` file.

### Cross-Section Geometry (`cross_sections.csv`)

Defines the cross-section geometry at multiple water elevations. The model interpolates to obtain area (A), wetted perimeter (P), and hydraulic radius (Rh) for any water elevation.

**Format:**
```csv
;x,z,eta,B,A,P,Rh,sigma,xl,xr,beta
0.0,-12.0,0.0,367.879,0.0,735.759,0.0,367.879,-183.94,183.94,1
0.0,-12.0,0.625,367.879,229.925,735.857,0.312,367.879,-183.94,183.94,1
0.0,-12.0,1.25,367.879,459.849,736.146,0.625,367.879,-183.94,183.94,1
...
```

**Columns:**
- `x`: Longitudinal coordinate along estuary (m)
- `z`: Bed elevation (m, positive upward)
- `eta`: Water elevation above bed (m)
- `B`: Top width (m)
- `A`: Cross-sectional area (m²)
- `P`: Wetted perimeter (m)
- `Rh`: Hydraulic radius (m)
- `sigma`: Storage width (m, for floodplain storage)
- `xl`, `xr`: Left/right bank positions (m)
- `beta`: Momentum correction coefficient

### Along-Channel Data (`along_channel_data.csv`)

Defines bed elevation and Manning's roughness coefficient at each computational section.

**Format:**
```csv
;x,z,nmann,xutm,yutm,AngMd,AngMi
0.0,-12.0,0.03,0,0,0,0
400.0,-12.0,0.03,0,0,0,0
800.0,-12.0,0.03,0,0,0,0
...
```

**Columns:**
- `x`: Longitudinal coordinate (m)
- `z`: Bed elevation (m)
- `nmann`: Manning's roughness coefficient (s/m^(1/3))
- `xutm`, `yutm`: UTM coordinates (optional, for georeferencing)
- `AngMd`, `AngMi`: Angle degrees/minutes (optional)

### Boundary Conditions

#### Tidal Boundary (`tides.csv`)
Time series of water elevation at the downstream boundary.

**Format:**
```csv
;,eta
0.0,0.5
605.4,0.498
1210.8,0.493
...
```

**Columns:**
- First column: Time (seconds from simulation start)
- `eta`: Water elevation (m)

#### Discharge Boundary (`hydro.csv`)
Time series of discharge at the upstream boundary or lateral inflows.

**Format:**
```csv
;t,Q
0.0,150.0
3600.0,180.0
7200.0,200.0
...
```

**Columns:**
- `t`: Time (seconds)
- `Q`: Discharge (m³/s)

---

## 🧩 Numerical Method

The model implements a **TVD (Total Variation Diminishing) McCormack scheme**, a second-order accurate predictor-corrector method with flux limiters for shock-capturing and monotonicity preservation.

### Main Simulation Loop

The core algorithm in `CSimulation::bDoSimulation()` follows these steps:

1. **Read Configuration** (`CDataReader::readConfigurationFile()`)
   - Parse `config.yaml` file
   - Read cross-section geometry and along-channel data
   - Load boundary conditions and initial conditions

2. **Initialize Simulation** 
   - Precompute geometric properties (`precomputeEstuaryData()`)
   - Calculate bed slopes (`calculateBedSlope()`)
   - Set initial conditions (`calculateAlongEstuaryInitialConditions()`)

3. **Main Time Loop** (until `m_dCurrentTime` reaches `m_dSimDuration`)
   ```
   while (t < T_final):
       a) Check Courant condition (adaptive time step)
       b) Update boundary conditions (interpolate from time series)
       c) PREDICTOR STEP (forward differences):
          - Calculate predictor values using forward finite differences
          - Apply source terms (friction, bed slope, density gradients)
       d) Update predictor boundaries
       e) CORRECTOR STEP (backward differences):
          - Calculate corrector values using backward finite differences
          - Apply source terms
       f) Update corrector boundaries
       g) MERGE PREDICTOR-CORRECTOR with TVD flux limiters:
          - Apply MinMod/Roe/Van Leer/Van Albada limiter
          - Ensure monotonicity and prevent spurious oscillations
       h) Optional: Smooth solution at boundaries (sponge layer)
       i) Apply dry bed treatment if enabled
       j) Update salinity and sediment transport if enabled
       k) Write output to NetCDF at specified intervals
   ```

4. **Write Final Output** (`CDataWriter::writeNetCDF()`)
   - Save simulation results to NetCDF file
   - Include metadata (simulation parameters, coordinates, time)

### Key Methods

| Method | Description |
|--------|-------------|
| `calculateBedSlope()` | Computes bed slope using central/forward/backward differences |
| `calculateBoundaryConditions()` | Interpolates boundary values (elevation → area, discharge) |
| `calculatePredictor()` | Predictor step: forward finite differences |
| `calculateCorrector()` | Corrector step: backward finite differences |
| `mergePredictorCorrector()` | Merges predictor/corrector with TVD flux limiters |
| `calculate_GS_A_terms()` | Bed slope source terms (gravity × area × slope) |
| `calculateSourceTerms()` | Friction source terms with adaptive upwind scheme |
| `smoothSolution()` | Selective smoothing at boundaries (sponge layer) |
| `calculate_sediment_transport()` | Sediment transport equations (if enabled) |
| `calculate_salinity()` | Salinity advection-diffusion (if enabled) |

### Boundary Conditions

- **Type 0 (Open)**: Zero gradient (∂Q/∂x = 0, ∂A/∂x = 0)
- **Type 1 (Discharge)**: Prescribed discharge Q(t) from time series
- **Type 2 (Elevation)**: Prescribed water elevation η(t), interpolated to area A(t)

### TVD Flux Limiters

Available limiters for `mergePredictorCorrector()`:
1. **MinMod**: Most diffusive, stable for strong shocks
2. **Roe**: Balanced diffusion/dispersion
3. **Van Leer**: Smooth, good for continuous flows
4. **Van Albada**: Similar to Van Leer, alternative formulation

---

## 🚀 Quick Start

### Prerequisites

**Ubuntu/Debian:**
```bash
sudo apt-get update
sudo apt-get install build-essential cmake libnetcdf-dev libnetcdf-c++4-dev libomp-dev
```

**CentOS/RHEL:**
```bash
sudo yum install gcc-c++ cmake netcdf-devel openmp-devel
```

**macOS:**
```bash
brew install cmake netcdf libomp
```

### Building the Project

```bash
# Clone the repository
git clone https://github.com/mcobosb/Barrier.git
cd Barrier

# Create and enter build directory
mkdir -p build && cd build

# Configure with CMake (choose build type)
cmake -DCMAKE_BUILD_TYPE=ReleaseFast ..

# Build with parallel compilation
make -j$(nproc)

# The executable 'barrier' will be created in the build directory
```

### Build Types

BARRIER supports multiple build configurations optimized for different purposes:

| Build Type | Optimization | Use Case |
|------------|-------------|----------|
| `Debug` | `-O0 -g` | Development, debugging with gdb, sanitizers |
| `Release` | `-O3 -march=native -flto` | Production use, maximum performance |
| `ReleaseFast` | `-O3 -march=native -ffast-math` | **Maximum speed** (default) |
| `RelWithDebInfo` | `-O3 -g` | Profiling with performance tools |
| `PGOGenerate` | Profile-guided optimization (step 1) | Generate profiling data |
| `PGOUse` | Profile-guided optimization (step 2) | Use profiling data for optimization |

**Example - Debug build:**
```bash
cmake -DCMAKE_BUILD_TYPE=Debug ..
make -j$(nproc)
```

### Running a Simulation

```bash
# Navigate to a test case directory
cd ../tests/1022

# Run the simulation
../../build/barrier .conf

# Or run from the build directory with full path
cd ../../build
./barrier ../tests/1022/.conf
```

**Output:**
- NetCDF file: `barrier_sim_YYYYMMDD_HHMM_[parameters].nc`
- Console output: Progress bar, timing information, diagnostics
- Log file: (if enabled in configuration)

### Example Output

```
Hydrodynamic modeling of main channel of estuaries (03 December 2025) - BARRiEr Project (https://barrier.dinamicambiental.es) for Intel x86/GNU C++ (12:22:12 Dec 15 2025 build)
(C) 2025 Manuel Cobos

-------------------------------------------------------------------------------
This program is distributed in the hope that it will be useful. but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details. You
should have received a copy of the GNU General Public License along with this
program; if not. contact the Free Software Foundation. Inc.. 675 Mass Ave.
Cambridge. MA 02139. USA.
-------------------------------------------------------------------------------

    - Started on mcobos at 12:22 on Monday 15 December 2025
    - Initializing
      - Detected configuration file
      - Loading initial NetCDF for continuing simulation: testGRE_barrier_sim_20251215_1040_CS806_T7d_dt3600_BC12_SAL_TVD_DRY_CFL15.nc
      - Smoothing bathymetry (3 passes)
      - Output file: testGRE_barrier_sim_20251215_1222_CS806_T7d_dt3600_BC12_SAL_TVD_DRY_CFL15.nc
    - Running
    - Elapsed[Remaining] Time: 72.000[ 0.000] s - Progress: 100.000%
    - Run ended at 12:23 on Monday 15 December 2025

Run successfully end
```

## 📈 Output

Results are saved in **NetCDF format** (Climate and Forecast conventions) in the same directory as the configuration file.

### Output File Naming Convention

```
barrier_sim_YYYYMMDD_HHMM_CS###_T#d_dt####_BC##_[options].nc
```

**Example:**
```
barrier_sim_20251203_0330_CS251_T7d_dt3600_BC02_SAL_TVD_DRY_CFL10.nc
```

**Filename components:**
- `YYYYMMDD_HHMM`: Simulation start date and time
- `CS###`: Number of cross-sections (e.g., CS251 = 251 sections)
- `T#d`: Simulation duration (e.g., T7d = 7 days)
- `dt####`: Time step in seconds (e.g., dt3600 = 1 hour)
- `BC##`: Boundary condition type (e.g., BC02 = downstream elevation)
- `SAL`: Salinity transport enabled
- `TVD`: TVD flux limiter enabled
- `DRY`: Dry bed handling enabled
- `CFL##`: Courant number × 100 (e.g., CFL10 = 0.1)

### Output Variables

The NetCDF file contains the following variables:

| Variable | Units | Description | Dimensions |
|----------|-------|-------------|------------|
| `time` | seconds | Time since simulation start | `(time)` |
| `x` | meters | Longitudinal coordinate | `(x)` |
| `z` | meters | Bed elevation | `(x)` |
| `eta` | meters | Water surface elevation | `(time, x)` |
| `h` | meters | Water depth | `(time, x)` |
| `A` | m² | Cross-sectional area | `(time, x)` |
| `Q` | m³/s | Discharge | `(time, x)` |
| `U` | m/s | Mean velocity (Q/A) | `(time, x)` |
| `Fr` | - | Froude number | `(time, x)` |
| `S` | psu | Salinity | `(time, x)` |
| `rho` | kg/m³ | Water density | `(time, x)` |
| `C` | kg/m³ | Sediment concentration | `(time, x)` |

*Note: Variables depend on enabled options in configuration (salinity, sediment transport, etc.)*

### Metadata

The NetCDF file includes comprehensive metadata:
- Simulation parameters (time step, Courant number, etc.)
- Boundary conditions
- Numerical scheme settings
- Creation date and model version
- CF-1.8 compliance for interoperability

### Visualization

Results can be visualized using:

- **ncview**: Quick NetCDF file viewer
  ```bash
  ncview barrier_sim_*.nc
  ```
- **Panoply**: NASA's NetCDF viewer (GUI)
- **ParaView**: Advanced 3D visualization

---

## 🧪 Mathematical Foundation

The model solves the **1D Saint-Venant equations** (shallow water equations) with source terms for estuarine flow:

### Governing Equations

**Continuity Equation:**

$$\frac{\partial A}{\partial t} + \frac{\partial Q}{\partial x} = q_l$$

**Momentum Equation:**

$$\frac{\partial Q}{\partial t} + \frac{\partial}{\partial x}\left(\beta\frac{Q^2}{A}\right) + gA\frac{\partial \eta}{\partial x} = -gAS_f - gAS_0 + S_\rho + M_l$$

Where:
- $A$ = Cross-sectional area (m²)
- $Q$ = Discharge (m³/s)
- $\eta$ = Water surface elevation (m)
- $h$ = Water depth (m)
- $q_l$ = Lateral inflow per unit length (m²/s)
- $M_l$ = Momentum added by lateral inflow (m³/s²)
- $g$ = Gravitational acceleration (9.81 m/s²)
- $\beta$ = Momentum correction coefficient
- $S_f$ = Friction slope
- $S_0$ = Bed slope ($-\partial z_b / \partial x$)
- $S_\rho$ = Density gradient source term (for salinity)

### Friction Slope

The friction slope is calculated using Manning's equation:

$$S_f = \frac{n^2 |Q| Q}{A^2 R_h^{4/3}}$$

Where:
- $n$ = Manning's roughness coefficient (s/m^(1/3))
- $R_h$ = Hydraulic radius ($A/P$, m)
- $P$ = Wetted perimeter (m)

### Density Gradient Term

For salinity-driven circulation:

$$S_\rho = -gA\frac{1}{\rho}\frac{\partial \rho}{\partial x}$$

Where $\rho = \rho(S)$ is the water density as a function of salinity.

### Numerical Scheme

The model uses the **TVD McCormack scheme**, a second-order predictor-corrector method:

1. **Predictor (forward differences):**
   $$U_i^* = U_i^n - \frac{\Delta t}{\Delta x}(F_{i+1}^n - F_i^n) + \Delta t \cdot S_i^n$$

2. **Corrector (backward differences):**
   $$U_i^{**} = U_i^* - \frac{\Delta t}{\Delta x}(F_i^* - F_{i-1}^*) + \Delta t \cdot S_i^*$$

3. **TVD limiter:**
   $$U_i^{n+1} = U_i^* + \frac{1}{2}\psi(r_i)(U_i^{**} - U_i^*)$$

Where $\psi(r)$ is the flux limiter function (MinMod, Van Leer, etc.) and $r_i$ is the smoothness indicator.

### Stability Condition

The Courant-Friedrichs-Lewy (CFL) condition:

$$\Delta t \leq C \frac{\Delta x}{\max(|U| + \sqrt{gA/B})}$$

Where $C$ is the Courant number (typically 0.1 - 0.5).

---

## 🔧 Development

### Building for Development

**Debug build with sanitizers:**
```bash
cmake -DCMAKE_BUILD_TYPE=Debug ..
make -j$(nproc)

# Run with address sanitizer and undefined behavior sanitizer
./barrier ../tests/1022/.conf
```

**Profile-Guided Optimization (PGO):**
```bash
# Step 1: Build with profiling instrumentation
cmake -DCMAKE_BUILD_TYPE=PGOGenerate ..
make -j$(nproc)

# Step 2: Run representative test cases to generate profiling data
./barrier ../tests/1022/.conf
./barrier ../tests/0004/.conf

# Step 3: Rebuild with profiling data
cd ..
rm -rf build && mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=PGOUse ..
make -j$(nproc)

# The resulting binary is optimized for your typical workload
```

### Performance Profiling

**Using `perf` (Linux):**
```bash
# Record performance data
perf record -g ./barrier ../tests/1022/.conf

# Analyze hotspots
perf report
```

**Using `gprof`:**
```bash
# Build with profiling
cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_CXX_FLAGS="-pg" ..
make -j$(nproc)

# Run simulation
./barrier ../tests/1022/.conf

# Analyze profile
gprof ./barrier gmon.out > analysis.txt
```

### Code Style

The project follows modern C++20 practices:
- **STL containers and algorithms**: Prefer `std::vector`, `std::array`, range-based for loops
- **RAII principles**: Automatic resource management
- **Const-correctness**: Use `const` wherever possible
- **Smart pointers**: For dynamic memory (if needed)
- **OpenMP parallelization**: Use `#pragma omp parallel for` for data-parallel loops
- **Meaningful names**: Descriptive variable/function names (e.g., `calculateBedSlope()`)

### Testing

```bash
# Run all test cases
cd tests
for dir in */; do
    echo "Testing $dir"
    cd "$dir"
    ../../build/barrier .conf
    cd ..
done
```

### Adding New Features

1. **New source term**: Implement in `calculateSourceTerms()` in `simulation.cpp`
2. **New boundary condition type**: Extend `calculateBoundaryConditions()`
3. **New output variable**: Add to `CDataWriter::writeNetCDF()`
4. **New input file format**: Extend `CDataReader::readConfigurationFile()`

---

## 📝 Documentation

### Code Documentation

- **Source files**: Comprehensive inline comments explaining algorithms
- **Header files**: Doxygen-style documentation for all classes and methods
- **CHANGELOG.md**: Detailed version history and changes

### References

- **TVD McCormack Scheme**: Toro, E.F. (2009). *Riemann Solvers and Numerical Methods for Fluid Dynamics*. Springer.
- **Saint-Venant Equations**: Chow, V.T. (1959). *Open-Channel Hydraulics*. McGraw-Hill.
- **Estuarine Dynamics**: Dyer, K.R. (1997). *Estuaries: A Physical Introduction*. Wiley.
- **NetCDF CF Conventions**: [Climate and Forecast Conventions](http://cfconventions.org/)

---

## 🤝 Contributing

Contributions are welcome! Please follow these guidelines:

1. **Fork the repository**
   ```bash
   git clone https://github.com/mcobosb/Barrier.git
   ```

2. **Create a feature branch**
   ```bash
   git checkout -b feature/amazing-feature
   ```

3. **Make your changes**
   - Follow the existing code style
   - Add comments for complex algorithms
   - Test your changes with existing test cases

4. **Commit your changes**
   ```bash
   git commit -m 'Add amazing feature: description'
   ```

5. **Push to your fork**
   ```bash
   git push origin feature/amazing-feature
   ```

6. **Open a Pull Request**
   - Describe your changes
   - Reference any related issues
   - Include test results

### Reporting Issues

If you find a bug or have a feature request:
- Check existing issues first
- Provide a minimal reproducible example
- Include system information (OS, compiler version, dependencies)
- Attach relevant error messages or logs

---

## 📄 License

This project is licensed under the **GNU General Public License v3.0** - see the [LICENSE](LICENSE) file for details.

**Key points:**
- ✅ You can use this software for any purpose
- ✅ You can modify and distribute modified versions
- ⚠️ Modified versions must also be licensed under GPL v3.0
- ⚠️ You must include the original copyright and license notices
- ❌ No warranty is provided

---

## 🙏 Acknowledgments

### BARRIER Project

**Building Analysis capacity to improve the Resilience of the Guadalquivir river estuary**

- **Project website**: [https://barrier.dinamicambiental.es](https://barrier.dinamicambiental.es)
- **Research Group**: [Environmental Fluid Dynamics Group](https://dinamicambiental.es)
- **Institution**: University of Granada, Spain

### Research Team

- **Dr. Manuel Cobos Budia** - Principal Investigator & Software Developer
- **Dr. Asunción Baquerizo Azofra** - Co-Investigator
- **Dr. Manuel Díez Minguito** - Co-Investigator
- **Dr. Agustín Millares Valenzuela** - Co-Investigator

### Funding

This work has been carried out within the framework of the project **BARRIER**:

- **Reference**: PID2023-148298OA-I00
- **Program**: Proyectos de Generación de Conocimiento - Call 2023
- **Agency**: Ministerio de Ciencia, Innovación y Universidades (Spain)

### Dependencies

This project uses the following open-source libraries:
- **NetCDF**: [Unidata NetCDF](https://www.unidata.ucar.edu/software/netcdf/)
- **OpenMP**: [OpenMP Architecture Review Board](https://www.openmp.org/)
- **CMake**: [Kitware CMake](https://cmake.org/)

---

## 📧 Contact

**Project Maintainer**: Manuel Cobos Budia  
**Email**: [mcobosb@ugr.es](mailto:mcobosb@ugr.es)  
**Institution**: Environmental Fluid Dynamics Group, University of Granada  
**Website**: [https://dinamicambiental.es](https://dinamicambiental.es)

### Support

- **Technical Issues**: Open an issue on [GitHub](https://github.com/mcobosb/Barrier/issues)
- **Scientific Collaboration**: Contact the research team via email
- **General Questions**: mcobosb@ugr.es

---

## 📊 Citation

If you use BARRIER in your research, please cite:

```bibtex
@software{barrier2026,
  title = {{BARRIER}: 1D Estuarine Hydrodynamic Model},
  author = {Cobos, Manuel},
  year = {2026},
  version = {0.11.0},
  url = {https://github.com/mcobosb/Barrier},
  note = {Building Analysis capacity to improve the Resilience of the Guadalquivir river estuary}
}
```

---

*For technical support, scientific collaboration, or general inquiries, please contact the development team.*

**Last updated**: December 2025 | **Version**: 0.11.0
