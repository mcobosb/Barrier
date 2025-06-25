# BARRIER - Building Analysis capacity to improve the Resilience of the Guadalquivir river estuary

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![C++](https://img.shields.io/badge/C%2B%2B-20-blue.svg)](https://isocpp.org/)
[![CMake](https://img.shields.io/badge/CMake-3.16%2B-green.svg)](https://cmake.org/)
[![OpenMP](https://img.shields.io/badge/OpenMP-parallel-orange.svg)](https://www.openmp.org/)

## 📖 Overview

BARRIER is a high-performance numerical simulation software designed to model hydrodynamic processes in the Guadalquivir river estuary. The project aims to **Build Analysis capacity to improve the Resilience** of this critical estuarine system through advanced computational modeling.

### Key Features

- **Saint-Venant Equations**: Solves 1D shallow water equations for river flow modeling
- **Hydrodynamic Simulation**: Models water flow, sediment transport, and water quality parameters
- **High Performance**: Parallel processing with OpenMP for efficient computation
- **NetCDF Output**: Standard scientific data format for results visualization and analysis
- **Cross-platform**: Runs on Linux, macOS, and Windows systems

## 🏗️ Project Structure

```
Barrier/
├── src/                    # Source code
│   ├── main.cpp           # Main application entry point
│   ├── simulation.cpp     # Core simulation engine
│   ├── data_reader.cpp    # Input data processing
│   ├── data_writer.cpp    # NetCDF output generation
│   ├── cross_section.cpp  # River cross-section modeling
│   └── hydrograph.cpp     # Hydrographic data handling
├── include/               # Header files
├── input/                 # Input data files
├── output/                # Simulation results
├── build/                 # Build directory
└── CMakeLists.txt         # Build configuration
```

## 🚀 Quick Start

### Prerequisites

```bash
# Ubuntu/Debian
sudo apt-get update
sudo apt-get install build-essential cmake libnetcdf-dev libomp-dev

# CentOS/RHEL
sudo yum install gcc-c++ cmake netcdf-devel openmp-devel

# macOS
brew install cmake netcdf libomp
```

### Building the Project

```bash
# Clone the repository
git clone https://github.com/your-username/Barrier.git
cd Barrier

# Create and enter build directory
mkdir build && cd build

# Configure and build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$(nproc)
```

### Running a Simulation

```bash
# Run the simulation
./barrier

# The program will read configuration from '.ini'
# and generate results in the specified output directory
```

## ⚙️ Configuration

The simulation is configured through the `.ini` file:

```ini
# Configuration file for Barrier simulation
input_path = "/path/to/input/data"
output_path = "/path/to/output/results"
simulation_duration = 86400.0  # seconds
time_step = 1.0               # seconds
cross_sections_file = "cross_sections.dat"
boundary_conditions = "boundaries.dat"
```

## 📊 Input Data Format

### Cross-Sections Data
River geometry defined by cross-sectional profiles:
```
# X-coordinate, Z-elevation, Width, Manning coefficient
0.0  -5.0  100.0  0.035
1000.0  -4.8  95.0   0.032
2000.0  -4.5  90.0   0.030
```

### Boundary Conditions
Time series of water levels or discharge:
```
# Time(s), Discharge(m³/s) or Water_Level(m)
0.0     150.0
3600.0  180.0
7200.0  200.0
```

## 📈 Output

Results are saved in **NetCDF format** containing:

- **Water elevation** (m)
- **Flow velocity** (m/s)
- **Discharge** (m³/s)
- **Cross-sectional area** (m²)
- **Hydraulic radius** (m)
- **Sediment concentration** (kg/m³) - if enabled

### Visualization

Results can be visualized using:
- **Python**: `xarray`, `matplotlib`, `pandas`
- **MATLAB**: Built-in NetCDF support
- **ParaView**: For advanced 3D visualization
- **ncview**: Quick NetCDF file viewer

Example Python visualization:
```python
import xarray as xr
import matplotlib.pyplot as plt

# Load simulation results
ds = xr.open_dataset('output/simulation_results.nc')

# Plot water elevation over time
ds.water_elevation.plot(x='time', y='x')
plt.title('Water Elevation Along Estuary')
plt.show()
```

## 🧪 Mathematical Foundation

The model solves the Saint-Venant equations for 1D unsteady flow:

**Continuity Equation:**
![Continuity](https://latex.codecogs.com/svg.latex?\frac{\partial%20A}{\partial%20t}%20+%20\frac{\partial%20Q}{\partial%20x}%20=%200)

**Momentum Equation:**
![Momentum](https://latex.codecogs.com/svg.latex?\frac{\partial%20Q}{\partial%20t}%20+%20\frac{\partial}{\partial%20x}\left(\frac{Q^2}{A}\right)%20+%20gA\frac{\partial%20h}{\partial%20x}%20+%20gAS_f%20=%200)

Where:
- A = Cross-sectional area (m²)
- Q = Discharge (m³/s)
- h = Water depth (m)
- Sf = Friction slope
- g = Gravitational acceleration (9.81 m/s²)

## 🔧 Development

### Building for Development

```bash
# Debug build with additional checks
cmake -DCMAKE_BUILD_TYPE=Debug ..
make

# With sanitizers for memory debugging
cmake -DCMAKE_BUILD_TYPE=Debug -DENABLE_SANITIZERS=ON ..
make
```

### Code Style

The project follows modern C++ practices:
- C++20 standard
- RAII principles
- STL containers and algorithms
- OpenMP for parallelization

## 📝 Documentation

- **Doxygen**: Generate API documentation with `doxygen Doxyfile`
- **User Manual**: Detailed usage instructions in `docs/manual.pdf`
- **Theory Guide**: Mathematical background in `docs/theory.pdf`

## 🤝 Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## 📄 License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- **BARRIER Project**: Building Analysis capacity to improve the Resilience of the Guadalquivir river estuary (https://barrier.dinamicambiental.es)
- **Research Team**: Asunción Baquerizo Azofra, Manuel Díez Minguito, Agustín Millares Valenzuela
- **Funding**: This work has been carried out within the framework of the project BARRIER. Ref.: PID2023-148298OA-I00. Call 2023 of <<Proyectos de Generación de Conocimiento>>. Ministerior de Ciencia, Innovación y Universidades

## 📧 Contact

**Project Maintainer**: Manuel Cobos Budia  
**Email**: mcobosb@ugr.es 
**Institution**: Environmental Fluid Dynamics Group (https://dinamicambiental.es, University of Granada)

---

*For technical support, please open an issue on GitHub or contact the development team.*
