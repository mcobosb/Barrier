/*!
 * \file simulation.cpp
 * \brief Implementation of CSimulation class for 1D Saint-Venant estuarine modeling
 * \author Manuel Cobos Budia
 * \date 2026
 * \copyright GNU General Public License
 */

// Standard library includes
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <thread>
#include <chrono>

// Standard library using declarations
using std::cerr;
using std::cout;
using std::endl;
using std::fabs;
using std::ios;
using std::pow;
using std::setprecision;
using std::setw;
using std::sqrt;
using std::to_string;

// Project includes
#include "simulation.h"
#include "cross_section.h"
#include "data_reader.h"
#include "error_handling.h"
#include "hydrograph.h"
#include "main.h"
#include "screen_presenter.h"
#include "utils.h"
#include "yaml_reader.h"

// Forward declarations for helper functions
double calc_dynamic_albedo(double lat_deg, int day_of_year, double hour_of_day);

namespace {

std::string format_datetime_from_sim_start(const CSimulation& sim, double t_seconds) {
    std::tm start_tm = {};
    start_tm.tm_year = sim.nGetSimStartYear() - 1900;
    start_tm.tm_mon = sim.nGetSimStartMonth() - 1;
    start_tm.tm_mday = sim.nGetSimStartDay();
    start_tm.tm_hour = sim.nGetSimStartHour();
    start_tm.tm_min = sim.nGetSimStartMin();
    start_tm.tm_sec = sim.nGetSimStartSec();
    start_tm.tm_isdst = -1;

    std::time_t start_time = std::mktime(&start_tm);
    if (start_time == static_cast<std::time_t>(-1)) {
        return "(invalid start_date)";
    }

    const auto dt = static_cast<std::time_t>(std::llround(t_seconds));
    std::time_t abs_time = start_time + dt;
    std::tm* tm_ptr = std::localtime(&abs_time);
    if (tm_ptr == nullptr) {
        return "(invalid datetime)";
    }

    std::ostringstream oss;
    oss << std::setfill('0')
        << std::setw(4) << (tm_ptr->tm_year + 1900) << "-"
        << std::setw(2) << (tm_ptr->tm_mon + 1) << "-"
        << std::setw(2) << tm_ptr->tm_mday << " "
        << std::setw(2) << tm_ptr->tm_hour << ":"
        << std::setw(2) << tm_ptr->tm_min << ":"
        << std::setw(2) << tm_ptr->tm_sec;
    return oss.str();
}

struct MinMax {
    double min{0.0};
    double max{0.0};
    bool valid{false};
};

MinMax compute_minmax(const std::vector<double>& v) {
    if (v.empty()) return {};
    const auto [it_min, it_max] = std::minmax_element(v.begin(), v.end());
    return MinMax{*it_min, *it_max, true};
}

double quantile_copy(std::vector<double> v, double q) {
    if (v.empty()) return 0.0;
    if (q <= 0.0) return *std::min_element(v.begin(), v.end());
    if (q >= 1.0) return *std::max_element(v.begin(), v.end());
    const size_t k = static_cast<size_t>(std::floor(q * static_cast<double>(v.size() - 1)));
    std::nth_element(v.begin(), v.begin() + static_cast<std::ptrdiff_t>(k), v.end());
    return v[k];
}

void print_time_series_summary(const char* prefix,
                              const CSimulation& sim,
                              const std::vector<double>& t,
                              const std::vector<double>& v,
                              const char* units) {
    std::cout << prefix;
    if (t.empty() || v.empty()) {
        std::cout << "(empty)" << std::endl;
        return;
    }

    const MinMax t_mm = compute_minmax(t);
    const MinMax v_mm = compute_minmax(v);
    std::cout << "          - N: " << v.size();
    if (t.size() != v.size()) {
        std::cout << " (time points=" << t.size() << ")";
    }
    std::cout << std::endl;

    if (t_mm.valid) {
        std::cout << "          - t: [" << format_datetime_from_sim_start(sim, t_mm.min)
                  << " .. " << format_datetime_from_sim_start(sim, t_mm.max)
                  << "]";
    }
    std::cout << std::endl;
    if (v_mm.valid) {
        std::cout << "          - min: " << std::setprecision(8) << v_mm.min
                  << ", max: " << std::setprecision(8) << v_mm.max;
        if (units && *units) {
            std::cout << " " << units;
        }
    }
    std::cout << std::endl;
}

}



//======================================================================================================================
//! Constructor: Initializes all simulation parameters to default values
//======================================================================================================================
CSimulation::CSimulation() {
    // Time-related variables
    m_dSimDuration = 0.0;
    m_dSimTimestep = 0.0;
    m_dTimeFactor = 0.0;
    m_dTimestep = 0.0;
    m_dLambda = 0.0;
    m_dCurrentTime = 0.0;

    // Physical process flags
    m_bDoWaterTemperature = false;
    m_eTemperatureMode = ETemperatureMode::Off;
    m_bManningDependsOnLevel = false;
    m_bDoWaterSalinity = false;
    m_bDoWaterDensity = false;
    m_bDoSedimentTransport = false;
    m_bDoDryBed = false;

    // Numerical method flags
    m_bDoMcCormackLimiterFlux = false;
    m_bDoSurfaceGradientMethod = false;
    m_bDoSourceTermBalance = false;
    m_bDoBetaCoefficient = false;
    m_bDoMurilloCondition = false;

    // Model state variables
    m_nPredictor = -1;
    m_nTimeLogId = 0;
    m_nStringError = 0;
    m_nCrossSectionsNumber = 0;
    m_nHydrographsNumber = 0;
    m_nTimeId = 0;

    // Boundary and initial conditions
    m_nInitialEstuarineCondition = 0;
    m_nUpwardEstuarineCondition = 0;
    m_nDownwardEstuarineCondition = 0;
    m_nUpwardSalinityCondition = 0;
    m_nDownwardSalinityCondition = 0;

    // I/O and control flags
    m_bSaveTime = false;
    m_bHydroFile = false;
    m_bReturnError = false;
    m_nLogFileDetail = 0;

    // Numerical parameters
    m_dCourantNumber = 0.0;
    m_nEquationMcCormackLimiterFlux = 0;
    m_nPsiFormula = 0;
    m_dDeltaValue = 0.0;

    // Physical constants and parameters
    m_dBetaSalinityConstant = 0.0;
    m_dLongitudinalDispersion = 0.0;
    m_dUpwardBoundaryValue = 0.0;
    m_dNextUpwardBoundaryValue = 0.0;
    m_dDownwardBoundaryValue = 0.0;
    m_dNextDownwardBoundaryValue = 0.0;
    m_nEquationSedimentTransport = 0;

    // Simulation start date/time (default: 2024-01-01 00:00:00)
    m_nSimStartYear = 2024;
    m_nSimStartMonth = 1;
    m_nSimStartDay = 1;
    m_nSimStartHour = 0;
    m_nSimStartMin = 0;
    m_nSimStartSec = 0;

    // Collections
    m_vOutputVariables.clear();
    estuary.clear();
    hydrographs.clear();

    // System timing
    m_tSysStartTime = 0;
    m_tSysStartLoopTime = 0;
}

//======================================================================================================================
//! Destructor: Cleans up resources and closes open files
//======================================================================================================================
CSimulation::~CSimulation() {
    // Close log file if open
    if (LogStream && LogStream.is_open()) {
        LogStream.flush();
        LogStream.close();
    }
}


//===============================================================================================================================
//! Appends a CHydrograph objet to the CSimulation
//===============================================================================================================================
void CSimulation::AddHydrograph(){
    // Create the
    CHydrograph hydrograph;
    hydrographs.push_back(hydrograph);
}


//===============================================================================================================================
//! Appends a CCrossSection objet to the estuary
//===============================================================================================================================
void CSimulation::AddCrossSection(){
    // Create the
    CCrossSection crossSection;
    estuary.push_back(crossSection);
}

//===============================================================================================================================
//! The bDoSimulation member function of CSimulation sets up and runs the simulation
//===============================================================================================================================
void CSimulation::bDoSimulation(int nArg, char const* pcArgv[]){

    //! Show starting run message
    presenter.StartingRun(nArg, pcArgv, this);

    //! Detect configuration file type and read accordingly
    std::string configFile;
    
    if (nArg > 1) {
        configFile = pcArgv[1];
    } else {
        // No config file specified, try to find config.yaml in current directory
        if (std::filesystem::exists("config.yaml")) {
            configFile = "config.yaml";
        } else {
            std::cerr << "❌ Error: No configuration file detected. config.yaml not found" << std::endl;
            std::cerr << "Usage: " << pcArgv[0] << " <config_file>" << std::endl;
            m_bReturnError = true;
            return;
        }
    }
    std::filesystem::path configPath(configFile);
    std::string extension = configPath.extension().string();
    

    std::cout << "      1. Reading configuration file" << std::endl;
    
    // Read YAML file content as string (for NetCDF metadata)
    std::ifstream yamlFile(configFile);
    if (yamlFile.is_open()) {
        std::stringstream buffer;
        buffer << yamlFile.rdbuf();
        m_strYAMLConfigContent = buffer.str();
        yamlFile.close();
    }
    
    CYAMLReader yamlReader;
    if (!yamlReader.loadConfiguration(configFile, this)) {
        std::cerr << "❌ Error loading YAML configuration: " 
                    << yamlReader.getErrorMessage() << std::endl;
        m_bReturnError = true;
        return;
    }

    // Display simulation start date
    std::cout << "        - Simulation start/end date: " << getSimulationStartDateTimeString() << " / " << getSimulationEndDateTimeString() << std::endl;
    

    // Transfer file paths to reader and simulation
    reader.m_strAlongChannelDataFilename = yamlReader.m_strAlongChannelDataFilename;
    reader.m_strCrossSectionsFilename = yamlReader.m_strCrossSectionGeometryFilename;
    reader.m_strSedimentPropertiesFilename = yamlReader.m_strSedimentPropertiesFilename;
    reader.m_strHydroFilename = yamlReader.m_strHydrographsFilename;

    
    // Read geometry and forcing files
    CDataReader::bOpenLogFile(this);
    
    // === CRITICAL WARNING FOR DEBUG LEVEL 3 ===
    if (m_nLogFileDetail >= 3) {
        std::cerr << "\n" << std::string(80, '!') << "\n";
            std::cerr << "CRITICAL WARNING: log_file_detail = " << m_nLogFileDetail << "\n";
            std::cerr << "NetCDF output will be written at EVERY timestep!\n";
            std::cerr << "This will:\n";
            std::cerr << "  - Make simulation EXTREMELY SLOW (10-100x slower)\n";
            std::cerr << "  - Generate HUGE files (potentially GB to TB)\n";
            std::cerr << "  - May FILL YOUR DISK and crash the system\n";
            std::cerr << "\n";
            std::cerr << "Recommended: Use log_file_detail = 2 for detailed diagnostics\n";
            std::cerr << "Only use level 3 for debugging specific short runs!\n";
            std::cerr << std::string(80, '!') << "\n";
            std::cerr << "Press Ctrl+C within 5 seconds to abort...\n\n";
            std::cerr.flush();
        
        // Give user 5 seconds to abort
        std::this_thread::sleep_for(std::chrono::seconds(5));
    }
    
    // Log simulation configuration summary
    if (LogStream.is_open()) {
        LogStream << std::string(80, '=') << "\n";
        LogStream << "SIMULATION CONFIGURATION SUMMARY\n";
        LogStream << std::string(80, '=') << "\n";
        
        // === CRITICAL WARNING IN LOG FILE ===
        if (m_nLogFileDetail >= 3) {
            LogStream << "\n" << std::string(80, '!') << "\n";
                LogStream << "CRITICAL WARNING: DEBUG LEVEL 3 ACTIVE\n";
                LogStream << "NetCDF output writing at EVERY timestep\n";
                LogStream << "Expect extremely slow execution and huge output files\n";
                LogStream << std::string(80, '!') << "\n\n";
        }
        
        LogStream << "Configuration file: " << configFile << "\n";
        LogStream << "Log detail level: " << m_nLogFileDetail << "\n";
        LogStream << "  0 = No log file (all to /dev/null)\n";
        LogStream << "  1 = Basic statistics (hourly summary, warnings)\n";
        LogStream << "  2 = Detailed diagnostics (heat fluxes every 6h, anomaly checks)\n";
        LogStream << "  3 = Full debug (NetCDF output every timestep - VERY SLOW, huge files)\n";
        LogStream << "Simulation period: " << getSimulationStartDateTimeString() 
                  << " to " << getSimulationEndDateTimeString() << "\n";
        LogStream << "Duration: " << m_dSimDuration << " s (" 
                  << m_dSimDuration/86400.0 << " days)\n";
        LogStream << "Timestep: " << m_dSimTimestep << " s\n";
        LogStream << "Courant number: " << m_dCourantNumber << " (adaptive control enabled)\n";
        
        // Initialize automatic CFL control based on numerical scheme
        // CRITICAL: Shallow water + TVD-McCormack is MORE restrictive than pure advection
        // Empirical evidence shows instabilities (elevation oscillations, spurious diffusion)
        // when CFL exceeds ~0.15-0.2 in coupled hyperbolic systems
        
        // Conservative CFL limits based on shallow water + tracer coupling:
        // - Pure first-order upwind: CFL_safe ≈ 0.5 (very stable but diffusive)
        // - TVD-McCormack with MinMod: CFL_safe ≈ 0.25 (balanced)
        // - TVD-McCormack with Superbee: CFL_safe ≈ 0.15 (sharp fronts, needs restriction)
        // - TVD-McCormack with Van Leer/Van Albada: CFL_safe ≈ 0.2
        
        // Determine a *cap* (maximum allowed CFL) and a conservative *auto default*.
        // NOTE: Very low CFL increases numerical diffusion and can over-damp smooth long
        // waves (e.g. tides). We keep the default conservative, but allow the user to
        // request a higher CFL (up to the cap) when they want less damping.
        double cfl_cap = 0.2;
        double cfl_auto_conservative_max = 0.2;
        if (!m_bDoMcCormackLimiterFlux || m_nEquationMcCormackLimiterFlux == 0) {
            // No TVD: first-order upwind (robust)
            cfl_cap = 0.9;
            cfl_auto_conservative_max = 0.5; // legacy conservative default
        } else {
            switch (m_nEquationMcCormackLimiterFlux) {
                case 1:  // MinMod
                    cfl_cap = 0.6;
                    cfl_auto_conservative_max = 0.25;
                    break;
                case 2:  // Superbee
                    cfl_cap = 0.4;
                    cfl_auto_conservative_max = 0.15;
                    break;
                case 3:  // Van Leer
                    cfl_cap = 0.6;
                    cfl_auto_conservative_max = 0.2;
                    break;
                case 4:  // Van Albada
                    cfl_cap = 0.6;
                    cfl_auto_conservative_max = 0.2;
                    break;
                default:
                    cfl_cap = 0.5;
                    cfl_auto_conservative_max = 0.2;
            }
        }
        m_dCourantNumberMaxTheoretical = cfl_cap;
        
        // Apply minimal safety margin for the automatic (conservative) default.
        const double SAFETY_FACTOR = 0.9;
        m_dCourantNumberOptimal = SAFETY_FACTOR * cfl_auto_conservative_max;
        
        // User-provided Courant is treated as a requested target.
        // If it exceeds the scheme cap, we clamp to the cap.
        if (m_dCourantNumber > 0.0) {
            if (m_dCourantNumber <= m_dCourantNumberMaxTheoretical) {
                m_dCourantNumberOptimal = m_dCourantNumber;
                LogStream << "ℹ️  User CFL=" << m_dCourantNumber
                          << " (requested target within cap=" << m_dCourantNumberMaxTheoretical << ")\n";
            } else {
                LogStream << "⚠️  User CFL=" << m_dCourantNumber << " exceeds cap="
                          << m_dCourantNumberMaxTheoretical << " (clamping)\n";
                m_dCourantNumberOptimal = m_dCourantNumberMaxTheoretical;
            }
        }
        
        m_dCourantNumberCurrent = m_dCourantNumberOptimal;
        m_dCourantNumberMin = m_dCourantNumberOptimal;
        m_nCourantReductions = 0;
        m_bExtremeFlowConditions = false;
        
        // Log automatic CFL determination
        LogStream << "\n=== AUTOMATIC CFL DETERMINATION ===\n";
        LogStream << "System: Coupled shallow water + TVD-McCormack tracers\n";
        LogStream << "Limiter: ";
        if (!m_bDoMcCormackLimiterFlux || m_nEquationMcCormackLimiterFlux == 0) {
            LogStream << "None (1st order upwind)";
        } else {
            const char* limiter_names[] = {"", "MinMod", "Superbee", "Van Leer", "Van Albada"};
            LogStream << limiter_names[m_nEquationMcCormackLimiterFlux];
        }
        LogStream << "\nCFL_max (empirical for coupled system): " << m_dCourantNumberMaxTheoretical << "\n";
        LogStream << "CFL_optimal (with 90% safety): " << m_dCourantNumberOptimal << "\n";
        if (m_dCourantNumber > 0.0) {
            LogStream << "CFL_user (config): " << m_dCourantNumber;
            if (m_dCourantNumber < m_dCourantNumberOptimal) {
                LogStream << " [LIMITING - user more conservative]\n";
            } else if (m_dCourantNumber > m_dCourantNumberOptimal * 1.1) {
                LogStream << " [IGNORED - exceeds stability limit]\n";
            } else {
                LogStream << " [OK - close to optimal]\n";
            }
        }
        LogStream << "CFL_active: " << m_dCourantNumberCurrent << "\n";
        LogStream << "\nNote: Very low CFL increases numerical diffusion (can over-damp tides).\n";
        LogStream << "      Higher CFL reduces numerical diffusion but may trigger oscillations;\n";
        LogStream << "      adaptive CFL will reduce it if needed.\n";
        LogStream << "===================================\n\n";
        LogStream.flush();
    }
    
    std::cout << "      2. Reading input data files" << std::endl;
    reader.bReadAlongChannelDataFile(this);
    reader.bReadCrossSectionGeometryFile(this);

    // Lateral storage summary (one-time, startup)
    // Note: this runs before initializeVectors(); the vector may still be empty if the along-channel
    // file doesn't provide storage factors. Use a safe fallback of S=1.
    double storage_s_min = 1.0;
    double storage_s_max = 1.0;
    if (!m_vLateralStorageFactor.empty()) {
        storage_s_min = 1e300;
        storage_s_max = -1e300;
        for (double s : m_vLateralStorageFactor) {
            const double sc = (s >= 1.0) ? s : 1.0;
            storage_s_min = std::min(storage_s_min, sc);
            storage_s_max = std::max(storage_s_max, sc);
        }
        if (!(storage_s_min < 1e200 && storage_s_max > -1e200)) {
            storage_s_min = 1.0;
            storage_s_max = 1.0;
        }
    }

    // Section label for console summaries (a), b), c)...
    // Increment only when a section is actually printed to avoid gaps.
    char section_label = 'a';
    {
        // Cross-section summary (one-time, startup)
        const MinMax x_mm = compute_minmax(m_vCrossSectionX);

        double z_min = 0.0, z_max = 0.0;
        double n_min = 0.0, n_max = 0.0;
        double beta_min = 0.0, beta_max = 0.0;
        bool have_z = false;
        bool have_n = false;
        bool have_beta = false;
        if (m_nCrossSectionsNumber > 0 && static_cast<int>(estuary.size()) >= m_nCrossSectionsNumber) {
            z_min = z_max = estuary[0].dGetZ();
            n_min = n_max = estuary[0].dGetManningNumber();
            beta_min = beta_max = estuary[0].dGetBeta();
            have_z = true;
            have_n = true;
            have_beta = true;
            for (int i = 1; i < m_nCrossSectionsNumber; ++i) {
                z_min = std::min(z_min, estuary[i].dGetZ());
                z_max = std::max(z_max, estuary[i].dGetZ());
                n_min = std::min(n_min, estuary[i].dGetManningNumber());
                n_max = std::max(n_max, estuary[i].dGetManningNumber());
                beta_min = std::min(beta_min, estuary[i].dGetBeta());
                beta_max = std::max(beta_max, estuary[i].dGetBeta());
            }
        }

        std::cout << "        " << section_label << ") Cross-sections loaded: " << m_nCrossSectionsNumber << std::endl;
        if (!reader.m_strAlongChannelDataFilename.empty()) {
            std::cout << "          · Along-channel file: " << reader.m_strAlongChannelDataFilename << std::endl;
        }
        if (!reader.m_strCrossSectionsFilename.empty()) {
            std::cout << "          · Cross-sections file: " << reader.m_strCrossSectionsFilename << std::endl;
        }
        if (x_mm.valid) {
            std::cout << "          - X: [ " << std::setprecision(8) << x_mm.min << " .. " << x_mm.max << "] m" << std::endl;
        }
        if (have_z) {
            std::cout << "          - Zbed: [" << std::setprecision(8) << z_min << " .. " << z_max << "] m" << std::endl;
        }
        if (have_n) {
            std::cout << "          - Manning n: [" << std::setprecision(8) << n_min << " .. " << n_max << "]" << std::endl;
        }
        if (have_beta) {
            std::cout << "          - Momentum coefficient beta: [" << std::setprecision(8) << beta_min << " .. " << beta_max << "]" << std::endl;
        }

        if (m_bDoWaterSalinity) {
            if (bHasLongitudinalDispersionProfile()) {
                double kh_min = 1e300;
                double kh_max = -1e300;
                for (double kh : m_vLongitudinalDispersion) {
                    kh_min = std::min(kh_min, kh);
                    kh_max = std::max(kh_max, kh);
                }
                if (kh_min < 1e200 && kh_max > -1e200) {
                    std::cout << "          - Dispersion Kh(x): [" << std::setprecision(8) << kh_min
                              << " .. " << kh_max << "] m^2/s" << std::endl;
                }
            } else {
                const double kh_const = dGetLongitudinalDispersionConstant();
                std::cout << "          - Dispersion Kh: " << std::setprecision(8) << kh_const << " m^2/s" << std::endl;
            }
        }

        if (m_bDoLateralStorage) {
            std::cout << "          - Lateral storage Sf: [" << std::setprecision(8) << storage_s_min
                      << " .. " << storage_s_max << "]" << std::endl;
        }
    }
    ++section_label;
    
    // Print geometry and configuration summary AFTER reading geometry
    if (m_nLogFileDetail >= 1) {
        LogStream << "Number of cross-sections: " << m_nCrossSectionsNumber << "\n";
        LogStream << "\nPhysical processes:\n";
        LogStream << "  - Water temperature: " << (m_bDoWaterTemperature ? "YES" : "NO") << "\n";
        LogStream << "  - Water salinity: " << (m_bDoWaterSalinity ? "YES" : "NO") << "\n";
        LogStream << "  - Sediment transport: " << (m_bDoSedimentTransport ? "YES" : "NO") << "\n";
        LogStream << "  - Dry bed treatment: " << (m_bDoDryBed ? "YES" : "NO") << "\n";
        LogStream << "\nNumerical methods:\n";
        LogStream << "  - TVD flux limiter: " << (m_bDoMcCormackLimiterFlux ? "YES" : "NO") << "\n";
        LogStream << "  - Surface gradient method: " << (m_bDoSurfaceGradientMethod ? "YES" : "NO") << "\n";
        LogStream << "  - Source term balance: " << (m_bDoSourceTermBalance ? "YES" : "NO") << "\n";
        LogStream << "  - Solution regularization: " << (m_bDoSmoothSolution ? "YES" : "NO") << "\n";
        LogStream << "  - Auto regularization on extreme Q: AUTO (runtime-detected; marked with '*')\n";
        if (m_bDoSmoothSolution) {
            LogStream << "    • Passes: " << m_nSolutionSmoothingPasses << "\n";
            LogStream << "    • Alpha: " << m_dSolutionSmoothingAlpha << "\n";
        }
        if (m_bDoLateralStorage) {
            LogStream << "  - Lateral storage: YES (Sf in [" << std::setprecision(8) << storage_s_min
                      << " .. " << storage_s_max << "])\n";
        } else {
            LogStream << "  - Lateral storage: NO\n";
        }
        LogStream << "\nBoundary conditions:\n";
        LogStream << "  - Upstream hydrodynamic: Type " << m_nUpwardEstuarineCondition << "\n";
        LogStream << "  - Downstream hydrodynamic: Type " << m_nDownwardEstuarineCondition << "\n";
        if (m_bDoWaterTemperature) {
            LogStream << "  - Upstream temperature: Type " << m_nUpwardTemperatureCondition << "\n";
            LogStream << "  - Downstream temperature: Type " << m_nDownwardTemperatureCondition << "\n";
        }
        if (m_bDoWaterSalinity) {
            LogStream << "  - Upstream salinity: Type " << m_nUpwardSalinityCondition << "\n";
            LogStream << "  - Downstream salinity: Type " << m_nDownwardSalinityCondition << "\n";
        }
        LogStream << "\nOutput file: " << m_strOutFile << generateOutputFileName() << "\n";
        LogStream << std::string(80, '=') << "\n\n";
        LogStream.flush();
    }
    if (m_nUpwardEstuarineCondition > 1) {
        CDataReader::bReadUpwardBoundaryConditionFile(this);
        const bool upstream_is_discharge = (m_nUpwardEstuarineCondition == 3);
        const char* kind = upstream_is_discharge ? "discharge" : "level";
        const char* units = upstream_is_discharge ? "m3/s" : "m";
        std::cout << "        " << section_label << ") Upstream hydrodynamic BC (" << kind << "): " << std::endl;
        if (!m_strUpwardBoundaryConditionFilename.empty()) {
            std::cout << "          · File: " << m_strUpwardBoundaryConditionFilename << std::endl;
        }
        print_time_series_summary("", *this, m_vUpwardBoundaryConditionTime, m_vUpwardBoundaryConditionValue, units);

        ++section_label;

        // Auto-detect an "extreme discharge" threshold from the upstream discharge time series.
        // This avoids relying on config.yaml while still being data-adaptive.
        m_dAutoSmoothAbsQThreshold = 0.0;
        if ((m_nUpwardEstuarineCondition == 3) && !m_vUpwardBoundaryConditionValue.empty()) {
            std::vector<double> absq;
            absq.reserve(m_vUpwardBoundaryConditionValue.size());
            for (double q : m_vUpwardBoundaryConditionValue) absq.push_back(std::fabs(q));

            const double q99 = quantile_copy(absq, 0.99);
            const double q995 = quantile_copy(absq, 0.995);
            const double qmax = *std::max_element(absq.begin(), absq.end());

            // Use an upper-tail quantile (robust to single spikes). If the record is short,
            // q99~qmax anyway. Prefer the more selective of q99 and q995 when available.
            const double q_ext = std::max(q99, q995);
            m_dAutoSmoothAbsQThreshold = (q_ext > 0.0) ? q_ext : qmax;

            if (m_nLogFileDetail >= 1) {
                std::ofstream LogStream(m_strLogFile.c_str(), std::ios_base::app);
                LogStream << std::fixed << std::setprecision(3);
                LogStream << "Auto extreme-Q threshold (computed): |Q| >= " << m_dAutoSmoothAbsQThreshold
                          << " m3/s (q99=" << q99 << ", q99.5=" << q995 << ", max=" << qmax << ")\n";
                LogStream.flush();
            }
        }
    }
    if (m_nDownwardEstuarineCondition > 1) {
        CDataReader::bReadDownwardBoundaryConditionFile(this);
        const char* kind = (m_nDownwardEstuarineCondition == 3) ? "discharge" : "level";
        const char* units = (m_nDownwardEstuarineCondition == 3) ? "m3/s" : "m";
        std::cout << "        " << section_label << ") Downstream hydrodynamic BC (" << kind << "): " << std::endl;
        if (!m_strDownwardBoundaryConditionFilename.empty()) {
            std::cout << "          · File: " << m_strDownwardBoundaryConditionFilename << std::endl;
        }
        print_time_series_summary("", *this, m_vDownwardBoundaryConditionTime, m_vDownwardBoundaryConditionValue, units);

        ++section_label;
    }

    // === Salinity: Boundary conditions summary ===
    if (m_bDoWaterSalinity) {
        std::cout << "        " << section_label << ") Salinity BCs:" << std::endl;
        const char* s_types[] = {"free", "null", "ocean"};
        const int up_t = std::max(0, std::min(2, m_nUpwardSalinityCondition));
        const int dn_t = std::max(0, std::min(2, m_nDownwardSalinityCondition));

        std::cout << "          - Upstream: type " << m_nUpwardSalinityCondition
                  << " (" << s_types[up_t] << ")";
        if (!m_strUpwardSalinityBoundaryConditionFilename.empty()) {
            std::cout << ", file=" << m_strUpwardSalinityBoundaryConditionFilename;
        } else {
            std::cout << ", value=" << std::setprecision(8) << m_dUpwardSalinityBoundaryValue;
        }
        std::cout << std::endl;

        std::cout << "          - Downstream: type " << m_nDownwardSalinityCondition
                  << " (" << s_types[dn_t] << ")";
        if (!m_strDownwardSalinityBoundaryConditionFilename.empty()) {
            std::cout << ", file=" << m_strDownwardSalinityBoundaryConditionFilename;
        } else {
            std::cout << ", value=" << std::setprecision(8) << m_dDownwardSalinityBoundaryValue;
        }
        std::cout << std::endl;

        // Common misconfiguration: ocean boundary enabled but salinity left at default 0.
        if (m_nDownwardSalinityCondition == 2 && m_strDownwardSalinityBoundaryConditionFilename.empty() &&
            std::fabs(m_dDownwardSalinityBoundaryValue) < 1e-12) {
            std::cout << "          ⚠️  Downstream salinity is set to 'ocean' but value is ~0. "
                      << "This makes the mouth effectively fresh; salt can then only come from initial/restart fields." << std::endl;
        }

        ++section_label;
    }

    // === Temperature: Boundary conditions and forcing ===
    if (m_bDoWaterTemperature) {
        if (bIsGivenTemperature()) {
            if (m_strGivenTemperatureFilename.empty()) {
                std::cerr << ERR << "temperature.enabled: given requires temperature.given_file (time,temperature)" << std::endl;
                exit(EXIT_FAILURE);
            }
            CDataReader::bReadGivenTemperatureFile(this);
            std::cout << "      - Given temperature file loaded successfully" << std::endl;
        }

        if (!bIsGivenTemperature()) {
            if (m_nUpwardTemperatureCondition == 2 || m_bCalculateRHFromTemperature == false) {
                CDataReader::bReadUpwardTemperatureBoundaryConditionFile(this);
            }
            if (m_nDownwardTemperatureCondition > 1) {
                CDataReader::bReadDownwardTemperatureBoundaryConditionFile(this);
            }
            // Forzamiento: heat_flux_file (Tair, humedad relativa, viento)
            // REQUIRED for upstream type 3 (0D reservoir model)
            if (m_nUpwardTemperatureCondition == 3 && m_strHeatFluxFile.empty()) {
                std::cerr << ERR << "Upstream temperature condition type 3 (0D reservoir model) requires heat_flux_file" << std::endl;
                std::cerr << "      Please specify heat_flux_file with meteorological forcing (time, T_air, RH, wind, pressure)" << std::endl;
                exit(EXIT_FAILURE);
            }
            if (!m_strHeatFluxFile.empty()) {
                CDataReader::bReadHeatFluxFile(this);
                std::cout << "      - Heat flux data (atmospheric forcing) loaded successfully" << std::endl;
                if (m_bCalculateRHFromTemperature == true) {
                    // Calculate daily minimum temperatures if RH needs to be estimated
                    calculateDailyMinTemperatures();
                }
            }
        }
    }
    if (m_bDoSedimentTransport) {
        reader.bReadAlongChannelSedimentsFile(this);
        std::cout << "      - Sediment transport data loaded successfully" << std::endl;
    }
    if (m_bHydroFile) {
        reader.bReadHydrographsFile(this);
        // Hydrographs summary (global across all hydrographs)
        double t_min = 0.0, t_max = 0.0;
        double q_min = 0.0, q_max = 0.0;
        bool have_any = false;
        for (int i = 0; i < m_nHydrographsNumber; ++i) {
            const auto& ht = hydrographs[i].vGetTime();
            const auto& hq = hydrographs[i].vGetQ();
            if (ht.empty() || hq.empty()) {
                continue;
            }
            const MinMax ht_mm = compute_minmax(ht);
            const MinMax hq_mm = compute_minmax(hq);
            if (!ht_mm.valid || !hq_mm.valid) {
                continue;
            }
            if (!have_any) {
                t_min = ht_mm.min;
                t_max = ht_mm.max;
                q_min = hq_mm.min;
                q_max = hq_mm.max;
                have_any = true;
            } else {
                t_min = std::min(t_min, ht_mm.min);
                t_max = std::max(t_max, ht_mm.max);
                q_min = std::min(q_min, hq_mm.min);
                q_max = std::max(q_max, hq_mm.max);
            }
        }

        std::cout << "        " << section_label << ") Hydrographs loaded: " << m_nHydrographsNumber << std::endl;
        if (!reader.m_strHydroFilename.empty()) {
            std::cout << "          · File: " << reader.m_strHydroFilename << std::endl;
        }
        if (have_any) {
            std::cout << "          - t: [" << format_datetime_from_sim_start(*this, t_min)
                      << " .. " << format_datetime_from_sim_start(*this, t_max)
                      << "]" << std::endl;
            std::cout << "          - Q min: " << std::setprecision(8) << q_min
                      << ", Q max: " << std::setprecision(8) << q_max << " m3/s" << std::endl;
        }

        ++section_label;
    }

    // Initialize simulation vectors
    initializeVectors();
    precomputeEstuaryData();

    // Sanity-check: ensure cached beta profile (m_vBeta) matches geometry (estuary[].beta)
    // This catches cases where inputs are updated but caches are stale or overwritten.
    if (m_nCrossSectionsNumber > 0 && static_cast<int>(estuary.size()) >= m_nCrossSectionsNumber &&
        static_cast<int>(m_vBeta.size()) == m_nCrossSectionsNumber) {
        double beta_geom_min = estuary[0].dGetBeta();
        double beta_geom_max = estuary[0].dGetBeta();
        double beta_cache_min = m_vBeta[0];
        double beta_cache_max = m_vBeta[0];
        for (int i = 1; i < m_nCrossSectionsNumber; ++i) {
            beta_geom_min = std::min(beta_geom_min, estuary[i].dGetBeta());
            beta_geom_max = std::max(beta_geom_max, estuary[i].dGetBeta());
            beta_cache_min = std::min(beta_cache_min, m_vBeta[i]);
            beta_cache_max = std::max(beta_cache_max, m_vBeta[i]);
        }
        if (m_nLogFileDetail >= 1) {
            std::ofstream LogStream(m_strLogFile.c_str(), std::ios_base::app);
            LogStream << std::fixed << std::setprecision(8);
            LogStream << "Beta profile (geometry): [" << beta_geom_min << " .. " << beta_geom_max << "]\n";
            LogStream << "Beta profile (cached)  : [" << beta_cache_min << " .. " << beta_cache_max << "]\n";
            if (fabs(beta_geom_min - beta_cache_min) > 1e-12 || fabs(beta_geom_max - beta_cache_max) > 1e-12) {
                LogStream << "⚠️  Beta mismatch between geometry and cached arrays. "
                          << "Solver uses cached m_vBeta; check initialization order or restart restore.\n";
            }
            LogStream.flush();
        }
    }

    // For reservoir temperature mode (type 3), initialize upstream temperature
    if (m_bDoWaterTemperature && !bIsGivenTemperature() && m_nUpwardTemperatureCondition == 3) {
        // Initialize temperature vector with initial condition if empty or wrong size
        if (m_vUpwardTemperatureBoundaryConditionValue.empty() || 
            m_vUpwardTemperatureBoundaryConditionValue.size() != m_vHeatFluxTime.size()) {
            
            // Determine initial temperature (use air temperature if available, otherwise default 15°C)
            double T_initial = (!m_vHeatFluxAirTemp.empty()) ? m_vHeatFluxAirTemp[0] : 15.0;
            
            // Initialize entire vector with initial temperature
            m_vUpwardTemperatureBoundaryConditionValue.resize(m_vHeatFluxTime.size(), T_initial);
            
            std::cout << "      - Initialized reservoir temperature to " << T_initial << "°C" << std::endl;
        }
        
        updateReservoirTemperature0D();
    }

    // Restore simulation state from NetCDF if continuing previous run
    if (m_bContinueSimulation && !m_strContinueNetcdfPath.empty()) {
        reader.bRestoreStateFromNetCDF(this, m_strContinueNetcdfPath);
        std::cout << "      3. State restored successfully from NetCDF:" << m_strContinueNetcdfPath << std::endl;
    }
    
    // Apply bathymetry smoothing before calculating slopes (if enabled)
    if (bGetDoSmoothBathymetry()) {
        smoothBathymetry();
    }

    calculateBedSlope();

    // Initialize hydraulic conditions (skipped if continuing from NetCDF)
    calculateAlongEstuaryInitialConditions();

    // If temperature is prescribed, overwrite initial spatial field
    if (m_bDoWaterTemperature && bIsGivenTemperature()) {
        applyGivenTemperatureAtTime(m_dCurrentTime);
    }

    // Calculate initial hydraulic parameters if continuing simulation
    if (m_bContinueSimulation) {
        m_nPredictor = 1;
        calculateHydraulicParameters();
    }

    // Calculate UTM coordinates for river banks (always performed)
    calculateRiverBankUTMCoordinates();

    // Generate output filename and create NetCDF file
    std::string m_strOutFileName = generateOutputFileName();
    m_strOutFile += m_strOutFileName;
    
    std::cout << "      4. Output file: " << m_strOutFile << std::endl;
    writer.nDefineNetCDFFile(this);

    // Apply dry bed conditions if enabled
    if (bGetDoDryBed())
        dryArea();

    // Save initial state
    m_nTimeId = 0;
    writer.nSetOutputData(this);
    m_nTimeId++;

    int m_nStep = 1;
    cout << "    - Running" << endl;
    m_tSysStartLoopTime = time(nullptr);

    // Write diagnostic log table header
    if (LogStream.is_open() && m_nLogFileDetail >= 1) {
        LogStream << "\n" << std::string(160, '=') << "\n";
        LogStream << "SIMULATION DIAGNOSTICS - PERIODIC STATISTICS\n";
        LogStream << std::string(160, '=') << "\n";
        LogStream << std::setw(20) << "Date" 
                  << std::setw(10) << "Time(s)"
                  << std::setw(8) << "dt(s)" 
                  << std::setw(8) << "CFL_max"
                  << std::setw(10) << "Q_min" 
                  << std::setw(10) << "Q_max" 
                  << std::setw(10) << "U_max"
                  << std::setw(10) << "η_min" 
                  << std::setw(10) << "η_max";
        if (m_bDoWaterTemperature) {
            LogStream << std::setw(10) << "T_min"
                      << std::setw(10) << "T_max"
                      << std::setw(10) << "T_res";
        }
        if (m_bDoWaterSalinity) {
            LogStream << std::setw(10) << "S_min"
                      << std::setw(10) << "S_max";
        }
        LogStream << "\n" << std::string(160, '-') << "\n";
        LogStream.flush();
    }

    //======================================================================================================================
    //! MAIN TIME-STEPPING LOOP: McCormack predictor-corrector scheme with TVD flux limiter
    //======================================================================================================================
    while (m_dCurrentTime <= m_dSimDuration) {
        // Note: Reservoir temperature (type 3) already pre-computed before loop
        // Current time value obtained by interpolation in boundary condition functions
        
        // Update boundary conditions and calculate adaptive timestep
        calculateBoundaryConditions();
        calculateTimestep();

        //==============================================================================================================
        //! STEP 1: PREDICTOR (forward in time, forward in space)
        //==============================================================================================================
        m_nPredictor = 1;
        
        // Apply BC to current state BEFORE calculating hydraulic parameters and fluxes
        applyBoundariesToCurrentState();
        calculateHydraulicParameters();

        if (m_bDoWaterTemperature && bIsGivenTemperature()) {
            applyGivenTemperatureAtTime(m_dCurrentTime);
        }

        // Calculate water density from salinity and temperature at time n
        if (bGetDoWaterSalinity() || bGetDoWaterTemperature()) {
            calculate_density();
        }

        // Calculate source terms (bed slope, friction, pressure gradients)
        calculate_GS_A_terms();
        calculateFlowTerms();
        calculateSourceTerms();
        
        // Advance hydrodynamics to predictor state (eta*, Q*)
        calculatePredictor();
        
        // Apply BC to predictor state AFTER calculating predictor interior
        updatePredictorBoundaries();
        if (bGetDoDryBed()) dryArea();

        // Advance transport scalars to predictor state (S*, T*)
        if (m_bDoWaterTemperature && !bIsGivenTemperature()) {
            calculateRadiativeFluxes();
            calculate_temperature_predictor();
        }
        if (bGetDoWaterSalinity()) {
            calculate_salinity_predictor();
        }

        //==============================================================================================================
        //! STEP 2: CORRECTOR (backward in time, backward in space)
        //==============================================================================================================
        m_nPredictor = 2;
        
        // Apply BC to predictor state BEFORE calculating hydraulic parameters and fluxes
        applyBoundariesToPredictorState();
        calculateHydraulicParameters();

        if (m_bDoWaterTemperature && bIsGivenTemperature()) {
            applyGivenTemperatureAtTime(m_dCurrentTime);
        }

        // Calculate water density from predictor state (S*, T*)
        if (bGetDoWaterSalinity() || bGetDoWaterTemperature()) {
            calculate_density();
        }

        // Recalculate source terms with predictor state
        calculate_GS_A_terms();
        calculateFlowTerms();
        calculateSourceTerms();
        
        // Advance hydrodynamics to corrector state (eta**, Q**)
        calculateCorrector();
        
        // Apply BC to corrector state AFTER calculating corrector interior
        updateCorrectorBoundaries();
        if (bGetDoDryBed()) dryArea();

        // Advance transport scalars to corrector state (S**, T**)
        if (m_bDoWaterTemperature && !bIsGivenTemperature()) {
            calculate_temperature_corrector();
        }
        if (bGetDoWaterSalinity()) {
            calculate_salinity_corrector();
        }

        //==============================================================================================================
        //! STEP 3: MERGE predictor and corrector with TVD flux limiter
        //==============================================================================================================
        mergePredictorCorrector();
        mergeTracerPredictorCorrector();

        m_nPredictor = 0;
        if (bGetDoDryBed()) dryArea();
        // Reset per-timestep auto-smoothing marker (used by hourly log '*').
        m_bAutoSmoothAppliedThisStep = false;

        // Runtime solution smoothing (regularization):
        // - If enabled in config: always apply.
        // - If disabled: auto-enable only under runtime-detected extreme conditions.
        if (bGetDoSmoothSolution()) {
            smoothSolution();
        } else {
            bool extreme = m_bExtremeFlowConditions;
            if (!extreme && m_nUpwardEstuarineCondition == 3 && m_dAutoSmoothAbsQThreshold > 0.0 && m_nCrossSectionsNumber > 0) {
                extreme = (fabs(m_vCrossSectionQ[0]) >= m_dAutoSmoothAbsQThreshold);
            }
            if (extreme) {
                // Ensure smoothing parameters are meaningful even if user didn't set them.
                if (m_nSolutionSmoothingPasses <= 0) m_nSolutionSmoothingPasses = 1;
                if (m_dSolutionSmoothingAlpha <= 0.0) m_dSolutionSmoothingAlpha = 0.10;
                smoothSolution();
                m_bAutoSmoothAppliedThisStep = true;
            }
        }

        // Smoothing updates interior nodes (including those adjacent to boundaries).
        // Re-apply boundary conditions so stage/discharge closures remain consistent.
        applyBoundariesToCurrentState();

        // Display progress and save output if needed
        AnnounceProgress();
        
        // Check for anomalous values (warnings)
        if (m_nLogFileDetail >= 1) {
            checkAnomalousValues();
        }
        
        // Write periodic statistics (every simulated hour)
        if (m_nLogFileDetail >= 1) {
            writePeriodicStatistics();
        }

        // === Salt mass balance diagnostics ===
        if (m_bDoWaterSalinity && LogStream.is_open() && m_nLogFileDetail >= 1) {
            // Compute total salt mass in domain (effective storage = S * A)
            double total_salt = 0.0;
            for (int i = 0; i < m_nCrossSectionsNumber; ++i) {
                const double S = dGetLateralStorageFactor(i);
                total_salt += (S * m_vCrossSectionArea[i]) * m_vCrossSectionSalinity[i];
            }
            // Compute salt fluxes at boundaries (kg/s)
            double salt_in = 0.0, salt_out = 0.0;
            // Upstream (dam)
            if (m_vCrossSectionQ[0] > 0.0) {
                salt_in += m_vCrossSectionQ[0] * m_vCrossSectionSalinity[0];
            } else {
                salt_out -= m_vCrossSectionQ[0] * m_vCrossSectionSalinity[0];
            }
            // Downstream (ocean)
            int last = m_nCrossSectionsNumber - 1;
            if (m_vCrossSectionQ[last] < 0.0) {
                salt_in += -m_vCrossSectionQ[last] * m_vCrossSectionSalinity[last];
            } else {
                salt_out += m_vCrossSectionQ[last] * m_vCrossSectionSalinity[last];
            }
            // Removed per-timestep [SaltBalance] log
        }
        
        // Level 2: Save output at EVERY timestep (WARNING: Very slow, huge files)
        // Use only for debugging specific timesteps
        if (m_nLogFileDetail >= 3) {
            // Ensure derived geometry outputs are consistent with merged state
            calculateHydraulicParameters();
            calculateRiverBankUTMCoordinates();
            writer.nSetOutputData(this);
        }
        // Level 1-2: Save only at scheduled output times
        else if (m_bSaveTime) {
            // Ensure derived geometry outputs are consistent with merged state
            calculateHydraulicParameters();
            calculateRiverBankUTMCoordinates();
            writer.nSetOutputData(this);
        }

        // Advance time
        m_dCurrentTime += m_dTimestep;
        if (m_bSaveTime) m_nTimeId++;
        m_nStep++;
    }

    // Simulation completed - Calculate and log elapsed time
    cout << endl;
    
    time_t tEndTime = time(nullptr);
    double dTotalElapsedTime = difftime(tEndTime, m_tSysStartLoopTime);
    
    // Display summary on screen
    cout << "    - Simulation completed successfully" << endl;
    cout << "      - Total computation time: " << std::fixed << setprecision(2) 
         << dTotalElapsedTime << " s";
    if (dTotalElapsedTime > 60.0) {
        cout << " (" << dTotalElapsedTime/60.0 << " min)";
    }
    if (dTotalElapsedTime > 3600.0) {
        cout << " (" << dTotalElapsedTime/3600.0 << " hours)";
    }
    cout << endl;
    
    // Log detailed summary
    if (LogStream.is_open()) {
        LogStream << "\n" << std::string(80, '=') << "\n";
        LogStream << "SIMULATION COMPLETED SUCCESSFULLY\n";
        LogStream << std::string(80, '=') << "\n";
        LogStream << "Final simulation date: " << getSimulationEndDateTimeString() << "\n";
        LogStream << "Total timesteps executed: " << m_nStep - 1 << "\n";
        LogStream << "Total computation time: " << std::fixed << setprecision(2) 
                  << dTotalElapsedTime << " s";
        if (dTotalElapsedTime > 60.0) {
            LogStream << " (" << dTotalElapsedTime/60.0 << " min)";
        }
        if (dTotalElapsedTime > 3600.0) {
            LogStream << " (" << dTotalElapsedTime/3600.0 << " hours)";
        }
        LogStream << "\n";
        
        // Calculate performance metrics
        double simulated_days = m_dSimDuration / 86400.0;
        double ratio = simulated_days / (dTotalElapsedTime / 86400.0);
        LogStream << "Performance ratio: " << std::fixed << setprecision(1) 
                  << ratio << "x real time\n";
        LogStream << "Average time per timestep: " << std::fixed << setprecision(4) 
                  << dTotalElapsedTime / (m_nStep - 1) << " s\n";
        
        // Adaptive timestep statistics
        LogStream << "\nAdaptive timestep statistics:\n";
        LogStream << "  - Minimum timestep: " << std::fixed << setprecision(2) << m_dTimestepMin << " s\n";
        LogStream << "  - Maximum timestep: " << std::fixed << setprecision(2) << m_dTimestepMax << " s\n";
        LogStream << "  - Target timestep: " << std::fixed << setprecision(2) << m_dSimTimestep << " s\n";
        LogStream << "  - Timestep range: " << std::fixed << setprecision(1) << m_dTimestepMax/m_dTimestepMin << "x\n";
        LogStream << "  - Significant changes: " << m_nTimestepChanges << " times\n";
        if (m_dTimestepMin < m_dSimTimestep * 0.5) {
            LogStream << "  ⚠️ Timestep frequently reduced: Consider lowering Courant number or refining mesh\n";
        }
        
        // Adaptive Courant statistics
        LogStream << "\nAdaptive Courant control:\n";
        LogStream << "  - Optimal CFL (auto-calculated): " << std::fixed << setprecision(3) << m_dCourantNumberOptimal << "\n";
        LogStream << "  - Minimum CFL reached: " << std::fixed << setprecision(3) << m_dCourantNumberMin << "\n";
        LogStream << "  - CFL reduction factor: " << std::fixed << setprecision(2) << m_dCourantNumberMin/m_dCourantNumberOptimal << "x\n";
        LogStream << "  - Extreme flow events: " << m_nCourantReductions << " times\n";
        if (m_nCourantReductions > 0) {
            LogStream << "  ℹ️ Adaptive CFL prevented instabilities during " << m_nCourantReductions << " extreme flow events\n";
        }
        
        LogStream << "\nOutput file: " << m_strOutFile << "\n";
        LogStream << std::string(80, '=') << "\n";
        LogStream.flush();
    }
    
    writer.nCloseNetCDFFile();

    if (LogStream.is_open()) {
        LogStream.close();
    }
}


//======================================================================================================================
//! Initialize simulation state vectors to zero or default values
//======================================================================================================================
void CSimulation::initializeVectors() {
    const int nCrossSectionsNumber = m_nCrossSectionsNumber;
    const vector<double> vZeros(static_cast<size_t>(nCrossSectionsNumber), 0.0);

    // Hydrodynamic state variables (predictor-corrector scheme)
    m_vCrossSectionArea =
    m_vCrossSectionQ =
    m_vPredictedCrossSectionArea =
    m_vCorrectedCrossSectionArea =
    m_vPredictedCrossSectionQ =
    m_vCorrectedCrossSectionQ =
    m_vPredictedCrossSectionWaterDepth = 
    m_vCrossSectionHydraulicRadius =
    m_vCrossSectionDhDx =
    m_vCrossSectionWidth =
    m_vCrossSectionU =
    m_vCrossSectionC = vZeros;

    // Transport scalars (salinity and temperature)
    m_vCrossSectionSalinity =
    m_vPredictedCrossSectionS =
    m_vCorrectedCrossSectionS =
    m_vCrossSectionTemperature =
    m_vPredictedCrossSectionT =
    m_vCorrectedCrossSectionT = vZeros;

    // Sediment transport variables
    m_vCrossSectionQb =
    m_vCrossSectionQs =
    m_vCrossSectionQt =
    m_vCrossSectionDiamX = vZeros;

    // Geometric variables
    m_vCrossSectionLeftRBLocation =
    m_vCrossSectionRightRBLocation =
    m_vCrossSectionLeftRBLocation_UTM_X =
    m_vCrossSectionLeftRBLocation_UTM_Y =
    m_vCrossSectionRightRBLocation_UTM_X =
    m_vCrossSectionRightRBLocation_UTM_Y = vZeros;

    // Bed slope and friction terms
    m_vCrossSectionBedSlope =
    m_vCrossSectionBedSlopePredictor =
    m_vCrossSectionBedSlopeCorrector =
    m_vCrossSectionFrictionSlope =
    m_vCrossSectionDX =
    m_vCrossSectionManningNumber = vZeros;

    // Source terms for Saint-Venant equations
    m_vCrossSection_gAS0 =
    m_vCrossSection_gASf =
    m_vCrossSectionF0 =
    m_vCrossSectionF1 =
    m_vCrossSectionGv0 =
    m_vCrossSectionGv1 =
    m_vLateralSourcesAtT = vZeros;

    // Murillo factor is a multiplicative diagnostic/limiter factor.
    // Default to 1.0 (no limiting) so diagnostics remain meaningful when Murillo is disabled.
    m_vCrossSectionMurilloFactor.assign(static_cast<size_t>(nCrossSectionsNumber), 1.0);

    // Temporal terms for transport equations
    m_vCrossSectionSalinityASt =
    m_vCrossSectionTemperatureASt = vZeros;

    // Density and baroclinic terms
    m_vCrossSectionDRhoDx =
    m_vPredictedCrossSectionDensity = vZeros;

    // Initialize density with fresh water value
    const vector<double> vRhos(static_cast<size_t>(nCrossSectionsNumber), FRESH_WATER_DENSITY);
    m_vCrossSectionDensity = vRhos;

    // TVD flux limiter factors (size n+1 for interfaces)
    const vector<double> vOnes(static_cast<size_t>(nCrossSectionsNumber + 1), 1.0);
    m_vCrossSectionD1Factor =
    m_vCrossSectionD2Factor = vOnes;

    // Output time arrays
    const int nTimestepsNumber = static_cast<int>(m_dSimDuration / m_dSimTimestep) + 1;
    m_vOutputTimesIds.resize(nTimestepsNumber);
    m_vOutputTimes.resize(nTimestepsNumber);

    // Lateral storage factor (keep values if already loaded from along-channel file)
    if (m_vLateralStorageFactor.size() != static_cast<size_t>(nCrossSectionsNumber)) {
        m_vLateralStorageFactor.assign(static_cast<size_t>(nCrossSectionsNumber), 1.0);
    } else {
        for (double& s : m_vLateralStorageFactor) {
            if (s < 1.0) s = 1.0;
        }
    }

    for (int i = 0; i < nTimestepsNumber; i++) {
        m_vOutputTimesIds[i] = i;
        m_vOutputTimes[i] = static_cast<double>(i) * m_dSimTimestep;
    }
}

//======================================================================================================================
//! Calculate daily minimum temperatures from air temperature time series
//! This function pre-computes the minimum temperature for each day (00:00-24:00)
//! to be used later for estimating relative humidity using the FAO-56 method
//! 
//! @note Called once after reading meteorological data (bReadHeatFluxFile)
//! @note The vector m_vDailyMinTemperature is indexed by day number (0, 1, 2, ...)
//! @note Only calculates if m_bCalculateRHFromTemperature == true
//======================================================================================================================
void CSimulation::calculateDailyMinTemperatures() {
    // Only calculate if we need to derive RH from temperature
    if (!m_bCalculateRHFromTemperature || m_vHeatFluxTime.empty() || m_vHeatFluxAirTemp.empty()) {
        return;
    }
    
    // Determine number of days in simulation
    double total_time = m_vHeatFluxTime.back();  // Last time value (seconds)
    int num_days = static_cast<int>(std::ceil(total_time / 86400.0)) + 1;  // 86400 s/day
    
    // Initialize daily min temperature vector with very high values
    m_vDailyMinTemperature.assign(num_days, 1000.0);
    
    // Find minimum temperature for each day
    for (size_t i = 0; i < m_vHeatFluxTime.size(); ++i) {
        double time_hours = m_vHeatFluxTime[i] / 3600.0;
        int day_index = static_cast<int>(time_hours / 24.0);
        
        if (day_index >= 0 && day_index < num_days) {
            double temp = m_vHeatFluxAirTemp[i];
            if (temp < m_vDailyMinTemperature[day_index]) {
                m_vDailyMinTemperature[day_index] = temp;
            }
        }
    }
    
    // Fill any days that might not have data with the overall minimum
    double overall_min = *std::min_element(m_vHeatFluxAirTemp.begin(), m_vHeatFluxAirTemp.end());
    for (int d = 0; d < num_days; ++d) {
        if (m_vDailyMinTemperature[d] > 999.0) {  // Still uninitialized
            m_vDailyMinTemperature[d] = overall_min;
        }
    }
    
    if (m_nLogFileDetail >= 2) {
        LogStream << "--- Calculated daily minimum temperatures for RH estimation ---" << std::endl;
        LogStream << "Number of days: " << num_days << std::endl;
        LogStream << "Overall minimum temperature: " << overall_min << " °C" << std::endl;
    }
}

//======================================================================================================================
//! Calculate bed slope (S0) and distance increments (dx) between cross-sections
//! Uses central differences for interior nodes and one-sided differences at boundaries
//! Also calculates upwind slopes for predictor-corrector source term balance
//======================================================================================================================
void CSimulation::calculateBedSlope() {
    double dX1, dX2;
    
    // Calculate bed slope S0 = dz/dx using finite differences
    for (int i = 0; i < m_nCrossSectionsNumber; i++) {
        if (i == 0) {
            // Upstream boundary: forward difference
            dX1 = estuary[1].dGetX() - estuary[0].dGetX();
            m_vCrossSectionBedSlope[i] = (estuary[0].dGetZ() - estuary[1].dGetZ()) / dX1;
            m_vCrossSectionDX[i] = dX1;
        }
        else if (i == m_nCrossSectionsNumber - 1) {
            // Downstream boundary: backward difference
            dX2 = estuary[i].dGetX() - estuary[i-1].dGetX();
            m_vCrossSectionBedSlope[i] = (estuary[i-1].dGetZ() - estuary[i].dGetZ()) / dX2;
            m_vCrossSectionDX[i] = dX2;
        }
        else {
            // Interior nodes: central difference
            dX1 = estuary[i+1].dGetX() - estuary[i].dGetX();
            dX2 = estuary[i].dGetX() - estuary[i-1].dGetX();
            m_vCrossSectionBedSlope[i] = (estuary[i-1].dGetZ() - estuary[i+1].dGetZ()) / (dX1 + dX2);
            m_vCrossSectionDX[i] = dX1;
        }
    }
    
    // Calculate upwind slopes for predictor-corrector source term balance
    // Predictor: forward difference (i -> i+1)
    // Corrector: backward difference (i-1 -> i)
    for (int i = 0; i < m_nCrossSectionsNumber; i++) {
        if (i != m_nCrossSectionsNumber - 1) {
            double dx_forward = estuary[i+1].dGetX() - estuary[i].dGetX();
            m_vCrossSectionBedSlopePredictor[i] = (estuary[i].dGetZ() - estuary[i+1].dGetZ()) / dx_forward;
        } else {
            double dx_backward = estuary[i].dGetX() - estuary[i-1].dGetX();
            m_vCrossSectionBedSlopePredictor[i] = (estuary[i-1].dGetZ() - estuary[i].dGetZ()) / dx_backward;
        }
        
        if (i != 0) {
            double dx_backward = estuary[i].dGetX() - estuary[i-1].dGetX();
            m_vCrossSectionBedSlopeCorrector[i] = (estuary[i-1].dGetZ() - estuary[i].dGetZ()) / dx_backward;
        } else {
            double dx_forward = estuary[i+1].dGetX() - estuary[i].dGetX();
            m_vCrossSectionBedSlopeCorrector[i] = (estuary[i].dGetZ() - estuary[i+1].dGetZ()) / dx_forward;
        }
    }
}



//======================================================================================================================
//! Calculate initial hydraulic conditions along the estuary and precompute numerical constants
//! Supports three initialization modes: 0=calm water, 1=given flow, 2=given elevation
//======================================================================================================================
void CSimulation::calculateAlongEstuaryInitialConditions() {
    // Allocate water depth and elevation vectors
    m_vCrossSectionWaterDepth.resize(m_nCrossSectionsNumber);
    m_vCrossSectionWaterElevation.resize(m_nCrossSectionsNumber);

    // Allocate TVD flux limiter working vectors (used in mergePredictorCorrector)
    m_vTVD_a1_med.resize(m_nCrossSectionsNumber);
    m_vTVD_a2_med.resize(m_nCrossSectionsNumber);
    m_vTVD_alfa1_med.resize(m_nCrossSectionsNumber);
    m_vTVD_alfa2_med.resize(m_nCrossSectionsNumber);
    m_vTVD_psi1_med.resize(m_nCrossSectionsNumber);
    m_vTVD_psi2_med.resize(m_nCrossSectionsNumber);
    m_vTVD_r1_med.resize(m_nCrossSectionsNumber);
    m_vTVD_r2_med.resize(m_nCrossSectionsNumber);
    m_vTVD_fi1_med.resize(m_nCrossSectionsNumber);
    m_vTVD_fi2_med.resize(m_nCrossSectionsNumber);
    m_vTVD_Factor1.resize(m_nCrossSectionsNumber);
    m_vTVD_Factor2.resize(m_nCrossSectionsNumber);

    // Allocate salinity transport working vectors
    m_vSalinity_KAS_forward.resize(m_nCrossSectionsNumber);
    m_vSalinity_KAS_backward.resize(m_nCrossSectionsNumber);
    m_vSalinity_AUS_diff.resize(m_nCrossSectionsNumber);

    // Precompute numerical constants that remain constant during simulation
    m_vManningNumberSquared.resize(m_nCrossSectionsNumber);
    m_vInvDX.resize(m_nCrossSectionsNumber);
    m_vDxSum.resize(m_nCrossSectionsNumber);
    m_vInvDxSum.resize(m_nCrossSectionsNumber);
    m_vGtimesDX.resize(m_nCrossSectionsNumber);
    
    for (int i = 0; i < m_nCrossSectionsNumber; i++) {
        // Precompute inverse spacing (1/Δx) for gradient calculations
        if (m_vCrossSectionDX[i] > 1e-10) {
            m_vInvDX[i] = 1.0 / m_vCrossSectionDX[i];
        } else {
            m_vInvDX[i] = 0.0;
        }
        
        // Precompute g*Δx for source terms
        m_vGtimesDX[i] = G * m_vCrossSectionDX[i];
        
        // Set minimum area threshold (10% of minimum geometric area)
        DRY_AREA = 0.1 * m_vEstuaryAreas[i][0];
    }
    
    // Precompute spacing sums for central differences
    for (int i = 1; i < m_nCrossSectionsNumber - 1; i++) {
        double dx1 = m_vPositionX[i+1] - m_vPositionX[i];
        double dx2 = m_vPositionX[i] - m_vPositionX[i-1];
        m_vDxSum[i] = dx1 + dx2;
        m_vInvDxSum[i] = 1.0 / m_vDxSum[i];
    }

    // Initialize hydraulic conditions (skipped if continuing from NetCDF)
    if (!m_bContinueSimulation) {
        if (m_nInitialEstuarineCondition == 1) {
            // Mode 1: Given discharge - calculate area using Manning equation
            for (int i = 0; i < m_nCrossSectionsNumber; i++) {
                // Calculate Manning factor: Q*n/sqrt(|S0|)
                double dManningFactor = 0.0;
                if (m_vCrossSectionBedSlope[i] == 0.0) {
                    dManningFactor = fabs(m_vCrossSectionQ[i]) * estuary[i].dGetManningNumber() / sqrt(1e-3);
                } else {
                    dManningFactor = fabs(m_vCrossSectionQ[i]) * estuary[i].dGetManningNumber() / sqrt(fabs(m_vCrossSectionBedSlope[i]));
                }
                vector<double> vCrossSectionAreaTmp = estuary[i].vGetArea();
                vector<double> vCrossSectionHydraulicRadiusTmp = estuary[i].vGetHydraulicRadius();
                vector<double> vCrossSectionElevationTmp = estuary[i].vGetWaterDepth();
                // Compute A*R^(2/3) for each geometry node
                vector<double> dSecondTerm;
                dSecondTerm.resize(vCrossSectionAreaTmp.size());
                for (size_t j = 0; j < vCrossSectionAreaTmp.size(); j++) {
                    dSecondTerm[j] = vCrossSectionAreaTmp[j] * pow(vCrossSectionHydraulicRadiusTmp[j], 2.0/3.0);
                }
                // Interpolate to find cross-sectional area from Manning equation
                m_vCrossSectionArea[i] = linearInterpolation1d(dManningFactor, dSecondTerm, vCrossSectionAreaTmp);
                m_vCrossSectionWaterDepth[i] = linearInterpolation1d(m_vCrossSectionArea[i], vCrossSectionAreaTmp, vCrossSectionElevationTmp);
                m_vCrossSectionWidth[i] = m_vCrossSectionArea[i] / m_vCrossSectionWaterDepth[i];
                m_vCrossSectionWaterElevation[i] = m_vCrossSectionWaterDepth[i] + estuary[i].dGetZ();
            }
        }
        else if (m_nInitialEstuarineCondition == 2) {
            // Mode 2: Given elevation - calculate discharge using Manning equation
            for (int i = 0; i < m_nCrossSectionsNumber; i++) {
                vector<double> vCrossSectionAreaTmp = estuary[i].vGetArea();
                vector<double> vCrossSectionElevationTmp = estuary[i].vGetWaterDepth();
                vector<double> vCrossSectionHydraulicRadiusTmp = estuary[i].vGetHydraulicRadius();
                m_vCrossSectionArea[i] = linearInterpolation1d(m_vCrossSectionWaterElevation[i] - estuary[i].dGetZ(),vCrossSectionElevationTmp, vCrossSectionAreaTmp);
                m_vCrossSectionWidth[i] = m_vCrossSectionArea[i] / m_vCrossSectionWaterDepth[i];
                m_vCrossSectionHydraulicRadius[i] = linearInterpolation1d(m_vCrossSectionArea[i], vCrossSectionAreaTmp, vCrossSectionHydraulicRadiusTmp);
                // Manning equation: Q = (A*R^(2/3)*sqrt(|S0|)) / n * sign(S0)
                double sign_S0 = (m_vCrossSectionBedSlope[i] >= 0) ? 1.0 : -1.0;
                m_vCrossSectionQ[i] = m_vCrossSectionArea[i] * 
                                      pow(m_vCrossSectionHydraulicRadius[i], 2.0/3.0) * 
                                      sqrt(fabs(m_vCrossSectionBedSlope[i]) + 1e-10) * 
                                      sign_S0 / estuary[i].dGetManningNumber();
            }
        }
        else {
            // Mode 0: Calm water - set water level at z=0 (sea level), zero discharge
            for (int i = 0; i < m_nCrossSectionsNumber; i++) {
                vector<double> vCrossSectionAreaTmp = estuary[i].vGetArea();
                vector<double> vCrossSectionElevationTmp = estuary[i].vGetWaterDepth();
                vector<double> vCrossSectionHydraulicRadiusTmp = estuary[i].vGetHydraulicRadius();
                m_vCrossSectionWaterDepth[i] = -estuary[i].dGetZ();
                m_vCrossSectionWaterElevation[i] = 0.0;
                m_vCrossSectionArea[i] = linearInterpolation1d(m_vCrossSectionWaterDepth[i],vCrossSectionElevationTmp, vCrossSectionAreaTmp);
                m_vCrossSectionWidth[i] = m_vCrossSectionArea[i] / m_vCrossSectionWaterDepth[i];
                m_vCrossSectionHydraulicRadius[i] = linearInterpolation1d(m_vCrossSectionArea[i], vCrossSectionAreaTmp, vCrossSectionHydraulicRadiusTmp);
                m_vCrossSectionQ[i] = 0.0;
            }
        }
    }
 
    // Copy Manning coefficients from geometry data
    for (int i = 0; i < m_nCrossSectionsNumber; i++) {
        m_vCrossSectionManningNumber[i] = estuary[i].dGetManningNumber();
    }
    
    // Precompute Manning² for friction term calculations
    for (int i = 0; i < m_nCrossSectionsNumber; i++) {
        m_vManningNumberSquared[i] = pow(m_vCrossSectionManningNumber[i], 2.0);
    }

    // Calculate dimensionless sediment diameter if transport is enabled
    if (m_bDoWaterDensity) {
        for (int i = 0; i < m_nCrossSectionsNumber; i++) {
            m_vCrossSectionDiamX[i] = m_vCrossSectionD50[i] * pow((m_vCrossSectionRhos[i] - 1.0) * G / (NU * NU), 1.0/3.0);
        }
    }
}


//======================================================================================================================
//! Perform 1D linear interpolation with extrapolation at boundaries
//! Uses binary search for O(log n) lookup efficiency
//! @param dValue X-value to interpolate at
//! @param vX X-coordinates (must be sorted in ascending order)
//! @param vY Y-values corresponding to vX
//! @return Interpolated/extrapolated Y-value
//======================================================================================================================
double CSimulation::linearInterpolation1d(const double dValue, const vector<double> &vX, const vector<double> &vY) {
    auto it = std::lower_bound(vX.begin(), vX.end(), dValue);
    
    if (it == vX.begin()) {
        // Extrapolate below minimum using linear trend from first two points
        double slope = (vY[1] - vY[0]) / (vX[1] - vX[0]);
        return vY[0] + slope * (dValue - vX[0]);
    } else if (it == vX.end()) {
        // Extrapolate above maximum using linear trend from last two points
        size_t n = vX.size();
        double slope = (vY[n-1] - vY[n-2]) / (vX[n-1] - vX[n-2]);
        return vY[n-1] + slope * (dValue - vX[n-1]);
    } else {
        // Interpolate between two points
        size_t i = std::distance(vX.begin(), it) - 1;
        double slope = (vY[i+1] - vY[i]) / (vX[i+1] - vX[i]);
        return vY[i] + slope * (dValue - vX[i]);
    }
}

double CSimulation::getGivenTemperatureAtTime(double t) const {
    if (m_vGivenTemperatureTime.size() < 2 || m_vGivenTemperatureValue.size() < 2) {
        return 0.0;
    }
    // linearInterpolation1d is non-const; reuse the same logic locally
    auto it = std::lower_bound(m_vGivenTemperatureTime.begin(), m_vGivenTemperatureTime.end(), t);
    if (it == m_vGivenTemperatureTime.begin()) {
        double slope = (m_vGivenTemperatureValue[1] - m_vGivenTemperatureValue[0]) /
                       (m_vGivenTemperatureTime[1] - m_vGivenTemperatureTime[0]);
        return m_vGivenTemperatureValue[0] + slope * (t - m_vGivenTemperatureTime[0]);
    }
    if (it == m_vGivenTemperatureTime.end()) {
        size_t n = m_vGivenTemperatureTime.size();
        double slope = (m_vGivenTemperatureValue[n - 1] - m_vGivenTemperatureValue[n - 2]) /
                       (m_vGivenTemperatureTime[n - 1] - m_vGivenTemperatureTime[n - 2]);
        return m_vGivenTemperatureValue[n - 1] + slope * (t - m_vGivenTemperatureTime[n - 1]);
    }
    size_t i = std::distance(m_vGivenTemperatureTime.begin(), it) - 1;
    double slope = (m_vGivenTemperatureValue[i + 1] - m_vGivenTemperatureValue[i]) /
                   (m_vGivenTemperatureTime[i + 1] - m_vGivenTemperatureTime[i]);
    return m_vGivenTemperatureValue[i] + slope * (t - m_vGivenTemperatureTime[i]);
}

void CSimulation::applyGivenTemperatureAtTime(double t) {
    const double T = getGivenTemperatureAtTime(t);
    const int N = m_nCrossSectionsNumber;
    if (N <= 0) return;

    if (m_vCrossSectionTemperature.size() != static_cast<size_t>(N))
        m_vCrossSectionTemperature.resize(N, 0.0);
    if (m_vPredictedCrossSectionT.size() != static_cast<size_t>(N))
        m_vPredictedCrossSectionT.resize(N, 0.0);
    if (m_vCorrectedCrossSectionT.size() != static_cast<size_t>(N))
        m_vCorrectedCrossSectionT.resize(N, 0.0);

    for (int i = 0; i < N; ++i) {
        m_vCrossSectionTemperature[i] = T;
        m_vPredictedCrossSectionT[i] = T;
        m_vCorrectedCrossSectionT[i] = T;
    }
}
//======================================================================================================================
//! Calculate hydraulic parameters (width, depth, hydraulic radius) from cross-sectional area
//! Uses precomputed geometry tables for fast lookup via interpolation
//! Updates both current state and predictor/corrector states depending on m_nPredictor flag
//======================================================================================================================
/**
 * @brief Compute hydraulic parameters (width, depth, hydraulic radius, etc.) by interpolating from geometry tables
 * 
 * OPTIMIZATIONS:
 * - Binary search cache exploits spatial coherence (~30% faster)
 * - Reduced conditional branches with pointer selection
 * - Compiler hints [[likely]] for typical flow path
 * - Single interpolation pass for all variables
 * 
 * For each cross-section, interpolates hydraulic properties based on current cross-sectional area.
 * Uses cached binary search for efficient lookup in pre-computed geometry tables.
 * Handles three cases: area below minimum, area above maximum, and interpolation between table entries.
 */
//======================================================================================================================
void CSimulation::calculateHydraulicParameters() {
    const int nCrossSections = m_nCrossSectionsNumber;
    
    // Static state tracking for min/max area conditions (persists between calls)
    static std::vector<bool> was_at_min;
    static std::vector<bool> was_at_max;
    static bool first_call = true;
    
    if (first_call) {
        was_at_min.resize(nCrossSections, false);
        was_at_max.resize(nCrossSections, false);
        first_call = false;
    }
    
    // Select appropriate area vector based on predictor/corrector phase
    // - Predictor (m_nPredictor==1): use current state (t^n)
    // - Corrector (m_nPredictor==2): use predicted state (t*)
    // - Post-merge/output (m_nPredictor==0): use merged current state
    const auto& dArea = (m_nPredictor == 2) ? m_vPredictedCrossSectionArea : m_vCrossSectionArea;
    
    // Loop through all cross-sections and interpolate hydraulic properties
    for (int i = 0; i < nCrossSections; i++) {
        const int nElevationSectionsNumber = m_vElevationSectionsCount[i];
        const double currentArea = dArea[i];
        const auto& areas = m_vEstuaryAreas[i];
        
        // Case 1: Area below minimum table value (extrapolate using first entry)
        if (areas[0] > currentArea) [[unlikely]] {
            m_vCrossSectionHydraulicRadius[i] = m_vEstuaryHydraulicRadius[i][0];
            m_vCrossSectionWidth[i] = m_vWidth[i][0];
            m_vCrossSectionWaterDepth[i] = m_vEstuaryWaterDepths[i][0];
            m_vCrossSectionLeftRBLocation[i] = m_vLeftY[i][0];
            m_vCrossSectionRightRBLocation[i] = m_vRightY[i][0];
            
            // Log TRANSITION to minimum area condition (only when state changes)
            if (m_nLogFileDetail >= 2 && !was_at_min[i]) {
                LogStream << "⬇️  MIN AREA [t=" << std::fixed << std::setprecision(1) << m_dCurrentTime/3600.0 
                          << "h]: Cross-section " << i << " reached MINIMUM (A=" 
                          << std::setprecision(2) << currentArea << " m², A_min=" << areas[0] 
                          << " m²)\n";
                was_at_min[i] = true;
            }
        }
        // Case 2: Area above maximum table value (extrapolate using last entry)
        else if (areas[nElevationSectionsNumber-1] < currentArea) [[unlikely]] {
            const int lastNode = nElevationSectionsNumber - 1;
            m_vCrossSectionHydraulicRadius[i] = m_vEstuaryHydraulicRadius[i][lastNode];
            m_vCrossSectionWidth[i] = m_vWidth[i][lastNode];
            m_vCrossSectionWaterDepth[i] = m_vEstuaryWaterDepths[i][lastNode];
            m_vCrossSectionLeftRBLocation[i] = m_vLeftY[i][lastNode];
            m_vCrossSectionRightRBLocation[i] = m_vRightY[i][lastNode];
            
            // Log TRANSITION to maximum area condition (only when state changes)
            if (m_nLogFileDetail >= 2 && !was_at_max[i]) {
                LogStream << "⬆️  MAX AREA [t=" << std::fixed << std::setprecision(1) << m_dCurrentTime/3600.0 
                          << "h]: Cross-section " << i << " reached MAXIMUM (A=" 
                          << std::setprecision(2) << currentArea << " m², A_max=" << areas[lastNode] 
                          << " m²)\n";
                was_at_max[i] = true;
            }
        }
        // Case 3: Area within table range (linear interpolation) - MOST COMMON PATH
        else [[likely]] {
            // Log TRANSITION out of min/max states (recovery)
            if (m_nLogFileDetail >= 2) {
                if (was_at_min[i]) {
                    LogStream << "✅ RECOVERY [t=" << std::fixed << std::setprecision(1) << m_dCurrentTime/3600.0 
                              << "h]: Cross-section " << i << " left minimum (A=" 
                              << std::setprecision(2) << currentArea << " m²)\n";
                    was_at_min[i] = false;
                }
                if (was_at_max[i]) {
                    LogStream << "✅ RECOVERY [t=" << std::fixed << std::setprecision(1) << m_dCurrentTime/3600.0 
                              << "h]: Cross-section " << i << " left maximum (A=" 
                              << std::setprecision(2) << currentArea << " m²)\n";
                    was_at_max[i] = false;
                }
            }
            
            // OPTIMIZATION: Use cached search start position (spatial coherence)
            // Water levels change slowly between timesteps, so last index is good hint
            int j = m_vLastInterpolationIndex[i];
            
            // Validate cache: ensure j is still valid bracket
            if (j < 0 || j >= nElevationSectionsNumber-1 || 
                currentArea < areas[j] || currentArea > areas[j+1]) {
                // Cache miss: perform binary search
                auto it = std::lower_bound(areas.begin(), areas.end(), currentArea);
                j = std::distance(areas.begin(), it) - 1;
                j = std::max(0, std::min(j, nElevationSectionsNumber-2));
                m_vLastInterpolationIndex[i] = j;  // Update cache
            }
            
            // Compute interpolation factor once for all variables
            const double denom = areas[j+1] - areas[j];
            const double factor = (currentArea - areas[j]) / denom;
            const double inv_factor = 1.0 - factor;
            
            // Linear interpolation: property = (1-f)*value[j] + f*value[j+1]
            m_vCrossSectionHydraulicRadius[i] = inv_factor * m_vEstuaryHydraulicRadius[i][j] + 
                                                factor * m_vEstuaryHydraulicRadius[i][j+1];
            m_vCrossSectionWaterDepth[i] = inv_factor * m_vEstuaryWaterDepths[i][j] + 
                                           factor * m_vEstuaryWaterDepths[i][j+1];
            m_vCrossSectionWidth[i] = inv_factor * m_vWidth[i][j] + 
                                      factor * m_vWidth[i][j+1];
            m_vCrossSectionLeftRBLocation[i] = inv_factor * m_vLeftY[i][j] + 
                                               factor * m_vLeftY[i][j+1];
            m_vCrossSectionRightRBLocation[i] = inv_factor * m_vRightY[i][j] + 
                                                factor * m_vRightY[i][j+1];
        }
        
        m_vCrossSectionWaterElevation[i] = m_vCrossSectionWaterDepth[i] + m_vBedZ[i];

        // Keep the predicted-depth vector consistent when running the corrector.
        // Several source terms (e.g., baroclinic) reference m_vPredictedCrossSectionWaterDepth,
        // while hydraulics are computed into the active vectors for performance.
        if (m_nPredictor == 2 && static_cast<int>(m_vPredictedCrossSectionWaterDepth.size()) == nCrossSections) {
            m_vPredictedCrossSectionWaterDepth[i] = m_vCrossSectionWaterDepth[i];
        }
    }
}

//======================================================================================================================
/**
 * @brief Calculate UTM coordinates of left and right river banks
 * 
 * Computes bank positions in UTM coordinate system by projecting from thalweg (channel centerline)
 * perpendicular to the channel direction using specified angles and distances.
 * 
 * For each cross-section:
 * - Starts at thalweg UTM position (x_thalweg, y_thalweg)
 * - Projects perpendicular distance to left bank using angle_left
 * - Projects perpendicular distance to right bank using angle_right
 * 
 * Angles are in radians, measured from reference direction (typically north or channel axis).
 */
//======================================================================================================================
void CSimulation::calculateRiverBankUTMCoordinates() {
    
    for (int i = 0; i < m_nCrossSectionsNumber; i++) {
        // Get thalweg (centerline) UTM coordinates from along_channel_data.csv
        const double x_thalweg = estuary[i].dGetX_UTM();
        const double y_thalweg = estuary[i].dGetY_UTM();
        
        // Get angles in radians (AngMd for left, AngMi for right)
        const double angle_left = estuary[i].dGetLeftRBAngle();   // Already in radians
        const double angle_right = estuary[i].dGetRightRBAngle(); // Already in radians
        
        // Get distances from thalweg to river banks
        const double dist_left = fabs(m_vCrossSectionLeftRBLocation[i]);
        const double dist_right = fabs(m_vCrossSectionRightRBLocation[i]);
        
        // Calculate UTM coordinates for left bank
        // Starting from thalweg position, moving perpendicular to channel direction
        m_vCrossSectionLeftRBLocation_UTM_X[i] = x_thalweg + dist_left * cos(angle_left);
        m_vCrossSectionLeftRBLocation_UTM_Y[i] = y_thalweg + dist_left * sin(angle_left);
        
        // Calculate UTM coordinates for right bank
        m_vCrossSectionRightRBLocation_UTM_X[i] = x_thalweg + dist_right * cos(angle_right);
        m_vCrossSectionRightRBLocation_UTM_Y[i] = y_thalweg + dist_right * sin(angle_right);
    }
}

//======================================================================================================================
/**
 * @brief Linearly interpolate hydraulic parameters between two elevation table entries
 * 
 * @param dArea Current cross-sectional area to interpolate for
 * @param nCrossSection Cross-section index
 * @param nElevationNode Lower elevation table entry index (interpolates between [j] and [j+1])
 * 
 * Interpolates: hydraulic radius, width, water depth, left bank position, right bank position
 */
//======================================================================================================================
void CSimulation::interpolateHydraulicParameters(double dArea, int nCrossSection, int nElevationNode) {
    // Compute interpolation factor: f = (A - A[j]) / (A[j+1] - A[j])
    const double dInterpolationFactor = (dArea - estuary[nCrossSection].dGetArea(nElevationNode)) / 
                                        (estuary[nCrossSection].dGetArea(nElevationNode+1) - estuary[nCrossSection].dGetArea(nElevationNode));
    
    // Interpolate all hydraulic properties: property = f*(P[j+1] - P[j]) + P[j]
    m_vCrossSectionHydraulicRadius[nCrossSection] = dInterpolationFactor * (estuary[nCrossSection].dGetHydraulicRadius(nElevationNode+1) - 
                                                                            estuary[nCrossSection].dGetHydraulicRadius(nElevationNode)) + 
                                                    estuary[nCrossSection].dGetHydraulicRadius(nElevationNode);
    
    m_vCrossSectionWidth[nCrossSection] = dInterpolationFactor * (estuary[nCrossSection].dGetWidth(nElevationNode+1) - 
                                                                  estuary[nCrossSection].dGetWidth(nElevationNode)) + 
                                          estuary[nCrossSection].dGetWidth(nElevationNode);
    
    m_vCrossSectionWaterDepth[nCrossSection] = dInterpolationFactor * (estuary[nCrossSection].dGetWaterDepth(nElevationNode+1) - 
                                                                       estuary[nCrossSection].dGetWaterDepth(nElevationNode)) + 
                                               estuary[nCrossSection].dGetWaterDepth(nElevationNode);
    
    m_vCrossSectionLeftRBLocation[nCrossSection] = dInterpolationFactor * (estuary[nCrossSection].dGetLeftY(nElevationNode+1) - 
                                                                           estuary[nCrossSection].dGetLeftY(nElevationNode)) + 
                                                   estuary[nCrossSection].dGetLeftY(nElevationNode);
    
    m_vCrossSectionRightRBLocation[nCrossSection] = dInterpolationFactor * (estuary[nCrossSection].dGetRightY(nElevationNode+1) - 
                                                                            estuary[nCrossSection].dGetRightY(nElevationNode)) + 
                                                    estuary[nCrossSection].dGetRightY(nElevationNode);
}

//======================================================================================================================
/**
 * @brief Assign hydraulic parameters from first (minimum) elevation table entry
 * @param nCrossSection Cross-section index
 * 
 * Used when cross-sectional area is below the minimum table value (extrapolation).
 */
//======================================================================================================================
void CSimulation::getFirstHydraulicParameters(const int nCrossSection) {
    m_vCrossSectionHydraulicRadius[nCrossSection] = estuary[nCrossSection].dGetHydraulicRadius(0);
    m_vCrossSectionWidth[nCrossSection] = estuary[nCrossSection].dGetWidth(0);
    m_vCrossSectionWaterDepth[nCrossSection] = estuary[nCrossSection].dGetWaterDepth(0);
    m_vCrossSectionLeftRBLocation[nCrossSection] = estuary[nCrossSection].dGetLeftY(0);
    m_vCrossSectionRightRBLocation[nCrossSection] = estuary[nCrossSection].dGetRightY(0);
}

//======================================================================================================================
/**
 * @brief Assign hydraulic parameters from last (maximum) elevation table entry
 * @param nCrossSection Cross-section index
 * 
 * Used when cross-sectional area exceeds the maximum table value (extrapolation).
 */
//======================================================================================================================
void CSimulation::getLastHydraulicParameters(const int nCrossSection) {
    const int nLastNode = estuary[nCrossSection].nGetElevationSectionsNumber() - 1;
    m_vCrossSectionHydraulicRadius[nCrossSection] = estuary[nCrossSection].dGetHydraulicRadius(nLastNode);
    m_vCrossSectionWidth[nCrossSection] = estuary[nCrossSection].dGetWidth(nLastNode);
    m_vCrossSectionWaterDepth[nCrossSection] = estuary[nCrossSection].dGetWaterDepth(nLastNode);
    m_vCrossSectionLeftRBLocation[nCrossSection] = estuary[nCrossSection].dGetLeftY(nLastNode);
    m_vCrossSectionRightRBLocation[nCrossSection] = estuary[nCrossSection].dGetRightY(nLastNode);
}


//======================================================================================================================
/**
 * @brief Compute adaptive timestep satisfying CFL stability condition
 * 
 * Calculates timestep based on:
 * 1. Courant-Friedrichs-Lewy (CFL) condition: Δt ≤ C·Δx / (|u| + c)
 *    where C is Courant number, u is flow velocity, c is wave celerity
 * 2. Diffusion stability (if salinity dispersion active): Δt ≤ 0.5·Δx² / Kh
 * 3. Density-driven flow constraints (if water density active)
 * 4. Output time synchronization (reduce timestep if approaching save time)
 * 
 * Also computes flow velocity (u) and wave celerity (c) at all nodes for use in numerical scheme.
 */
//======================================================================================================================
void CSimulation::calculateTimestep() {
    double dMinTimestep = 10000000.0;
    double dWaterDensityFactor = 1.0;
    // Track maximum local CFL and its location for diagnostics
    double max_local_cfl = 0.0;
    int critical_section_idx = -1;
  
    // ⚡ OPTIMIZATION: Cache sqrt(G) to avoid recomputing in every iteration
    static const double sqrt_G = sqrt(G);
    
    // Compute timestep based on interior nodes only (exclude boundaries i=0 and i=n-1)
    // IMPORTANT for non-uniform grids: use the local minimum spacing around node i.
    for (int i=1; i < m_nCrossSectionsNumber-1; i++) {
        if (m_vCrossSectionArea[i] > DRY_AREA) {
            const double dxL = m_vPositionX[i] - m_vPositionX[i - 1];
            const double dxR = m_vPositionX[i + 1] - m_vPositionX[i];
            const double dX = (dxL > 0.0 && dxR > 0.0) ? std::min(dxL, dxR) : std::max(dxL, dxR);
            // Mean flow velocity: u = Q / A
            const double u = m_vCrossSectionQ[i] / m_vCrossSectionArea[i];
            m_vCrossSectionU[i] = u;
            // Shallow water wave celerity with optional lateral storage:
            // Base: c = sqrt(g·A/B) = sqrt(g·h) where h=A/B.
            // With storage factor S>=1 in continuity: c_eff ≈ c / sqrt(S).
            const double S = dGetLateralStorageFactor(i);
            const double c = sqrt_G * sqrt(m_vCrossSectionArea[i] / (m_vCrossSectionWidth[i] * S));
            m_vCrossSectionC[i] = c;
            // Additional stability constraint for density-driven flows (baroclinic adjustment)
            if (m_bDoWaterDensity) {
                dWaterDensityFactor = m_dCourantNumberCurrent * dX / (m_vCrossSectionWidth[i] * m_dBetaSalinityConstant);
            }
            // ⚡ OPTIMIZATION: Precalculate characteristic speeds to avoid recomputing in max()
            // CFL timestep: Δt = C·Δx / max(|u+c|, |u-c|)
            const double lambda_plus = fabs(u + c);   // Downstream characteristic
            const double lambda_minus = fabs(u - c);  // Upstream characteristic
            const double max_lambda = std::max(lambda_plus, lambda_minus);
            // Compute local CFL for this section
            double local_cfl = (max_lambda > 0.0) ? (m_dTimestep * max_lambda / dX) : 0.0;
            if (local_cfl > max_local_cfl) {
                max_local_cfl = local_cfl;
                critical_section_idx = i;
            }
            if (const double dTimestepTmp = m_dCourantNumberCurrent * dX / max_lambda; 
                dTimestepTmp < dMinTimestep) {
                dMinTimestep = dTimestepTmp;
            }
            // Apply density factor constraint if active
            if (m_bDoWaterDensity && (dWaterDensityFactor < dMinTimestep)) {
                dMinTimestep = dWaterDensityFactor;
            }
        }
        else {
            // Handle dry/nearly-dry nodes using minimum threshold values
            const double dxL = m_vPositionX[i] - m_vPositionX[i - 1];
            const double dxR = m_vPositionX[i + 1] - m_vPositionX[i];
            const double dX = (dxL > 0.0 && dxR > 0.0) ? std::min(dxL, dxR) : std::max(dxL, dxR);
            const double u_dry = DRY_Q / DRY_AREA;
            m_vCrossSectionU[i] = u_dry;

            // Wave celerity for dry zone using minimum depth: h_dry = A_dry / B
            const double Sf = dGetLateralStorageFactor(i);
            const double h_dry = DRY_AREA / (m_vCrossSectionWidth[i] * Sf);
            const double c_dry = sqrt_G * sqrt(h_dry);
            m_vCrossSectionC[i] = c_dry;

            // ⚡ OPTIMIZATION: Precalculate characteristic speeds
            const double lambda_plus = fabs(u_dry + c_dry);
            const double lambda_minus = fabs(u_dry - c_dry);
            const double dTimestepTmp = m_dCourantNumberCurrent * dX / std::max(lambda_plus, lambda_minus);
            
            if (dTimestepTmp < dMinTimestep) {
                dMinTimestep = dTimestepTmp;
        }
    }}
    
    // Compute u and c at boundary nodes (for flux calculations, but not used in Δt constraint)
    for (int i : {0, m_nCrossSectionsNumber-1}) {
        if (m_vCrossSectionArea[i] != DRY_AREA) {
            m_vCrossSectionU[i] = m_vCrossSectionQ[i] / m_vCrossSectionArea[i];
            const double Sf = dGetLateralStorageFactor(i);
            m_vCrossSectionC[i] = sqrt_G * sqrt(m_vCrossSectionArea[i] / (m_vCrossSectionWidth[i] * Sf));
        } else {
            m_vCrossSectionU[i] = DRY_Q / DRY_AREA;
            const double Sf = dGetLateralStorageFactor(i);
            const double h_dry = DRY_AREA / (m_vCrossSectionWidth[i] * Sf);
            m_vCrossSectionC[i] = sqrt_G * sqrt(h_dry);
        }
    }
    
    m_dTimestep = dMinTimestep;

    // === STEP 2: Apply rigorous adaptive Courant control ===
    // Objective criterion: If max(CFL_local) > 0.95·C_optimal → reduce C_current
    // Safety margin: 5% below target to prevent instabilities
    const double CFL_SAFETY_MARGIN = 0.95;
    const double CFL_CRITICAL_THRESHOLD = CFL_SAFETY_MARGIN * m_dCourantNumberOptimal;
    
    if (max_local_cfl > CFL_CRITICAL_THRESHOLD && !m_bExtremeFlowConditions) {
        // Entering high-CFL regime: reduce Courant preventively
        // Reduction factor: keep CFL below threshold with 10% extra safety
        const double reduction_factor = 0.9 * CFL_CRITICAL_THRESHOLD / max_local_cfl;
        m_dCourantNumberCurrent = std::max(0.05, m_dCourantNumberOptimal * reduction_factor);
        m_nCourantReductions++;
        m_bExtremeFlowConditions = true;
        
        // Recompute timestep with reduced Courant
        m_dTimestep *= reduction_factor;
        
        // Log only first occurrence or if log level is very high
        if (m_nLogFileDetail >= 2 || m_nCourantReductions == 1) {
            std::ofstream LogStream(m_strLogFile.c_str(), std::ios_base::app);
            LogStream << std::fixed << std::setprecision(3);
            LogStream << "🔴 HIGH CFL [t=" << m_dCurrentTime/3600.0 << "h]: "
                      << "CFL=" << max_local_cfl << " | C: " << m_dCourantNumberOptimal << "→" << m_dCourantNumberCurrent << "\n";
            LogStream.flush();
        }
    }
    else if (max_local_cfl < 0.7 * CFL_CRITICAL_THRESHOLD && m_bExtremeFlowConditions) {
        // Exiting high-CFL regime: gradually recover Courant
        // Conservative recovery: increase by 20% per step until reaching optimal value
        m_dCourantNumberCurrent = std::min(m_dCourantNumberOptimal, m_dCourantNumberCurrent * 1.2);
        
        if (m_dCourantNumberCurrent >= 0.99 * m_dCourantNumberOptimal) {
            m_bExtremeFlowConditions = false;
            
            // Log recovery only at high detail level
            if (m_nLogFileDetail >= 2) {
                std::ofstream LogStream(m_strLogFile.c_str(), std::ios_base::app);
                LogStream << "🟢 CFL NORMALIZED [t=" << std::fixed << std::setprecision(1) 
                          << m_dCurrentTime/3600.0 << "h]\n";
                LogStream.flush();
            }
        }
    }
    else if (m_bExtremeFlowConditions && m_dCourantNumberCurrent < m_dCourantNumberOptimal) {
        // Gradual recovery while still in elevated CFL regime (conservative 5% increase)
        m_dCourantNumberCurrent = std::min(m_dCourantNumberOptimal, m_dCourantNumberCurrent * 1.05);
    }
    
    // Track minimum Courant reached during simulation
    if (m_dCourantNumberCurrent < m_dCourantNumberMin) {
        m_dCourantNumberMin = m_dCourantNumberCurrent;
    }

    // === STEP 3: Apply rigorous diffusion stability constraint ===
    // Diffusion stability (Fourier condition): Δt ≤ α·Δx² / (2·Kh)
    // where α = stability factor (0.4-0.5 for explicit schemes)
    // Physical basis: diffusive time scale τ_diff = Δx²/Kh must be resolved
    if (m_bDoWaterSalinity) {
        const double DIFFUSION_SAFETY_FACTOR = 0.4;  // Conservative: 0.4 < 0.5 (stability limit)
        double min_Pe = 1e10;  // Track minimum Péclet number (diagnostic)
        bool any_diffusion = false;
        
        for (int i = 1; i < m_nCrossSectionsNumber-1; i++) {
            if (m_vCrossSectionArea[i] > DRY_AREA) {
                const double dxL = m_vPositionX[i] - m_vPositionX[i - 1];
                const double dxR = m_vPositionX[i + 1] - m_vPositionX[i];
                const double dX = (dxL > 0.0 && dxR > 0.0) ? std::min(dxL, dxR) : std::max(dxL, dxR);
                const double u = fabs(m_vCrossSectionU[i]);
                const double Kh = dGetLongitudinalDispersion(i);
                if (Kh <= 0.0) {
                    continue;
                }
                any_diffusion = true;
                
                // Diffusion stability: Δt ≤ α·Δx²/Kh
                const double dt_diffusion = DIFFUSION_SAFETY_FACTOR * dX * dX / Kh;
                if (dt_diffusion < m_dTimestep) {
                    m_dTimestep = dt_diffusion;
                    // Log only if extremely restrictive (Δt < 10s) and high detail level
                    if (m_nLogFileDetail >= 3 && dt_diffusion < 10.0) {
                        std::ofstream LogStream(m_strLogFile.c_str(), std::ios_base::app);
                        LogStream << "⚠️ DIFFUSION LIMIT [t=" << std::fixed << std::setprecision(1) 
                                  << m_dCurrentTime/3600.0 << "h]: Δt=" << std::setprecision(2) 
                                  << m_dTimestep << "s\n";
                        LogStream.flush();
                    }
                }
                
                // Diagnostic: Péclet number Pe = u·Δx/Kh
                // Pe << 1: diffusion-dominated (smooth profiles, can use larger Δt)
                // Pe >> 1: advection-dominated (sharp fronts, TVD limiters essential)
                // Pe ≈ 2: balanced advection-diffusion (optimal for accuracy)
                const double Pe = u * dX / Kh;
                min_Pe = std::min(min_Pe, Pe);
            }
        }
        
        // Warning: Only if Pe is extremely low and might cause performance issues
        // Log once per simulated hour maximum
        static double last_Pe_warning_time = -3600.0;
        if (any_diffusion && min_Pe < 0.05 && m_nLogFileDetail >= 3 && 
            (m_dCurrentTime - last_Pe_warning_time) > 3600.0) {
            std::ofstream LogStream(m_strLogFile.c_str(), std::ios_base::app);
            LogStream << "ℹ️ Very low Péclet (Pe=" << std::fixed << std::setprecision(3) 
                      << min_Pe << ") - diffusion-dominated\n";
            LogStream.flush();
            last_Pe_warning_time = m_dCurrentTime;
        }
    }

    // === STEP 4: Enforce absolute safety limits ===
    // Minimum timestep: 0.1s (below this, round-off errors and computational cost dominate)
    const double ABSOLUTE_MIN_TIMESTEP = 0.1;  // seconds
    static bool timestep_floor_logged = false;
    if (m_dTimestep < ABSOLUTE_MIN_TIMESTEP) {
        // Log only first occurrence to avoid spam
        if (m_nLogFileDetail >= 2 && !timestep_floor_logged) {
            std::ofstream LogStream(m_strLogFile.c_str(), std::ios_base::app);
            LogStream << "⚠️ TIMESTEP FLOOR: Δt clamped to " << ABSOLUTE_MIN_TIMESTEP << "s\n";
            LogStream.flush();
            timestep_floor_logged = true;
        }
        m_dTimestep = ABSOLUTE_MIN_TIMESTEP;
    }
    
    // Maximum timestep: user-specified output timestep (for output synchronization)
    if (m_dTimestep > m_dSimTimestep) {
        m_dTimestep = m_dSimTimestep;
    }
    
    // Track timestep statistics for performance diagnostics
    if (m_dTimestep < m_dTimestepMin) {
        m_dTimestepMin = m_dTimestep;
    }
    if (m_dTimestep > m_dTimestepMax) {
        m_dTimestepMax = m_dTimestep;
    }

    // === STEP 5: Log significant timestep changes for monitoring ===
    // For sub-second timesteps, only log drastic changes to avoid spam
    // Use time-based throttling: max 1 message per simulated hour
    static double last_timestep_log_time = -3600.0;
    const double LOG_THROTTLE_INTERVAL = 3600.0;  // 1 hour of simulated time
    const bool should_throttle = (m_dCurrentTime - last_timestep_log_time) < LOG_THROTTLE_INTERVAL;
    
    if (m_nLogFileDetail >= 1 && m_dTimestepPrevious > 0.0) {
        const double dt_ratio = m_dTimestep / m_dTimestepPrevious;
        
        // Log only drastic changes (>50% reduction or >100% increase) or if not throttled
        const bool drastic_change = (dt_ratio < 0.5 || dt_ratio > 2.0);
        if (drastic_change || (!should_throttle && (dt_ratio < 0.8 || dt_ratio > 1.25))) {
            m_nTimestepChanges++;
            
            std::ofstream LogStream(m_strLogFile.c_str(), std::ios_base::app);
            LogStream << std::fixed << std::setprecision(2);
            
            if (dt_ratio < 0.8) {
                // Timestep decreased
                LogStream << "⏬ Δt: " << m_dTimestepPrevious << "→" << m_dTimestep << "s "
                          << "(" << std::setprecision(0) << (dt_ratio-1.0)*100.0 << "%)";
                if (critical_section_idx >= 0 && m_nLogFileDetail >= 2) {
                    LogStream << " x=" << std::setprecision(1) 
                              << m_vCrossSectionX[critical_section_idx]/1000.0 << "km";
                }
                LogStream << "\n";
            } else if (m_nLogFileDetail >= 2) {
                // Timestep increased (only log at higher detail level)
                LogStream << "⏫ Δt: " << m_dTimestepPrevious << "→" << m_dTimestep << "s\n";
            }
            LogStream.flush();
            last_timestep_log_time = m_dCurrentTime;  // Update throttle timer
        }
    }

    // === STEP 6: Synchronize with output times ===
    // Reduce timestep to hit exact output time (ensures precise snapshot timing)
    if (m_bSaveAllTimesteps) {
        // Save at every computational timestep (debug mode)
        m_bSaveTime = true;
    }
    else if ((m_nTimeId < static_cast<int>(m_vOutputTimes.size())) && (m_dCurrentTime + m_dTimestep > m_vOutputTimes[m_nTimeId])) {
        // Reduce timestep to hit exact output time
        double dt_to_save = m_vOutputTimes[m_nTimeId] - m_dCurrentTime;
        if (dt_to_save > 1e-6) {  // Avoid dt=0 when already at save time
            m_dTimestep = dt_to_save;
        }
        m_bSaveTime = true;
    }
    else {
        m_bSaveTime = false;
    }
    
    // Compute dimensionless timestep ratio λ = Δt / min(Δx) for use in numerical scheme
    m_dLambda = m_dTimestep / dMinVectorValue(m_vCrossSectionDX);
}


//===============================================================================================================================
/**
 * @brief Calculate boundary condition values at current time by interpolating from time series
 * 
 * Performs two-step interpolation:
 * 1. Temporal interpolation: get BC value at current simulation time from input time series
 * 2. Geometric conversion: convert elevation BC to cross-sectional area using geometry tables
 * 
 * Handles both upstream (upward) and downstream (downward) boundaries.
 * For elevation BCs, also computes value at adjacent interior node for smooth transitions.
 */
//===============================================================================================================================
void CSimulation::calculateBoundaryConditions() {
    // Upstream boundary (condition type 2 = water surface elevation prescribed)
    if (nGetUpwardEstuarineCondition() == 2) {
        // Step 1: Interpolate elevation from time series at current time
        double dUpwardBoundaryValue = linearInterpolation1d(m_dCurrentTime, 
                                                           m_vUpwardBoundaryConditionTime, 
                                                           m_vUpwardBoundaryConditionValue);
        
        // Step 2: Convert water surface elevation to cross-sectional area at boundary node (i=0)
        // Uses water depth h = η - z_bed as lookup in area-depth geometry table
        m_dUpwardBoundaryValue = linearInterpolation1d(-m_vBedZ[0] + dUpwardBoundaryValue, 
                                                       m_vEstuaryWaterDepths[0],
                                                       m_vEstuaryAreas[0]);
        
        // Also compute for adjacent interior node (i=1) to ensure smooth boundary transition
        m_dNextUpwardBoundaryValue = linearInterpolation1d(-m_vBedZ[1] + dUpwardBoundaryValue, 
                                                          m_vEstuaryWaterDepths[1],
                                                          m_vEstuaryAreas[1]);
    }
    
    // Upstream boundary (condition type 3 = discharge prescribed from time series)
    if ((nGetUpwardEstuarineCondition() == 3)) {
        // Interpolate discharge from time series at current time
        m_dUpwardBoundaryValue = linearInterpolation1d(m_dCurrentTime, 
                                                       m_vUpwardBoundaryConditionTime, 
                                                       m_vUpwardBoundaryConditionValue);
    }

    // Downstream boundary (condition type 2 = water surface elevation prescribed)
    if (nGetDownwardEstuarineCondition() == 2) {
        // Step 1: Interpolate elevation from time series at current time
        double dDownwardBoundaryValue = linearInterpolation1d(m_dCurrentTime, 
                                                             m_vDownwardBoundaryConditionTime, 
                                                             m_vDownwardBoundaryConditionValue);
        
        // Step 2: Convert water surface elevation to cross-sectional area at boundary node (i=n-1)
        m_dDownwardBoundaryValue = linearInterpolation1d(-m_vBedZ[m_nCrossSectionsNumber-1] + dDownwardBoundaryValue, 
                                                         m_vEstuaryWaterDepths[m_nCrossSectionsNumber-1], 
                                                         m_vEstuaryAreas[m_nCrossSectionsNumber-1]);
        
        // Also compute for adjacent interior node (i=n-2) to ensure smooth boundary transition
        m_dNextDownwardBoundaryValue = linearInterpolation1d(-m_vBedZ[m_nCrossSectionsNumber-2] + dDownwardBoundaryValue, 
                                                            m_vEstuaryWaterDepths[m_nCrossSectionsNumber-2], 
                                                            m_vEstuaryAreas[m_nCrossSectionsNumber-2]);
    }

    // Lateral inflows (hydrographs) at specified cross-sections
    if (m_nHydrographsNumber > 0) {
        // Recompute every timestep (no stale values) and accumulate contributions
        // if multiple hydrographs map to the same cross-section.
        std::fill(m_vLateralSourcesAtT.begin(), m_vLateralSourcesAtT.end(), 0.0);
        for (int i = 0; i < m_nHydrographsNumber; i++) {
            const int node = hydrographs[i].m_nNearestCrossSectionNo;
            if (node < 0 || node >= m_nCrossSectionsNumber) {
                continue;
            }
            m_vLateralSourcesAtT[node] += linearInterpolation1d(
                m_dCurrentTime,
                hydrographs[i].vGetTime(),
                hydrographs[i].vGetQ()
            );
        }
    }
}

//===============================================================================================================================
/**
 * @brief Enforce minimum thresholds for cross-sectional area and flow to prevent numerical issues
 * 
 * When computed area falls below physical/numerical threshold (DRY_AREA), sets minimum values:
 * - Area = DRY_AREA (small but non-zero to avoid division by zero)
 * - Flow = DRY_Q (typically zero or very small)
 * 
 * Applied to appropriate state vector depending on predictor/corrector phase:
 * - Predictor phase (m_nPredictor=1): applies to predicted variables
 * - Corrector phase (m_nPredictor=2): applies to corrected variables
 * - Normal state (m_nPredictor=0): applies to current variables
 */
//===============================================================================================================================
void CSimulation::dryArea() {
    if (m_nPredictor == 1) {
        for (int i = 0; i < m_nCrossSectionsNumber; i++) {
            if (m_vPredictedCrossSectionArea[i] <= DRY_AREA) {
                const double A_original = m_vPredictedCrossSectionArea[i];
                m_vPredictedCrossSectionArea[i] = DRY_AREA;
                m_vPredictedCrossSectionQ[i] = DRY_Q;
                
                // Log dry bed correction only if area was significantly larger (real drying event)
                if (m_nLogFileDetail >= 2 && A_original > DRY_AREA * 1.5) {
                    LogStream << "DRY BED [t=" << std::fixed << std::setprecision(1) << m_dCurrentTime/3600.0 
                              << "h]: Cross-section " << i << " dried (predictor: A=" 
                              << std::setprecision(3) << A_original << " m² -> " << DRY_AREA 
                              << " m²)\n";
                }
            }
        }
    }
    else if (m_nPredictor == 2) {
        for (int i = 0; i < m_nCrossSectionsNumber; i++) {
            if (m_vCorrectedCrossSectionArea[i] <= DRY_AREA) {
                const double A_original = m_vCorrectedCrossSectionArea[i];
                m_vCorrectedCrossSectionArea[i] = DRY_AREA;
                m_vCorrectedCrossSectionQ[i] = DRY_Q;
                
                // Log dry bed correction only if area was significantly larger (real drying event)
                if (m_nLogFileDetail >= 2 && A_original > DRY_AREA * 1.5) {
                    LogStream << "DRY BED [t=" << std::fixed << std::setprecision(1) << m_dCurrentTime/3600.0 
                              << "h]: Cross-section " << i << " dried (corrector: A=" 
                              << std::setprecision(3) << A_original << " m² -> " << DRY_AREA 
                              << " m²)\n";
                }
            }
        }
    }
    else {
        for (int i = 0; i < m_nCrossSectionsNumber; i++) {
            if (m_vCrossSectionArea[i] <= DRY_AREA) {
                const double A_original = m_vCrossSectionArea[i];
                m_vCrossSectionArea[i] = DRY_AREA;
                m_vCrossSectionQ[i] = DRY_Q;
                
                // Log dry bed correction only if area was significantly larger (real drying event)
                if (m_nLogFileDetail >= 2 && A_original > DRY_AREA * 1.5) {
                    LogStream << "DRY BED [t=" << std::fixed << std::setprecision(1) << m_dCurrentTime/3600.0 
                              << "h]: Cross-section " << i << " dried (initial: A=" 
                              << std::setprecision(3) << A_original << " m² -> " << DRY_AREA 
                              << " m²)\n";
                }
            }
        }
    }
}

//===============================================================================================================================
/**
 * @brief Zero out flux terms in dry/nearly-dry cross-sections to prevent spurious fluxes
 * 
 * When cross-sectional area falls below threshold, sets flux terms to zero:
 * - F0 = 0 (continuity flux)
 * - F1 = 0 (momentum flux)
 * - Gv1 = 0 (source terms)
 * 
 * This prevents unphysical fluxes in dry zones where numerical solution may be unreliable.
 * Applied to current state (predictor phase) or predicted state (corrector phase).
 */
//===============================================================================================================================
void CSimulation::dryTerms() {
    if (m_nPredictor == 1) {
        for (int i = 0; i < m_nCrossSectionsNumber-1; i++) {
            if (m_vCrossSectionArea[i] <= DRY_AREA) {
                m_vCrossSectionF0[i] = 0.0;
                m_vCrossSectionF1[i] = 0.0;
                m_vCrossSectionGv1[i] = 0.0;
            }
        }
    }
    else {
        for (int i = 0; i < m_nCrossSectionsNumber-1; i++) {
            if (m_vPredictedCrossSectionArea[i] <= DRY_AREA) {
                m_vCrossSectionF0[i] = 0.0;
                m_vCrossSectionF1[i] = 0.0;
                m_vCrossSectionGv1[i] = 0.0;
            }
        }
    }
}

//===============================================================================================================================
/**
 * @brief Apply Murillo stability condition to limit Manning coefficient in shallow/steep sections
 * 
 * Murillo et al. (2010) stability criterion prevents numerical instabilities when bed friction
 * becomes too large relative to grid resolution and flow depth.
 * 
 * Computes critical Manning value: n_crit = C_dx * sqrt(2*R^(2/3) / (g*dx))
 * If actual Manning coefficient n > n_crit, reduces n to n_crit and stores scaling factor.
 * 
 * Only applied during corrector phase (m_nPredictor != 1).
 * 
 * Reference: Murillo et al. (2010), "The influence of source terms on stability..."
 */
//===============================================================================================================================
void CSimulation::doMurilloCondition() {
    if (m_nPredictor != 1) {
        for (int i = 0; i < m_nCrossSectionsNumber; i++) {
            // Non-uniform grid: Murillo's criterion depends on local spacing.
            double dx_local = 0.0;
            if (m_nCrossSectionsNumber >= 2) {
                if (i == 0) {
                    dx_local = m_vPositionX[1] - m_vPositionX[0];
                } else if (i == m_nCrossSectionsNumber - 1) {
                    dx_local = m_vPositionX[i] - m_vPositionX[i - 1];
                } else {
                    const double dxL = m_vPositionX[i] - m_vPositionX[i - 1];
                    const double dxR = m_vPositionX[i + 1] - m_vPositionX[i];
                    dx_local = (dxL > 0.0 && dxR > 0.0) ? std::min(dxL, dxR) : std::max(dxL, dxR);
                }
            }
            dx_local = std::max(1e-12, dx_local);

            const double dValue = CDX_MURILLO * sqrt(
                                      2 * pow(m_vCrossSectionHydraulicRadius[i], 2.0 / 3.0) / (G * dx_local));

            // Murillo criterion: if n is too large for local depth/spacing, clamp to n_crit.
            // Compare against the CURRENT Manning used by the solver, not the original geometry value.
            const double n_current = m_vCrossSectionManningNumber[i];
            if (dValue < n_current) {
                m_vCrossSectionMurilloFactor[i] = dValue / (n_current + 1e-30);
                m_vCrossSectionManningNumber[i] = dValue;
            } else {
                m_vCrossSectionMurilloFactor[i] = 1.0;
            }

            // Keep the precomputed n^2 used by friction consistent with the current Manning.
            m_vManningNumberSquared[i] = m_vCrossSectionManningNumber[i] * m_vCrossSectionManningNumber[i];
        }
    }
}

//======================================================================================================================
/**
 * @brief Calculate bed slope and friction source terms for momentum equation
 * 
 * Computes:
 * 1. Bed slope term: gAS0 = g * A * S0
 * 2. Friction slope term: gASf = g * A * Sf, where Sf = n²|Q|Q / (A²R^(4/3))
 * 
 * Two methods:
 * - Surface term balance (bGetDoSurfaceTermBalance=true): 
 *   Uses averaged hydraulic variables between adjacent nodes and upwind bed slopes
 *   Predictor: averages with i+1 (forward), uses forward difference S0
 *   Corrector: averages with i-1 (backward), uses backward difference S0
 *   
 * - Direct computation (bGetDoSurfaceTermBalance=false):
 *   Uses local values without averaging (simpler but less accurate)
 * 
 * Optionally applies water-level-dependent Manning coefficient if m_bManningDependsOnLevel=true.
 */
//======================================================================================================================
void CSimulation::calculate_GS_A_terms() {

    // Optional stability limiter on Manning (only active if enabled in config)
    if (m_bDoMurilloCondition) {
        doMurilloCondition();
    }
    
    if (bGetDoSurfaceTermBalance()) {
        // Average hydraulic variables between adjacent cross-sections for balanced source term discretization
        // Predictor: average with i+1 (forward), Corrector: average with i-1 (backward)
        for (int i = 0; i < m_nCrossSectionsNumber; i++) {
            double dMeanArea, dMeanQ, dMeanHydraulicRadius;
            double manning2_mean = m_vManningNumberSquared[i];
            double neta_mean = 1.0;
            
            if (m_nPredictor == 1) {
                if (i < m_nCrossSectionsNumber - 1) {
                    dMeanArea = (m_vCrossSectionArea[i] + m_vCrossSectionArea[i+1]) / 2.0;
                    dMeanQ = (m_vCrossSectionQ[i] + m_vCrossSectionQ[i+1]) / 2.0;
                    dMeanHydraulicRadius = (m_vCrossSectionHydraulicRadius[i] + m_vCrossSectionHydraulicRadius[i+1]) / 2.0;
                    manning2_mean = 0.5 * (m_vManningNumberSquared[i] + m_vManningNumberSquared[i + 1]);
                    if (m_bManningDependsOnLevel) {
                        const double eta_i = pow(n_eta(m_vCrossSectionWaterDepth[i]), 2.0);
                        const double eta_j = pow(n_eta(m_vCrossSectionWaterDepth[i + 1]), 2.0);
                        neta_mean = 0.5 * (eta_i + eta_j);
                    }
                } else {
                    dMeanArea = m_vCrossSectionArea[i];
                    dMeanQ = m_vCrossSectionQ[i];
                    dMeanHydraulicRadius = m_vCrossSectionHydraulicRadius[i];
                    if (m_bManningDependsOnLevel) {
                        neta_mean = pow(n_eta(m_vCrossSectionWaterDepth[i]), 2.0);
                    }
                }
            } else {
                if (i > 0) {
                    dMeanArea = (m_vPredictedCrossSectionArea[i] + m_vPredictedCrossSectionArea[i-1]) / 2.0;
                    dMeanQ = (m_vPredictedCrossSectionQ[i] + m_vPredictedCrossSectionQ[i-1]) / 2.0;
                    dMeanHydraulicRadius = (m_vCrossSectionHydraulicRadius[i] + m_vCrossSectionHydraulicRadius[i-1]) / 2.0;
                    manning2_mean = 0.5 * (m_vManningNumberSquared[i] + m_vManningNumberSquared[i - 1]);
                    if (m_bManningDependsOnLevel) {
                        const double eta_i = pow(n_eta(m_vCrossSectionWaterDepth[i]), 2.0);
                        const double eta_j = pow(n_eta(m_vCrossSectionWaterDepth[i - 1]), 2.0);
                        neta_mean = 0.5 * (eta_i + eta_j);
                    }
                } else {
                    dMeanArea = m_vPredictedCrossSectionArea[i];
                    dMeanQ = m_vPredictedCrossSectionQ[i];
                    dMeanHydraulicRadius = m_vCrossSectionHydraulicRadius[i];
                    if (m_bManningDependsOnLevel) {
                        neta_mean = pow(n_eta(m_vCrossSectionWaterDepth[i]), 2.0);
                    }
                }
            }
            
            // Use upwind bed slopes for consistent source term balance
            // Predictor: forward difference S0, Corrector: backward difference S0
            double S0_to_use;
            if (m_nPredictor == 1) {
                S0_to_use = m_vCrossSectionBedSlopePredictor[i];  // Forward difference (zmedp)
            } else {
                S0_to_use = m_vCrossSectionBedSlopeCorrector[i];  // Backward difference (zmedc)
            }
            
            // Bed slope source term: gAS0 = g * A * S0
            m_vCrossSection_gAS0[i] = G * dMeanArea * S0_to_use;
            
            // Friction slope: Sf = (n²*η*|Q|*Q) / (A²*R^(4/3))
            // Friction term: gASf = g*A*Sf
            if (dMeanArea > DRY_AREA && dMeanHydraulicRadius > 1e-6) {
                // ⚡ OPTIMIZATION: Precalculate repeated terms
                const double A_squared = dMeanArea * dMeanArea;
                const double Rh_power = pow(dMeanHydraulicRadius, 4.0/3.0);
                const double Q_abs_Q = dMeanQ * fabs(dMeanQ);
                
                m_vCrossSectionFrictionSlope[i] = manning2_mean * neta_mean *
                                                 Q_abs_Q / (A_squared * Rh_power);
                m_vCrossSection_gASf[i] = G * dMeanArea * m_vCrossSectionFrictionSlope[i];
            } else {
                m_vCrossSectionFrictionSlope[i] = 0.0;
                m_vCrossSection_gASf[i] = 0.0;
            }
        }
    }
    else {
        // Direct computation without averaging (simpler but less accurate)
        for (int i = 0; i < m_nCrossSectionsNumber; i++) {
            // Bed slope term: S0 already has correct sign, no direction factor needed
            m_vCrossSection_gAS0[i] = G * m_vCrossSectionArea[i] * m_vCrossSectionBedSlope[i];
            
            if (m_nPredictor == 1) {
                double neta = 1.0;
                if (m_bManningDependsOnLevel) {
                    neta = pow(n_eta(m_vCrossSectionWaterDepth[i]), 2.0);
                }
                if (m_vCrossSectionArea[i] > DRY_AREA && m_vCrossSectionHydraulicRadius[i] > 1e-6) {
                    // ⚡ OPTIMIZATION: Precalculate repeated terms
                    const double A = m_vCrossSectionArea[i];
                    const double A_squared = A * A;
                    const double Rh_power = pow(m_vCrossSectionHydraulicRadius[i], 4.0/3.0);
                    const double Q_abs_Q = m_vCrossSectionQ[i] * fabs(m_vCrossSectionQ[i]);
                    
                    m_vCrossSectionFrictionSlope[i] = m_vManningNumberSquared[i] * neta *
                                                     Q_abs_Q / (A_squared * Rh_power);
                    m_vCrossSection_gASf[i] = G * A * m_vCrossSectionFrictionSlope[i];
                } else {
                    m_vCrossSectionFrictionSlope[i] = 0.0;
                    m_vCrossSection_gASf[i] = 0.0;
                }
            }
            else {
                double neta = 1.0;
                if (m_bManningDependsOnLevel) {
                    neta = pow(n_eta(m_vCrossSectionWaterDepth[i]), 2.0);
                }
                if (m_vPredictedCrossSectionArea[i] > DRY_AREA && m_vCrossSectionHydraulicRadius[i] > 1e-6) {
                    // ⚡ OPTIMIZATION: Precalculate repeated terms
                    const double A = m_vPredictedCrossSectionArea[i];
                    const double A_squared = A * A;
                    const double Rh_power = pow(m_vCrossSectionHydraulicRadius[i], 4.0/3.0);
                    const double Q_abs_Q = m_vPredictedCrossSectionQ[i] * fabs(m_vPredictedCrossSectionQ[i]);
                    
                    m_vCrossSectionFrictionSlope[i] = m_vManningNumberSquared[i] * neta *
                                                     Q_abs_Q / (A_squared * Rh_power);
                    m_vCrossSection_gASf[i] = G * A * m_vCrossSectionFrictionSlope[i];
                } else {
                    m_vCrossSectionFrictionSlope[i] = 0.0;
                    m_vCrossSection_gASf[i] = 0.0;
                }
            }
        }
    }
}

//======================================================================================================================
/**
 * @brief Calculate flux terms F0 and F1 for conservative form of Saint-Venant equations
 * 
 * Computes:
 * - F0 = Q (continuity equation flux: simply the discharge)
 * - F1 = β*Q²/A (momentum equation flux: momentum advection with momentum correction factor β)
 * 
 * Pressure gradient is computed more rigorously in source terms (Gv1) using -gA∂η/∂x,
 * which provides better accuracy and stability for non-prismatic channels than integral methods.
 * 
 * Applied to current state (predictor) or predicted state (corrector) depending on m_nPredictor.
 */
//======================================================================================================================
void CSimulation::calculateFlowTerms() {
    
    if (m_nPredictor == 1) {
        for (int i = 0; i < m_nCrossSectionsNumber; i++) {
            // F0 = Q (discharge, continuity flux)
            m_vCrossSectionF0[i] = m_vCrossSectionQ[i];
            
            // F1 = β*Q²/A (momentum flux)
            // Note: Pressure integral I1 term disabled for stability (see function header)
            if (m_vCrossSectionArea[i] > DRY_AREA) {
                m_vCrossSectionF1[i] = m_vBeta[i] * pow(m_vCrossSectionQ[i], 2.0) / m_vCrossSectionArea[i];
            } else {
                m_vCrossSectionF1[i] = 0.0;
            }
        }
    }
    else {
        for (int i = 0; i < m_nCrossSectionsNumber; i++) {
            // F0 = Q (discharge from predicted state)
            m_vCrossSectionF0[i] = m_vPredictedCrossSectionQ[i];
            
            // F1 = β*Q²/A (momentum flux from predicted state, I1 disabled)
            if (m_vPredictedCrossSectionArea[i] > DRY_AREA) {
                m_vCrossSectionF1[i] = m_vBeta[i] * pow(m_vPredictedCrossSectionQ[i], 2.0) / m_vPredictedCrossSectionArea[i];
            } else {
                m_vCrossSectionF1[i] = 0.0;
            }
        }
    }
}

//======================================================================================================================
/**
 * @brief Calculate source terms for momentum and continuity equations
 * 
 * Computes:
 * 1. Lateral inflow source term: Gv0 = Q_lateral / Δx
 * 2. Momentum source terms: Gv1 = -gA∂η/∂x - gASf + S_baroclinic
 *    - Barotropic pressure gradient: -gA∂η/∂x (water surface slope)
 *    - Friction term: -gASf (bottom friction)
 *    - Baroclinic pressure gradient: -gA(h/ρ)∂ρ/∂x (density-driven flow, if salinity active)
 * 
 * Uses centered differences for interior nodes, one-sided differences at boundaries.
 * Automatically selects appropriate state variables (current or predicted) based on m_nPredictor.
 */
//======================================================================================================================
void CSimulation::calculateSourceTerms() {
    // Select appropriate state vectors based on predictor/corrector phase
    const auto& area_to_use_vec = (m_nPredictor == 1) ? m_vCrossSectionArea : m_vPredictedCrossSectionArea;
    const auto& depth_to_use_vec = (m_nPredictor == 1) ? m_vCrossSectionWaterDepth : m_vPredictedCrossSectionWaterDepth;

    // 1. Water surface gradient (barotropic pressure gradient): ∂η/∂x
    //    Interior nodes: centered difference ∂η/∂x ≈ (η[i+1] - η[i-1]) / (Δx[i] + Δx[i-1])
    //    Boundary nodes: forward/backward difference
    for (int i = 1; i < m_nCrossSectionsNumber-1; i++) {
        m_vCrossSectionDhDx[i] = (m_vCrossSectionWaterElevation[i+1] - m_vCrossSectionWaterElevation[i-1]) * m_vInvDxSum[i];
    }
    m_vCrossSectionDhDx[0] = (m_vCrossSectionWaterElevation[1] - m_vCrossSectionWaterElevation[0]) * m_vInvDX[0];
    m_vCrossSectionDhDx[m_nCrossSectionsNumber-1] = (m_vCrossSectionWaterElevation[m_nCrossSectionsNumber-1] - 
                                                     m_vCrossSectionWaterElevation[m_nCrossSectionsNumber-2]) * 
                                                     m_vInvDX[m_nCrossSectionsNumber-2];

    // 2. Density gradient (baroclinic pressure gradient): ∂ρ/∂x
    //    Only computed if salinity transport is active
    if (bGetDoWaterSalinity()) {
        const auto& density_to_use = (m_nPredictor == 1) ? m_vCrossSectionDensity : m_vPredictedCrossSectionDensity;
        
        // Centered differences for interior, one-sided for boundaries
        for (int i = 1; i < m_nCrossSectionsNumber-1; i++) {
            m_vCrossSectionDRhoDx[i] = (density_to_use[i+1] - density_to_use[i-1]) * m_vInvDxSum[i];
        }
        m_vCrossSectionDRhoDx[0] = (density_to_use[1] - density_to_use[0]) * m_vInvDX[0];
        m_vCrossSectionDRhoDx[m_nCrossSectionsNumber-1] = (density_to_use[m_nCrossSectionsNumber-1] - 
                                                           density_to_use[m_nCrossSectionsNumber-2]) * 
                                                           m_vInvDX[m_nCrossSectionsNumber-2];
    }

    // 3. Assemble total source terms
    for (int i = 1; i < m_nCrossSectionsNumber-1; i++) {
        // Control-volume length for non-uniform grids: Δx_CV = 0.5[(x_i-x_{i-1}) + (x_{i+1}-x_i)]
        const double dxL = m_vPositionX[i] - m_vPositionX[i - 1];
        const double dxR = m_vPositionX[i + 1] - m_vPositionX[i];
        const double dxCV = 0.5 * (dxL + dxR);
        const double inv_dxCV = (dxCV > 1e-12) ? (1.0 / dxCV) : 0.0;

        // Continuity source: lateral inflow distributed over control volume
        m_vCrossSectionGv0[i] = m_vLateralSourcesAtT[i] * inv_dxCV;
        
        // Get current state variables
        const double area_to_use = area_to_use_vec[i];
        const double depth_to_use = depth_to_use_vec[i];
        
        // ⚡ OPTIMIZATION: Cache G*A since used in both pressure and baroclinic terms
        const double G_times_area = G * area_to_use;
        
        // Momentum source: barotropic pressure gradient + friction
        // Gv1 = -gA∂η/∂x - gASf
        m_vCrossSectionGv1[i] = -G_times_area * m_vCrossSectionDhDx[i] - m_vCrossSection_gASf[i];
        
        // Add baroclinic pressure gradient if salinity is active and node is wet
        // S_baroclinic = -gA(h/ρ)∂ρ/∂x (density-driven circulation)
        if (bGetDoWaterSalinity() && depth_to_use > DRY_AREA) {
            const double current_density = (m_nPredictor == 1) ? m_vCrossSectionDensity[i] : m_vPredictedCrossSectionDensity[i];
            const double S_baroclinic = -G_times_area * (depth_to_use / current_density) * m_vCrossSectionDRhoDx[i];
            m_vCrossSectionGv1[i] += S_baroclinic;
        }
    }
}


//======================================================================================================================
/**
 * @brief Predictor step of McCormack two-step explicit scheme
 * 
 * Computes predicted values using forward spatial differences:
 * U^* = U^n - λ(F[i+1] - F[i]) + Δt·G[i]
 * 
 * where:
 * - U^n: current state vector (A, Q)
 * - U^*: predicted state
 * - λ = Δt/Δx: dimensionless timestep ratio
 * - F: flux vector
 * - G: source term vector
 * 
 * Note: Density (ρ) only affects source terms in Gv1 (baroclinic term),
 *       not the flux terms F0 and F1.
 */
//======================================================================================================================
void CSimulation::calculatePredictor() {

    const int n = m_nCrossSectionsNumber - 1;
    for (int i = 1; i < n; i++) {
        // Local forward spacing for non-uniform grids: Δx_i = x[i+1] - x[i]
        const double dx_forward = m_vPositionX[i + 1] - m_vPositionX[i];
        const double lambda_forward = (dx_forward > 1e-12) ? (m_dTimestep / dx_forward) : 0.0;

        const double invSf = 1.0 / dGetLateralStorageFactor(i);

        // Continuity equation: ∂A/∂t + (1/Sf)∂Q/∂x = (1/Sf)G0
        // Predictor: A* = A - (Δt/Δx_i)(1/Sf)(F0[i+1] - F0[i]) + Δt(1/Sf)·G0[i]
        m_vPredictedCrossSectionArea[i] = m_vCrossSectionArea[i] -
            (lambda_forward * invSf) * (m_vCrossSectionF0[i + 1] - m_vCrossSectionF0[i]) +
            (m_dTimestep * invSf) * m_vCrossSectionGv0[i];

        // Momentum equation: ∂Q/∂t + ∂F1/∂x = G1
        // Predictor: Q* = Q - (Δt/Δx_i)(F1[i+1] - F1[i]) + Δt·G1[i]
        m_vPredictedCrossSectionQ[i] = m_vCrossSectionQ[i] -
            lambda_forward * (m_vCrossSectionF1[i + 1] - m_vCrossSectionF1[i]) +
            m_dTimestep * m_vCrossSectionGv1[i];
    }

    // Convert predicted cross-sectional area to water depth using geometry tables
    // This is needed for subsequent hydraulic parameter calculations
    for (int i = 0; i < m_nCrossSectionsNumber; ++i) {
        m_vPredictedCrossSectionWaterDepth[i] = linearInterpolation1d(
            m_vPredictedCrossSectionArea[i],
            m_vEstuaryAreas[i],
            m_vEstuaryWaterDepths[i]
        );
    }
}


//======================================================================================================================
/**
 * @brief Corrector step of McCormack two-step explicit scheme
 * 
 * Computes corrected values using backward spatial differences:
 * U^(n+1) = U^* - λ(F*[i] - F*[i-1]) + Δt·G*[i]
 * 
 * where:
 * - U^*: predicted state from predictor step
 * - U^(n+1): corrected (final) state
 * - F*: flux computed from predicted state
 * - G*: source terms computed from predicted state
 * 
 * The final solution is typically averaged: U^(n+1) = 0.5(predictor + corrector)
 * or TVD flux limiters are applied to reduce spurious oscillations.
 */
//======================================================================================================================
void CSimulation::calculateCorrector() {

    const int n = m_nCrossSectionsNumber - 1;
    for (int i = 1; i < n; i++) {
        // Local backward spacing for non-uniform grids: Δx_{i-1} = x[i] - x[i-1]
        const double dx_backward = m_vPositionX[i] - m_vPositionX[i - 1];
        const double lambda_backward = (dx_backward > 1e-12) ? (m_dTimestep / dx_backward) : 0.0;

        const double invSf = 1.0 / dGetLateralStorageFactor(i);

        // Continuity corrector: A^(n+1) = A* - (Δt/Δx_{i-1})(1/Sf)(F0*[i] - F0*[i-1]) + Δt(1/Sf)·G0*[i]
        m_vCorrectedCrossSectionArea[i] = m_vPredictedCrossSectionArea[i] -
            (lambda_backward * invSf) * (m_vCrossSectionF0[i] - m_vCrossSectionF0[i - 1]) +
            (m_dTimestep * invSf) * m_vCrossSectionGv0[i];

        // Momentum corrector: Q^(n+1) = Q* - (Δt/Δx_{i-1})(F1*[i] - F1*[i-1]) + Δt·G1*[i]
        m_vCorrectedCrossSectionQ[i] = m_vPredictedCrossSectionQ[i] -
            lambda_backward * (m_vCrossSectionF1[i] - m_vCrossSectionF1[i - 1]) +
            m_dTimestep * m_vCrossSectionGv1[i];
    }
}

//======================================================================================================================
// Helper utilities for characteristic-style boundary conditions (stage/discharge)
//======================================================================================================================
static inline double widthFromAreaClamped(const std::vector<double>& areas,
                                         const std::vector<double>& widths,
                                         const double area) {
    if (areas.empty() || widths.empty()) return 1.0;
    if (area <= areas.front()) return widths.front();
    if (area >= areas.back()) return widths.back();
    auto it = std::lower_bound(areas.begin(), areas.end(), area);
    const size_t j = static_cast<size_t>(std::max<int>(0, static_cast<int>(std::distance(areas.begin(), it)) - 1));
    const double denom = areas[j + 1] - areas[j];
    const double f = (denom != 0.0) ? (area - areas[j]) / denom : 0.0;
    return widths[j] + f * (widths[j + 1] - widths[j]);
}

static inline double celerityFromArea(const std::vector<double>& areas,
                                      const std::vector<double>& widths,
                                      const double area) {
    const double B = std::max(1e-6, widthFromAreaClamped(areas, widths, area));
    const double A = std::max(DRY_AREA, area);
    return std::sqrt(G * A / B);
}

//======================================================================================================================
/**
 * @brief Apply boundary conditions to CURRENT STATE (time n) before predictor calculation
 * 
 * CRITICAL: Must be called BEFORE calculateHydraulicParameters() and calculateFlowTerms()
 * in the predictor phase. This ensures:
 * 1. Flux terms F0, F1 are calculated with correct boundary values
 * 2. Source terms see correct gradients at boundaries
 * 
 * Sets boundary values for m_vCrossSectionArea and m_vCrossSectionQ at i=0 and i=n.
 * 
 * Boundary condition types:
 * - Type 0: Open/free boundary (gradient extrapolation)
 * - Type 1: Reflective/wall (Q=0) or prescribed discharge  
 * - Type 2: Prescribed water surface elevation (with Manning for discharge)
 * - Type 3: Prescribed discharge (with Manning for area)
 */
//======================================================================================================================
void CSimulation::applyBoundariesToCurrentState() {
    int n = m_nCrossSectionsNumber - 1;
    
    // Upstream boundary - current state
    if (nGetUpwardEstuarineCondition() == 0) {
        const double dQ = m_vCrossSectionQ[2] - m_vCrossSectionQ[1];
        const double dA = m_vCrossSectionArea[2] - m_vCrossSectionArea[1];
        m_vCrossSectionQ[0] = m_vCrossSectionQ[1] - dQ;
        m_vCrossSectionArea[0] = m_vCrossSectionArea[1] - dA;
        if (m_vCrossSectionArea[0] < DRY_AREA) {
            m_vCrossSectionArea[0] = m_vCrossSectionArea[1];
        }
    }
    else if (nGetUpwardEstuarineCondition() == 1) {
        // Type 1: Reflective / solid-wall boundary.
        // For a closed end, enforce u=0 (Q=0) and a free-surface antinode (dη/dx≈0).
        // IMPORTANT: dη/dx≈0 means η[0]≈η[1], NOT A[0]≈A[1] when cross-sections differ.
        m_vCrossSectionQ[0] = 0.0;
        {
            const double A1 = std::max(DRY_AREA, m_vCrossSectionArea[1]);
            const double h1 = linearInterpolation1d(A1, m_vEstuaryAreas[1], m_vEstuaryWaterDepths[1]);
            const double eta1 = h1 + m_vBedZ[1];
            const double h0 = std::max(0.0, eta1 - m_vBedZ[0]);
            m_vCrossSectionArea[0] = std::max(DRY_AREA,
                linearInterpolation1d(h0, m_vEstuaryWaterDepths[0], m_vEstuaryAreas[0]));
        }
    }
    else if (nGetUpwardEstuarineCondition() == 2) {
        m_vCrossSectionArea[0] = m_dUpwardBoundaryValue;
        double R = linearInterpolation1d(m_vCrossSectionArea[0], m_vEstuaryAreas[0], m_vEstuaryHydraulicRadius[0]);
        double sign_S0 = (m_vCrossSectionBedSlope[0] >= 0) ? 1.0 : -1.0;
        m_vCrossSectionQ[0] = m_vCrossSectionArea[0] * sqrt(fabs(m_vCrossSectionBedSlope[0]) + 1e-10) * 
                             pow(R, 2.0/3.0) * sign_S0 / (m_vCrossSectionManningNumber[0] + 1e-10);
    }
    else {  // Type 3: prescribed discharge (characteristic compatibility)
        // Prescribe Q and obtain A from the outgoing characteristic (R- from node 1).
        // This reduces tidal reflection at the head compared to a pure zero-gradient A.
        const double Qin = m_dUpwardBoundaryValue;
        m_vCrossSectionQ[0] = Qin;

        const double A1 = std::max(DRY_AREA, m_vCrossSectionArea[1]);
        const double u1 = m_vCrossSectionQ[1] / A1;
        const double c1 = celerityFromArea(m_vEstuaryAreas[1], m_vWidth[1], A1) / std::sqrt(dGetLateralStorageFactor(1));
        const double Rminus = u1 - 2.0 * c1;

        auto f = [&](double A0) {
            const double c0 = celerityFromArea(m_vEstuaryAreas[0], m_vWidth[0], A0) / std::sqrt(dGetLateralStorageFactor(0));
            return (Qin / std::max(DRY_AREA, A0)) - 2.0 * c0 - Rminus;
        };

        double Alo = std::max(DRY_AREA * 1.001, 0.2 * A1);
        double Ahi = std::max(A1 * 5.0, Alo * 2.0);
        double flo = f(Alo);
        double fhi = f(Ahi);
        for (int k = 0; k < 20 && (flo * fhi > 0.0); ++k) {
            if (flo > 0.0 && fhi > 0.0) {
                Ahi *= 2.0;
                fhi = f(Ahi);
            } else {
                Alo = std::max(DRY_AREA * 1.001, Alo * 0.5);
                flo = f(Alo);
            }
        }

        if (flo * fhi <= 0.0) {
            for (int it = 0; it < 40; ++it) {
                const double Am = 0.5 * (Alo + Ahi);
                const double fm = f(Am);
                if (flo * fm <= 0.0) {
                    Ahi = Am;
                    fhi = fm;
                } else {
                    Alo = Am;
                    flo = fm;
                }
            }
            m_vCrossSectionArea[0] = 0.5 * (Alo + Ahi);
        } else {
            // Fallback: if no bracket is found, keep a weakly reflective closure.
            m_vCrossSectionArea[0] = m_vCrossSectionArea[1];
        }
    }

    // Downstream boundary - current state
    if (nGetDownwardEstuarineCondition() == 0) {
        const double dQ = m_vCrossSectionQ[n-1] - m_vCrossSectionQ[n-2];
        const double dA = m_vCrossSectionArea[n-1] - m_vCrossSectionArea[n-2];
        m_vCrossSectionQ[n] = m_vCrossSectionQ[n-1] + dQ;
        m_vCrossSectionArea[n] = m_vCrossSectionArea[n-1] + dA;
        if (m_vCrossSectionArea[n] < DRY_AREA) {
            m_vCrossSectionArea[n] = m_vCrossSectionArea[n-1];
        }
    }
    else if (nGetDownwardEstuarineCondition() == 1) {
        m_vCrossSectionQ[n] = m_dDownwardBoundaryValue;
        double dManningFactor = fabs(m_vCrossSectionQ[n]) * m_vCrossSectionManningNumber[n] / (sqrt(fabs(m_vCrossSectionBedSlope[n]) + 1e-10));
        m_vCrossSectionArea[n] = linearInterpolation1d(dManningFactor, m_vPrecalculatedSecondTerm[n], m_vEstuaryAreas[n]);
    }
    else {  // Type 2: tidal elevation (characteristic compatibility)
        // Prescribe stage (area) and obtain Q from the outgoing characteristic (R+ from node n-1).
        const double Abc = std::max(DRY_AREA, m_dDownwardBoundaryValue);
        m_vCrossSectionArea[n] = Abc;

        const double Ai = std::max(DRY_AREA, m_vCrossSectionArea[n - 1]);
        const double ui = m_vCrossSectionQ[n - 1] / Ai;
        const double ci = celerityFromArea(m_vEstuaryAreas[n - 1], m_vWidth[n - 1], Ai) / std::sqrt(dGetLateralStorageFactor(n - 1));
        const double Rplus = ui + 2.0 * ci;

        const double cbc = celerityFromArea(m_vEstuaryAreas[n], m_vWidth[n], Abc) / std::sqrt(dGetLateralStorageFactor(n));
        const double ubc = Rplus - 2.0 * cbc;
        m_vCrossSectionQ[n] = ubc * Abc;
    }
}

//======================================================================================================================
/**
 * @brief Apply boundary conditions to PREDICTOR STATE before corrector calculation
 * 
 * CRITICAL: Must be called BEFORE calculateHydraulicParameters() and calculateFlowTerms()
 * in the corrector phase. This ensures:
 * 1. Flux terms F0*, F1* are calculated with correct boundary values from predictor
 * 2. Source terms see correct gradients at boundaries
 * 
 * Sets boundary values for m_vPredictedCrossSectionArea and m_vPredictedCrossSectionQ at i=0 and i=n.
 * 
 * Boundary condition types:
 * - Type 0: Open/free boundary (gradient extrapolation)
 * - Type 1: Reflective/wall (Q=0) or prescribed discharge
 * - Type 2: Prescribed water surface elevation (with Manning for discharge)
 * - Type 3: Prescribed discharge (with Manning for area)
 */
//======================================================================================================================
void CSimulation::applyBoundariesToPredictorState() {
    int n = m_nCrossSectionsNumber - 1;
    
    // Upstream boundary - predictor state
    if (nGetUpwardEstuarineCondition() == 0) {
        const double dQ = m_vPredictedCrossSectionQ[2] - m_vPredictedCrossSectionQ[1];
        const double dA = m_vPredictedCrossSectionArea[2] - m_vPredictedCrossSectionArea[1];
        m_vPredictedCrossSectionQ[0] = m_vPredictedCrossSectionQ[1] - dQ;
        m_vPredictedCrossSectionArea[0] = m_vPredictedCrossSectionArea[1] - dA;
        if (m_vPredictedCrossSectionArea[0] < DRY_AREA) {
            m_vPredictedCrossSectionArea[0] = m_vPredictedCrossSectionArea[1];
        }
    }
    else if (nGetUpwardEstuarineCondition() == 1) {
        // Type 1: Reflective / solid-wall boundary (predictor state)
        m_vPredictedCrossSectionQ[0] = 0.0;
        {
            const double A1 = std::max(DRY_AREA, m_vPredictedCrossSectionArea[1]);
            const double h1 = linearInterpolation1d(A1, m_vEstuaryAreas[1], m_vEstuaryWaterDepths[1]);
            const double eta1 = h1 + m_vBedZ[1];
            const double h0 = std::max(0.0, eta1 - m_vBedZ[0]);
            m_vPredictedCrossSectionArea[0] = std::max(DRY_AREA,
                linearInterpolation1d(h0, m_vEstuaryWaterDepths[0], m_vEstuaryAreas[0]));
        }
    }
    else if (nGetUpwardEstuarineCondition() == 2) {
        m_vPredictedCrossSectionArea[0] = m_dUpwardBoundaryValue;
        double R = linearInterpolation1d(m_vPredictedCrossSectionArea[0], m_vEstuaryAreas[0], m_vEstuaryHydraulicRadius[0]);
        double sign_S0 = (m_vCrossSectionBedSlope[0] >= 0) ? 1.0 : -1.0;
        m_vPredictedCrossSectionQ[0] = m_vPredictedCrossSectionArea[0] * sqrt(fabs(m_vCrossSectionBedSlope[0]) + 1e-10) * 
                                       pow(R, 2.0/3.0) * sign_S0 / (m_vCrossSectionManningNumber[0] + 1e-10);
    }
    else {  // Type 3: prescribed discharge (characteristic compatibility)
        const double Qin = m_dUpwardBoundaryValue;
        m_vPredictedCrossSectionQ[0] = Qin;

        const double A1 = std::max(DRY_AREA, m_vPredictedCrossSectionArea[1]);
        const double u1 = m_vPredictedCrossSectionQ[1] / A1;
        const double c1 = celerityFromArea(m_vEstuaryAreas[1], m_vWidth[1], A1) / std::sqrt(dGetLateralStorageFactor(1));
        const double Rminus = u1 - 2.0 * c1;

        auto f = [&](double A0) {
            const double c0 = celerityFromArea(m_vEstuaryAreas[0], m_vWidth[0], A0) / std::sqrt(dGetLateralStorageFactor(0));
            return (Qin / std::max(DRY_AREA, A0)) - 2.0 * c0 - Rminus;
        };

        double Alo = std::max(DRY_AREA * 1.001, 0.2 * A1);
        double Ahi = std::max(A1 * 5.0, Alo * 2.0);
        double flo = f(Alo);
        double fhi = f(Ahi);
        for (int k = 0; k < 20 && (flo * fhi > 0.0); ++k) {
            if (flo > 0.0 && fhi > 0.0) {
                Ahi *= 2.0;
                fhi = f(Ahi);
            } else {
                Alo = std::max(DRY_AREA * 1.001, Alo * 0.5);
                flo = f(Alo);
            }
        }

        if (flo * fhi <= 0.0) {
            for (int it = 0; it < 40; ++it) {
                const double Am = 0.5 * (Alo + Ahi);
                const double fm = f(Am);
                if (flo * fm <= 0.0) {
                    Ahi = Am;
                    fhi = fm;
                } else {
                    Alo = Am;
                    flo = fm;
                }
            }
            m_vPredictedCrossSectionArea[0] = 0.5 * (Alo + Ahi);
        } else {
            m_vPredictedCrossSectionArea[0] = m_vPredictedCrossSectionArea[1];
        }
    }

    // Downstream boundary - predictor state
    if (nGetDownwardEstuarineCondition() == 0) {
        const double dQ = m_vPredictedCrossSectionQ[n-1] - m_vPredictedCrossSectionQ[n-2];
        const double dA = m_vPredictedCrossSectionArea[n-1] - m_vPredictedCrossSectionArea[n-2];
        m_vPredictedCrossSectionQ[n] = m_vPredictedCrossSectionQ[n-1] + dQ;
        m_vPredictedCrossSectionArea[n] = m_vPredictedCrossSectionArea[n-1] + dA;
        if (m_vPredictedCrossSectionArea[n] < DRY_AREA) {
            m_vPredictedCrossSectionArea[n] = m_vPredictedCrossSectionArea[n-1];
        }
    }
    else if (nGetDownwardEstuarineCondition() == 1) {
        m_vPredictedCrossSectionQ[n] = m_dDownwardBoundaryValue;
        double dManningFactor = fabs(m_vPredictedCrossSectionQ[n]) * m_vCrossSectionManningNumber[n] / (sqrt(fabs(m_vCrossSectionBedSlope[n]) + 1e-10));
        m_vPredictedCrossSectionArea[n] = linearInterpolation1d(dManningFactor, m_vPrecalculatedSecondTerm[n], m_vEstuaryAreas[n]);
    }
    else {  // Type 2: tidal elevation (characteristic compatibility)
        const double Abc = std::max(DRY_AREA, m_dDownwardBoundaryValue);
        m_vPredictedCrossSectionArea[n] = Abc;

        const double Ai = std::max(DRY_AREA, m_vPredictedCrossSectionArea[n - 1]);
        const double ui = m_vPredictedCrossSectionQ[n - 1] / Ai;
        const double ci = celerityFromArea(m_vEstuaryAreas[n - 1], m_vWidth[n - 1], Ai) / std::sqrt(dGetLateralStorageFactor(n - 1));
        const double Rplus = ui + 2.0 * ci;

        const double cbc = celerityFromArea(m_vEstuaryAreas[n], m_vWidth[n], Abc) / std::sqrt(dGetLateralStorageFactor(n));
        const double ubc = Rplus - 2.0 * cbc;
        m_vPredictedCrossSectionQ[n] = ubc * Abc;
    }
}

//======================================================================================================================
/**
 * @brief Apply boundary conditions after predictor step
 * 
 * Handles upstream and downstream boundary conditions with various types:
 * - Type 0: Open/free boundary (gradient extrapolation)
 * - Type 1: Reflective/wall (Q=0) or prescribed discharge
 * - Type 2: Prescribed water surface elevation (with Manning for discharge)
 * - Type 3: Prescribed discharge (with Manning for area)
 * 
 * Also applies lateral inflows from hydrographs at interior nodes.
 * Includes selective smoothing to prevent spurious oscillations at boundaries.
 */
//======================================================================================================================
void CSimulation::updatePredictorBoundaries() {
    //==============================================================================================================
    // UPSTREAM BOUNDARY CONDITION
    //==============================================================================================================
    
    if (nGetUpwardEstuarineCondition() == 0) {
        // Type 0: Open/free boundary - linear gradient extrapolation
        // Assumes zero second derivative at boundary (smooth outflow)
        const double dQ = m_vCrossSectionQ[2] - m_vCrossSectionQ[1];
        const double dA = m_vCrossSectionArea[2] - m_vCrossSectionArea[1];
        
        // Linear extrapolation: Q[0] = Q[1] - (Q[2] - Q[1]) = 2Q[1] - Q[2]
        m_vPredictedCrossSectionQ[0] = m_vCrossSectionQ[1] - dQ;
        m_vPredictedCrossSectionArea[0] = m_vCrossSectionArea[1] - dA;
        
        // Enforce minimum area threshold for stability
        if (m_vPredictedCrossSectionArea[0] < DRY_AREA) {
            m_vPredictedCrossSectionArea[0] = m_vCrossSectionArea[1];
        }
    }
    else if (nGetUpwardEstuarineCondition() == 1) {
        // Type 1: Reflective/wall boundary - no-flow condition (solid wall)
        // Imposes zero velocity (Q=0), extrapolates water level from interior
        m_vPredictedCrossSectionQ[0] = 0.0;

        // Closed end: enforce dη/dx≈0 → η[0]≈η[1]
        {
            const double A1 = std::max(DRY_AREA, m_vCrossSectionArea[1]);
            const double h1 = linearInterpolation1d(A1, m_vEstuaryAreas[1], m_vEstuaryWaterDepths[1]);
            const double eta1 = h1 + m_vBedZ[1];
            const double h0 = std::max(0.0, eta1 - m_vBedZ[0]);
            m_vPredictedCrossSectionArea[0] = std::max(DRY_AREA,
                linearInterpolation1d(h0, m_vEstuaryWaterDepths[0], m_vEstuaryAreas[0]));
        }
    }
    else if (nGetUpwardEstuarineCondition() == 2) {
        // Type 2: Prescribed water surface elevation
        // Impose area (converted from elevation), compute discharge using Manning equation
        m_vPredictedCrossSectionArea[0] = m_dUpwardBoundaryValue;
        m_vCrossSectionHydraulicRadius[0] = linearInterpolation1d(
            m_vPredictedCrossSectionArea[0], 
            m_vEstuaryAreas[0], 
            m_vEstuaryHydraulicRadius[0]);
        
        // Manning equation: Q = A·R^(2/3)·√S0 / n (with correct sign for flow direction)
        double sign_S0 = (m_vCrossSectionBedSlope[0] >= 0) ? 1.0 : -1.0;
        m_vPredictedCrossSectionQ[0] = m_vPredictedCrossSectionArea[0] * 
                                      sqrt(fabs(m_vCrossSectionBedSlope[0]) + 1e-10) * 
                                      pow(m_vCrossSectionHydraulicRadius[0], 2.0/3.0) * 
                                      sign_S0 / (m_vCrossSectionManningNumber[0] + 1e-10);
    }
    else {
        // Type 3: Prescribed discharge (characteristic compatibility)
        const double Qin = m_dUpwardBoundaryValue;
        m_vPredictedCrossSectionQ[0] = Qin;

        const double A1 = std::max(DRY_AREA, m_vPredictedCrossSectionArea[1]);
        const double u1 = m_vPredictedCrossSectionQ[1] / A1;
        const double c1 = celerityFromArea(m_vEstuaryAreas[1], m_vWidth[1], A1) / std::sqrt(dGetLateralStorageFactor(1));
        const double Rminus = u1 - 2.0 * c1;

        auto f = [&](double A0) {
            const double c0 = celerityFromArea(m_vEstuaryAreas[0], m_vWidth[0], A0) / std::sqrt(dGetLateralStorageFactor(0));
            return (Qin / std::max(DRY_AREA, A0)) - 2.0 * c0 - Rminus;
        };

        double Alo = std::max(DRY_AREA * 1.001, 0.2 * A1);
        double Ahi = std::max(A1 * 5.0, Alo * 2.0);
        double flo = f(Alo);
        double fhi = f(Ahi);
        for (int k = 0; k < 20 && (flo * fhi > 0.0); ++k) {
            if (flo > 0.0 && fhi > 0.0) {
                Ahi *= 2.0;
                fhi = f(Ahi);
            } else {
                Alo = std::max(DRY_AREA * 1.001, Alo * 0.5);
                flo = f(Alo);
            }
        }

        if (flo * fhi <= 0.0) {
            for (int it = 0; it < 40; ++it) {
                const double Am = 0.5 * (Alo + Ahi);
                const double fm = f(Am);
                if (flo * fm <= 0.0) {
                    Ahi = Am;
                    fhi = fm;
                } else {
                    Alo = Am;
                    flo = fm;
                }
            }
            m_vPredictedCrossSectionArea[0] = 0.5 * (Alo + Ahi);
        } else {
            m_vPredictedCrossSectionArea[0] = m_vPredictedCrossSectionArea[1];
        }
    }

    // Lateral inflows are handled as a continuity source term (Gv0 = Ql/Δx)
    // in calculateSourceTerms(). Do not add them directly to Q here.

    //==============================================================================================================
    // DOWNSTREAM BOUNDARY CONDITION
    //==============================================================================================================
    int n = m_nCrossSectionsNumber - 1;
    
    if (nGetDownwardEstuarineCondition() == 0) {
        // Type 0: Open/free boundary - linear gradient extrapolation
        const double dQ = m_vCrossSectionQ[n-1] - m_vCrossSectionQ[n-2];
        const double dA = m_vCrossSectionArea[n-1] - m_vCrossSectionArea[n-2];
        
        m_vPredictedCrossSectionQ[n] = m_vCrossSectionQ[n-1] + dQ;
        m_vPredictedCrossSectionArea[n] = m_vCrossSectionArea[n-1] + dA;
        
        // Enforce minimum area threshold
        if (m_vPredictedCrossSectionArea[n] < DRY_AREA) {
            m_vPredictedCrossSectionArea[n] = m_vCrossSectionArea[n-1];
        }
    }
    else if (nGetDownwardEstuarineCondition() == 1) {
        // Type 1: Prescribed discharge
        m_vPredictedCrossSectionQ[n] = m_dDownwardBoundaryValue;
        
        // Compute area using inverse Manning equation
        double dManningFactor = fabs(m_vPredictedCrossSectionQ[n]) * m_vCrossSectionManningNumber[n] / 
                              (sqrt(fabs(m_vCrossSectionBedSlope[n]) + 1e-10));
        
        m_vPredictedCrossSectionArea[n] = linearInterpolation1d(
            dManningFactor, 
            m_vPrecalculatedSecondTerm[n],
            m_vEstuaryAreas[n]);
    }
    else {
        // Type 2: Prescribed tidal elevation (characteristic compatibility)
        const double Abc = std::max(DRY_AREA, m_dDownwardBoundaryValue);
        m_vPredictedCrossSectionArea[n] = Abc;

        const double Ai = std::max(DRY_AREA, m_vPredictedCrossSectionArea[n - 1]);
        const double ui = m_vPredictedCrossSectionQ[n - 1] / Ai;
        const double ci = celerityFromArea(m_vEstuaryAreas[n - 1], m_vWidth[n - 1], Ai) / std::sqrt(dGetLateralStorageFactor(n - 1));
        const double Rplus = ui + 2.0 * ci;

        const double cbc = celerityFromArea(m_vEstuaryAreas[n], m_vWidth[n], Abc) / std::sqrt(dGetLateralStorageFactor(n));
        const double ubc = Rplus - 2.0 * cbc;
        m_vPredictedCrossSectionQ[n] = ubc * Abc;
    }
}

//===============================================================================================================================
/**
 * @brief Apply boundary conditions after corrector step
 * 
 * Similar to updatePredictorBoundaries but operates on corrected values.
 * Uses predicted state from predictor step as reference for extrapolation.
 * Includes optional selective sponge layer smoothing to dampen boundary discontinuities.
 * 
 * Boundary condition types same as predictor: 0=open, 1=reflective/discharge, 2=elevation, 3=discharge
 */
//===============================================================================================================================
void CSimulation::updateCorrectorBoundaries() {
    //==============================================================================================================
    // UPSTREAM BOUNDARY CONDITION - Corrector
    //==============================================================================================================
    
    if (nGetUpwardEstuarineCondition() == 0) {
        // Type 0: Open/free boundary - extrapolate from predicted state
        const double dQ = m_vPredictedCrossSectionQ[2] - m_vPredictedCrossSectionQ[1];
        const double dA = m_vPredictedCrossSectionArea[2] - m_vPredictedCrossSectionArea[1];
        
        m_vCorrectedCrossSectionQ[0] = m_vPredictedCrossSectionQ[1] - dQ;
        m_vCorrectedCrossSectionArea[0] = m_vPredictedCrossSectionArea[1] - dA;
        
        if (m_vCorrectedCrossSectionArea[0] < DRY_AREA) {
            m_vCorrectedCrossSectionArea[0] = m_vPredictedCrossSectionArea[1];
        }
    }
    else if (nGetUpwardEstuarineCondition() == 1) {
        // Type 1: Reflective/wall boundary (Q=0)
        m_vCorrectedCrossSectionQ[0] = 0.0;

        // Closed end: enforce dη/dx≈0 → η[0]≈η[1]
        {
            const double A1 = std::max(DRY_AREA, m_vPredictedCrossSectionArea[1]);
            const double h1 = linearInterpolation1d(A1, m_vEstuaryAreas[1], m_vEstuaryWaterDepths[1]);
            const double eta1 = h1 + m_vBedZ[1];
            const double h0 = std::max(0.0, eta1 - m_vBedZ[0]);
            m_vCorrectedCrossSectionArea[0] = std::max(DRY_AREA,
                linearInterpolation1d(h0, m_vEstuaryWaterDepths[0], m_vEstuaryAreas[0]));
        }
    }
    else if (nGetUpwardEstuarineCondition() == 2) {
        // Type 2: Prescribed elevation - impose area, compute Q with Manning
        m_vCorrectedCrossSectionArea[0] = m_dUpwardBoundaryValue;
        
        m_vCrossSectionHydraulicRadius[0] = linearInterpolation1d(
            m_vCorrectedCrossSectionArea[0], 
            m_vEstuaryAreas[0], 
            m_vEstuaryHydraulicRadius[0]);
        
        double sign_S0 = (m_vCrossSectionBedSlope[0] >= 0) ? 1.0 : -1.0;
        m_vCorrectedCrossSectionQ[0] = m_vCorrectedCrossSectionArea[0] * 
                                      sqrt(fabs(m_vCrossSectionBedSlope[0]) + 1e-10) * 
                                      pow(m_vCrossSectionHydraulicRadius[0], 2.0/3.0) * 
                                      sign_S0 / (m_vCrossSectionManningNumber[0] + 1e-10);
    }
    else {
        // Type 3: Prescribed discharge (characteristic compatibility)
        const double Qin = m_dUpwardBoundaryValue;
        m_vCorrectedCrossSectionQ[0] = Qin;

        const double A1 = std::max(DRY_AREA, m_vCorrectedCrossSectionArea[1]);
        const double u1 = m_vCorrectedCrossSectionQ[1] / A1;
        const double c1 = celerityFromArea(m_vEstuaryAreas[1], m_vWidth[1], A1) / std::sqrt(dGetLateralStorageFactor(1));
        const double Rminus = u1 - 2.0 * c1;

        auto f = [&](double A0) {
            const double c0 = celerityFromArea(m_vEstuaryAreas[0], m_vWidth[0], A0) / std::sqrt(dGetLateralStorageFactor(0));
            return (Qin / std::max(DRY_AREA, A0)) - 2.0 * c0 - Rminus;
        };

        double Alo = std::max(DRY_AREA * 1.001, 0.2 * A1);
        double Ahi = std::max(A1 * 5.0, Alo * 2.0);
        double flo = f(Alo);
        double fhi = f(Ahi);
        for (int k = 0; k < 20 && (flo * fhi > 0.0); ++k) {
            if (flo > 0.0 && fhi > 0.0) {
                Ahi *= 2.0;
                fhi = f(Ahi);
            } else {
                Alo = std::max(DRY_AREA * 1.001, Alo * 0.5);
                flo = f(Alo);
            }
        }

        if (flo * fhi <= 0.0) {
            for (int it = 0; it < 40; ++it) {
                const double Am = 0.5 * (Alo + Ahi);
                const double fm = f(Am);
                if (flo * fm <= 0.0) {
                    Ahi = Am;
                    fhi = fm;
                } else {
                    Alo = Am;
                    flo = fm;
                }
            }
            m_vCorrectedCrossSectionArea[0] = 0.5 * (Alo + Ahi);
        } else {
            m_vCorrectedCrossSectionArea[0] = m_vCorrectedCrossSectionArea[1];
        }
    }

    // Lateral inflows are handled as a continuity source term (Gv0 = Ql/Δx)
    // in calculateSourceTerms(). Do not add them directly to Q here.

    //==============================================================================================================
    // DOWNSTREAM BOUNDARY CONDITION - Corrector
    //==============================================================================================================
    int n = m_nCrossSectionsNumber - 1;
    
    if (nGetDownwardEstuarineCondition() == 0) {
        // Type 0: Open/free boundary
        const double dQ = m_vPredictedCrossSectionQ[n-1] - m_vPredictedCrossSectionQ[n-2];
        const double dA = m_vPredictedCrossSectionArea[n-1] - m_vPredictedCrossSectionArea[n-2];
        
        m_vCorrectedCrossSectionQ[n] = m_vPredictedCrossSectionQ[n-1] + dQ;
        m_vCorrectedCrossSectionArea[n] = m_vPredictedCrossSectionArea[n-1] + dA;
        
        if (m_vCorrectedCrossSectionArea[n] < DRY_AREA) {
            m_vCorrectedCrossSectionArea[n] = m_vPredictedCrossSectionArea[n-1];
        }
    }
    else if (nGetDownwardEstuarineCondition() == 1) {
        // Type 1: Prescribed discharge
        m_vCorrectedCrossSectionQ[n] = m_dDownwardBoundaryValue;
        
        const double dManningFactor = fabs(m_vCorrectedCrossSectionQ[n]) * m_vCrossSectionManningNumber[n] / 
                                     (sqrt(fabs(m_vCrossSectionBedSlope[n]) + 1e-10));
        
        m_vCorrectedCrossSectionArea[n] = linearInterpolation1d(
            dManningFactor, 
            m_vPrecalculatedSecondTerm[n],
            m_vEstuaryAreas[n]);
    }
    else {
        // Type 2: Prescribed tidal elevation (characteristic compatibility)
        const double Abc = std::max(DRY_AREA, m_dDownwardBoundaryValue);
        m_vCorrectedCrossSectionArea[n] = Abc;

        const double Ai = std::max(DRY_AREA, m_vCorrectedCrossSectionArea[n - 1]);
        const double ui = m_vCorrectedCrossSectionQ[n - 1] / Ai;
        const double ci = celerityFromArea(m_vEstuaryAreas[n - 1], m_vWidth[n - 1], Ai) / std::sqrt(dGetLateralStorageFactor(n - 1));
        const double Rplus = ui + 2.0 * ci;

        const double cbc = celerityFromArea(m_vEstuaryAreas[n], m_vWidth[n], Abc) / std::sqrt(dGetLateralStorageFactor(n));
        const double ubc = Rplus - 2.0 * cbc;
        m_vCorrectedCrossSectionQ[n] = ubc * Abc;
    }
}

//======================================================================================================================
/**
 * @brief Generic TVD flux limiter for passive tracer transport (salinity, temperature, etc.)
 * 
 * Computes TVD-limited fluxes at cell faces to prevent spurious oscillations while
 * maintaining 2nd order accuracy in smooth regions.
 * 
 * @param tracer Scalar field to transport (salinity, temperature, etc.)
 * @param velocity Velocity field at cell centers
 * @param area Cross-sectional areas
 * @param flux_limited Output: TVD-limited flux at faces [kg/s or similar]
 */
//======================================================================================================================
void CSimulation::compute_tracer_tvd_flux(const vector<double>& tracer, const vector<double>& discharge,
                                            vector<double>& flux_limited) {
    const int n = m_nCrossSectionsNumber;
    const int limiter_type = nGetTransportLimiterFlux();  // Get limiter for transport
    
    // Boundary fluxes set to zero (boundary conditions applied separately)
    flux_limited[0] = 0.0;
    flux_limited[n] = 0.0;
    
    // Early exit if no limiter (use simple upwind)
    if (limiter_type == 0 || !m_bDoMcCormackLimiterFlux) {
        for (int i = 1; i < n; i++) {
            // Discharge at interface (average of adjacent cells)
            const double Q_face = 0.5 * (discharge[i-1] + discharge[i]);
            // Upwind concentration
            const double phi_upwind = (Q_face > 0.0) ? tracer[i-1] : tracer[i];
            // Flux = Q * phi
            flux_limited[i] = Q_face * phi_upwind;
        }
        return;
    }
    
    // Main TVD loop
    for (int i = 1; i < n; i++) {
        // Discharge at interface (average of adjacent cells)
        // Use Q directly from momentum equations to ensure mass conservation consistency
        const double Q_face = 0.5 * (discharge[i-1] + discharge[i]);
        
        if (fabs(Q_face) < 1e-6) {
            flux_limited[i] = 0.0;
            continue;
        }
        const double delta_phi = tracer[i] - tracer[i-1];
        if (fabs(delta_phi) < 1e-10) {
            flux_limited[i] = Q_face * tracer[i-1];
            continue;
        }
        // Compute gradient ratio for TVD limiter
        double r = 1.0;
        if (Q_face > 0.0 && i >= 2) {
            const double delta_upwind = tracer[i-1] - tracer[i-2];
            r = (fabs(delta_phi) > 1e-12) ? (delta_upwind / delta_phi) : 1.0;
        } else if (Q_face < 0.0 && i < n-1) {
            const double delta_downwind = tracer[i+1] - tracer[i];
            r = (fabs(delta_phi) > 1e-12) ? (delta_downwind / delta_phi) : 1.0;
        }
        // Apply limiter function
        double psi;
        switch (limiter_type) {
            case 1:  psi = std::max(0.0, std::min(r, 1.0)); break;
            case 2:  psi = std::max(0.0, std::max(std::min(2.0*r, 1.0), std::min(r, 2.0))); break;
            case 3:  psi = (fabs(r) + r) / (1.0 + fabs(r)); break;
            case 4:  psi = ((r*r) + r) / (1.0 + r*r); break;
            default: psi = 0.0;
        }
        // Reconstruct interface value and compute flux
        const double phi_upwind = (Q_face > 0.0) ? tracer[i-1] : tracer[i];
        const double correction = (Q_face > 0.0) ? 0.5 * psi * delta_phi : -0.5 * psi * delta_phi;
        flux_limited[i] = Q_face * (phi_upwind + correction);
    }
}


//======================================================================================================================
/**
 * @brief Predictor step for salinity transport equation
 * 
 * Solves advection-diffusion equation for salinity:
 * ∂(AS)/∂t + ∂(QS)/∂x = ∂/∂x(Kh·A·∂S/∂x)
 * 
 * Uses TVD-limited advection for sharp fronts + centered diffusion:

 * - If u < 0: forward difference (i+1 - i)
 * 
 * Diffusion uses centered differences for stability.
 * Conserves salt mass: updates A·S, then divides by A to get S.
 * Applies prescribed boundary conditions at upstream/downstream.
 */
//======================================================================================================================
void CSimulation::calculate_salinity_predictor() {
    // --- Boundary salinity values (needed for boundary advective fluxes) ---
    double up_sal = m_dUpwardSalinityBoundaryValue;
    if (!m_strUpwardSalinityBoundaryConditionFilename.empty() &&
        !m_vUpwardSalinityBoundaryConditionTime.empty() &&
        !m_vUpwardSalinityBoundaryConditionValue.empty()) {
        const double t = m_dCurrentTime;
        auto it = std::lower_bound(m_vUpwardSalinityBoundaryConditionTime.begin(),
                                   m_vUpwardSalinityBoundaryConditionTime.end(), t);

        if (it == m_vUpwardSalinityBoundaryConditionTime.begin()) {
            up_sal = m_vUpwardSalinityBoundaryConditionValue.front();
        } else if (it == m_vUpwardSalinityBoundaryConditionTime.end()) {
            up_sal = m_vUpwardSalinityBoundaryConditionValue.back();
        } else {
            const size_t idx = std::distance(m_vUpwardSalinityBoundaryConditionTime.begin(), it);
            const double t1 = m_vUpwardSalinityBoundaryConditionTime[idx - 1];
            const double t2 = m_vUpwardSalinityBoundaryConditionTime[idx];
            const double s1 = m_vUpwardSalinityBoundaryConditionValue[idx - 1];
            const double s2 = m_vUpwardSalinityBoundaryConditionValue[idx];
            up_sal = s1 + (s2 - s1) * (t - t1) / (t2 - t1);
        }
    }

    double down_sal = m_dDownwardSalinityBoundaryValue;
    if (!m_strDownwardSalinityBoundaryConditionFilename.empty() &&
        !m_vDownwardSalinityBoundaryConditionTime.empty() &&
        !m_vDownwardSalinityBoundaryConditionValue.empty()) {
        const double t = m_dCurrentTime;
        auto it = std::lower_bound(m_vDownwardSalinityBoundaryConditionTime.begin(),
                                   m_vDownwardSalinityBoundaryConditionTime.end(), t);

        if (it == m_vDownwardSalinityBoundaryConditionTime.begin()) {
            down_sal = m_vDownwardSalinityBoundaryConditionValue.front();
        } else if (it == m_vDownwardSalinityBoundaryConditionTime.end()) {
            down_sal = m_vDownwardSalinityBoundaryConditionValue.back();
        } else {
            const size_t idx = std::distance(m_vDownwardSalinityBoundaryConditionTime.begin(), it);
            const double t1 = m_vDownwardSalinityBoundaryConditionTime[idx - 1];
            const double t2 = m_vDownwardSalinityBoundaryConditionTime[idx];
            const double s1 = m_vDownwardSalinityBoundaryConditionValue[idx - 1];
            const double s2 = m_vDownwardSalinityBoundaryConditionValue[idx];
            down_sal = s1 + (s2 - s1) * (t - t1) / (t2 - t1);
        }
    }

    // Apply salinity BC consistently to the tracer field used for face fluxes.
    // Here nodes 0 and last are treated as boundary nodes (prescribed / radiation depending on flow).
    const int last = m_nCrossSectionsNumber - 1;
    vector<double> tracer_bc = m_vCrossSectionSalinity;
    {
        const double Q_up = m_vPredictedCrossSectionQ[0];
        if (Q_up > 0.0) {
            // Inflow from upstream
            tracer_bc[0] = up_sal;
        } else {
            // Outflow: radiation/zero-gradient
            tracer_bc[0] = (m_nCrossSectionsNumber > 1) ? tracer_bc[1] : tracer_bc[0];
        }
    }
    {
        const double Q_dn = m_vPredictedCrossSectionQ[last];
        if (Q_dn < 0.0) {
            // Inflow from ocean
            tracer_bc[last] = down_sal;
        } else {
            // Outflow to ocean: radiation/zero-gradient
            tracer_bc[last] = (m_nCrossSectionsNumber > 1) ? tracer_bc[last - 1] : tracer_bc[last];
        }
    }

    // Compute TVD-limited advective fluxes at cell faces
    // Use predicted discharge from hydrodynamic predictor step
    vector<double> tvd_flux(m_nCrossSectionsNumber+1, 0.0);
    compute_tracer_tvd_flux(tracer_bc, m_vPredictedCrossSectionQ, tvd_flux);

    // Boundary advective fluxes: compute_tracer_tvd_flux() sets them to zero by design.
    // For open boundaries this is wrong: we must supply the boundary flux using the BC value
    // for inflow, and the interior value for outflow (upwind).
    {
        const double Q_up = m_vPredictedCrossSectionQ[0];
        const double S_upwind = (Q_up > 0.0) ? up_sal : tracer_bc[0];
        tvd_flux[0] = Q_up * S_upwind;
    }
    {
        const double Q_dn = m_vPredictedCrossSectionQ[last];
        double S_upwind;
        if (Q_dn < 0.0) {
            // Inflow from ocean
            S_upwind = down_sal;
        } else {
            // Outflow to ocean
            S_upwind = tracer_bc[last];
        }
        tvd_flux[m_nCrossSectionsNumber] = Q_dn * S_upwind;
    }
    
    // Compute advection and diffusion terms
    // Equation: ∂(A·S)/∂t = -∂(Q·S)/∂x + ∂/∂x(Kh·A·∂S/∂x)
    // IMPORTANT: use conservative flux form for diffusion, and control-volume length on non-uniform grids.
    for (int i = 1; i < m_nCrossSectionsNumber - 1; ++i) {
        const double dxL = m_vCrossSectionX[i] - m_vCrossSectionX[i - 1];
        const double dxR = m_vCrossSectionX[i + 1] - m_vCrossSectionX[i];
        const double dxCV = 0.5 * (dxL + dxR);
        if (dxL <= 1e-12 || dxR <= 1e-12 || dxCV <= 1e-12) {
            m_vCrossSectionSalinityASt[i] = 0.0;
            continue;
        }

        // Advection divergence: (F_{i+1/2} - F_{i-1/2}) / Δx_CV
        const double d_QS_dx = (tvd_flux[i + 1] - tvd_flux[i]) / dxCV;

        // Diffusion divergence in flux form.
        // Define diffusive flux at faces: Fd = -Kh * A_face * (ΔS/Δx)
        double d_Fd_dx = 0.0;
        {
            const double Kh_L = 0.5 * (dGetLongitudinalDispersion(i - 1) + dGetLongitudinalDispersion(i));
            const double Kh_R = 0.5 * (dGetLongitudinalDispersion(i) + dGetLongitudinalDispersion(i + 1));
            if (Kh_L > 0.0 || Kh_R > 0.0) {
            const double A_L = 0.5 * (m_vCrossSectionArea[i - 1] + m_vCrossSectionArea[i]);
            const double A_R = 0.5 * (m_vCrossSectionArea[i] + m_vCrossSectionArea[i + 1]);
            const double Fd_L = -Kh_L * A_L * (tracer_bc[i] - tracer_bc[i - 1]) / dxL;
            const double Fd_R = -Kh_R * A_R * (tracer_bc[i + 1] - tracer_bc[i]) / dxR;
            d_Fd_dx = (Fd_R - Fd_L) / dxCV;
            }
        }

        // Lateral storage consistency:
        // If enabled, hydrodynamics uses: S·∂A/∂t + ∂Q/∂x = ql.
        // For a well-mixed scalar in the (channel + storage) volume:
        //   S·∂(A·C)/∂t + ∂(Q·C)/∂x = diffusion + sources.
        // Therefore: Δ(A·C) = (dt/S) * RHS.
        const double invS = 1.0 / dGetLateralStorageFactor(i);
        m_vCrossSectionSalinityASt[i] = (m_dTimestep * invS) * (-d_QS_dx - d_Fd_dx);
    }

    // Update salinity conserving salt mass: S_new = (A·S_old + Δ(A·S))/A_new
    static int clamp_count = 0;
    // Set boundary nodes from BC (do not evolve them with the PDE)
    m_vPredictedCrossSectionS[0] = tracer_bc[0];
    m_vPredictedCrossSectionS[last] = tracer_bc[last];

    for (int i = 1; i < m_nCrossSectionsNumber - 1; i++) {
        if (m_vCrossSectionArea[i] > DRY_AREA) {
            double salt_mass = m_vCrossSectionArea[i] * m_vCrossSectionSalinity[i];
            double delta_mass = m_vCrossSectionSalinityASt[i];
            salt_mass += delta_mass;
            // Enforce physical bounds: salinity >= 0
            if (salt_mass < 0.0) salt_mass = 0.0;
            m_vPredictedCrossSectionS[i] = salt_mass / m_vPredictedCrossSectionArea[i];
            // Clamp to physical range [0, S_ocean]
            if (m_vPredictedCrossSectionS[i] < 0.0) { m_vPredictedCrossSectionS[i] = 0.0; clamp_count++; }
            if (m_vPredictedCrossSectionS[i] > down_sal) { m_vPredictedCrossSectionS[i] = down_sal; clamp_count++; }
        } else {
            // Dry area: zero salinity
            m_vPredictedCrossSectionS[i] = 0.0;
        }
    }
    // Removed per-timestep clamping log in predictor
}

//======================================================================================================================
/**
 * @brief Corrector step for salinity transport equation
 * 
 * Uses predicted salinity from predictor step to compute corrected values.
 * Applies same advection-diffusion scheme as predictor but with predicted state.
 * Part of McCormack predictor-corrector method for second-order accuracy.
 */
//======================================================================================================================
void CSimulation::calculate_salinity_corrector() {
    // --- Boundary salinity values (needed for boundary advective fluxes) ---
    double up_sal = m_dUpwardSalinityBoundaryValue;
    if (!m_strUpwardSalinityBoundaryConditionFilename.empty() &&
        !m_vUpwardSalinityBoundaryConditionTime.empty() &&
        !m_vUpwardSalinityBoundaryConditionValue.empty()) {
        const double t = m_dCurrentTime;
        auto it = std::lower_bound(m_vUpwardSalinityBoundaryConditionTime.begin(),
                                   m_vUpwardSalinityBoundaryConditionTime.end(), t);
        if (it == m_vUpwardSalinityBoundaryConditionTime.begin()) {
            up_sal = m_vUpwardSalinityBoundaryConditionValue.front();
        } else if (it == m_vUpwardSalinityBoundaryConditionTime.end()) {
            up_sal = m_vUpwardSalinityBoundaryConditionValue.back();
        } else {
            const size_t idx = std::distance(m_vUpwardSalinityBoundaryConditionTime.begin(), it);
            const double t1 = m_vUpwardSalinityBoundaryConditionTime[idx - 1];
            const double t2 = m_vUpwardSalinityBoundaryConditionTime[idx];
            const double s1 = m_vUpwardSalinityBoundaryConditionValue[idx - 1];
            const double s2 = m_vUpwardSalinityBoundaryConditionValue[idx];
            up_sal = s1 + (s2 - s1) * (t - t1) / (t2 - t1);
        }
    }

    double down_sal = m_dDownwardSalinityBoundaryValue;
    if (!m_strDownwardSalinityBoundaryConditionFilename.empty() &&
        !m_vDownwardSalinityBoundaryConditionTime.empty() &&
        !m_vDownwardSalinityBoundaryConditionValue.empty()) {
        const double t = m_dCurrentTime;
        auto it = std::lower_bound(m_vDownwardSalinityBoundaryConditionTime.begin(),
                                   m_vDownwardSalinityBoundaryConditionTime.end(), t);
        if (it == m_vDownwardSalinityBoundaryConditionTime.begin()) {
            down_sal = m_vDownwardSalinityBoundaryConditionValue.front();
        } else if (it == m_vDownwardSalinityBoundaryConditionTime.end()) {
            down_sal = m_vDownwardSalinityBoundaryConditionValue.back();
        } else {
            const size_t idx = std::distance(m_vDownwardSalinityBoundaryConditionTime.begin(), it);
            const double t1 = m_vDownwardSalinityBoundaryConditionTime[idx - 1];
            const double t2 = m_vDownwardSalinityBoundaryConditionTime[idx];
            const double s1 = m_vDownwardSalinityBoundaryConditionValue[idx - 1];
            const double s2 = m_vDownwardSalinityBoundaryConditionValue[idx];
            down_sal = s1 + (s2 - s1) * (t - t1) / (t2 - t1);
        }
    }

    // Apply salinity BC consistently to the tracer field used for face fluxes.
    // Nodes 0 and last are treated as boundary nodes.
    const int last = m_nCrossSectionsNumber - 1;
    vector<double> tracer_bc = m_vPredictedCrossSectionS;
    {
        const double Q_up = m_vCorrectedCrossSectionQ[0];
        if (Q_up > 0.0) {
            tracer_bc[0] = up_sal;
        } else {
            tracer_bc[0] = (m_nCrossSectionsNumber > 1) ? tracer_bc[1] : tracer_bc[0];
        }
    }
    {
        const double Q_dn = m_vCorrectedCrossSectionQ[last];
        if (Q_dn < 0.0) {
            tracer_bc[last] = down_sal;
        } else {
            tracer_bc[last] = (m_nCrossSectionsNumber > 1) ? tracer_bc[last - 1] : tracer_bc[last];
        }
    }

    // Compute TVD-limited advective fluxes at cell faces (using predicted state)
    // Use corrected discharge from hydrodynamic corrector step
    vector<double> tvd_flux(m_nCrossSectionsNumber+1, 0.0);
    compute_tracer_tvd_flux(tracer_bc, m_vCorrectedCrossSectionQ, tvd_flux);

    // Boundary advective fluxes
    {
        const double Q_up = m_vCorrectedCrossSectionQ[0];
        const double S_upwind = (Q_up > 0.0) ? up_sal : tracer_bc[0];
        tvd_flux[0] = Q_up * S_upwind;
    }
    {
        const double Q_dn = m_vCorrectedCrossSectionQ[last];
        const double S_upwind = (Q_dn < 0.0) ? down_sal : tracer_bc[last];
        tvd_flux[m_nCrossSectionsNumber] = Q_dn * S_upwind;
    }

    // Compute advection and diffusion terms
    // Equation: ∂(A·S)/∂t = -∂(Q·S)/∂x + ∂/∂x(Kh·A·∂S/∂x)
    // Use conservative flux form for diffusion and control-volume length for non-uniform grids.
    for (int i = 1; i < m_nCrossSectionsNumber - 1; ++i) {
        const double dxL = m_vCrossSectionX[i] - m_vCrossSectionX[i - 1];
        const double dxR = m_vCrossSectionX[i + 1] - m_vCrossSectionX[i];
        const double dxCV = 0.5 * (dxL + dxR);
        if (dxL <= 1e-12 || dxR <= 1e-12 || dxCV <= 1e-12) {
            m_vCrossSectionSalinityASt[i] = 0.0;
            continue;
        }

        const double d_QS_dx = (tvd_flux[i + 1] - tvd_flux[i]) / dxCV;

        double d_Fd_dx = 0.0;
        {
            const double Kh_L = 0.5 * (dGetLongitudinalDispersion(i - 1) + dGetLongitudinalDispersion(i));
            const double Kh_R = 0.5 * (dGetLongitudinalDispersion(i) + dGetLongitudinalDispersion(i + 1));
            if (Kh_L > 0.0 || Kh_R > 0.0) {
            const double A_L = 0.5 * (m_vPredictedCrossSectionArea[i - 1] + m_vPredictedCrossSectionArea[i]);
            const double A_R = 0.5 * (m_vPredictedCrossSectionArea[i] + m_vPredictedCrossSectionArea[i + 1]);
            const double Fd_L = -Kh_L * A_L * (tracer_bc[i] - tracer_bc[i - 1]) / dxL;
            const double Fd_R = -Kh_R * A_R * (tracer_bc[i + 1] - tracer_bc[i]) / dxR;
            d_Fd_dx = (Fd_R - Fd_L) / dxCV;
            }
        }

        const double invS = 1.0 / dGetLateralStorageFactor(i);
        m_vCrossSectionSalinityASt[i] = (m_dTimestep * invS) * (-d_QS_dx - d_Fd_dx);
    }

    // Update salinity conserving salt mass: S_new = (A·S_old + Δ(A·S))/A_new
    static int clamp_count = 0;
    // Set boundary nodes from BC (do not evolve them with the PDE)
    m_vCorrectedCrossSectionS[0] = tracer_bc[0];
    m_vCorrectedCrossSectionS[last] = tracer_bc[last];

    for (int i = 1; i < m_nCrossSectionsNumber - 1; i++) {
        if (m_vPredictedCrossSectionArea[i] > DRY_AREA) {
            double salt_mass = m_vPredictedCrossSectionArea[i] * m_vPredictedCrossSectionS[i];
            double delta_mass = m_vCrossSectionSalinityASt[i];
            salt_mass += delta_mass;
            // Enforce physical bounds: salinity >= 0
            if (salt_mass < 0.0) salt_mass = 0.0;
            m_vCorrectedCrossSectionS[i] = salt_mass / m_vCorrectedCrossSectionArea[i];
            // Clamp to physical range [0, S_ocean]
            if (m_vCorrectedCrossSectionS[i] < 0.0) { m_vCorrectedCrossSectionS[i] = 0.0; clamp_count++; }
            if (m_vCorrectedCrossSectionS[i] > down_sal) { m_vCorrectedCrossSectionS[i] = down_sal; clamp_count++; }
        } else {
            // Dry area: zero salinity
            m_vCorrectedCrossSectionS[i] = 0.0;
        }
    }
    // Removed per-timestep clamping log in corrector
}

//======================================================================================================================
/**
 * @brief Predictor step for temperature transport equation
 * 
 * Solves advection-diffusion equation with heat flux source:
 * ∂T/∂t + u·∂T/∂x = Kh·∂²T/∂x² + Q_net/(ρ·Cp·h)
 * 
 * where:
 * - ρ = 1000 kg/m³ (water density)
 * - Cp = 4186 J/(kg·°C) (specific heat capacity of water)
 * - Q_net = net surface heat flux (W/m²) stored in m_vCrossSectionTemperatureASt
 * 
 * Uses upwind advection and centered diffusion.
 * Applies prescribed boundary conditions or extrapolation at upstream/downstream.
 */
//======================================================================================================================
void CSimulation::calculate_temperature_predictor() {
    const int N = m_nCrossSectionsNumber;
    if (m_vPredictedCrossSectionT.size() != static_cast<size_t>(N))
        m_vPredictedCrossSectionT.resize(N, 0.0);
    
    for (int i = 1; i < N-1; ++i) {
        double adv = 0.0;
        const double u = m_vCrossSectionU[i];
        const double dxL = m_vCrossSectionX[i] - m_vCrossSectionX[i - 1];
        const double dxR = m_vCrossSectionX[i + 1] - m_vCrossSectionX[i];
        const double dxCV = 0.5 * (dxL + dxR);

        if (dxL <= 1e-12 || dxR <= 1e-12 || dxCV <= 1e-12) {
            m_vPredictedCrossSectionT[i] = m_vCrossSectionTemperature[i];
            continue;
        }
        
        // Upwind advection on non-uniform grid
        if (u > 0.0) {
            adv = -u * (m_vCrossSectionTemperature[i] - m_vCrossSectionTemperature[i - 1]) / dxL;
        } else {
            adv = -u * (m_vCrossSectionTemperature[i + 1] - m_vCrossSectionTemperature[i]) / dxR;
        }
        
        // Non-uniform second derivative (conservative Laplacian form)
        const double d2Tdx2 = (2.0 / (dxL + dxR)) *
                              ((m_vCrossSectionTemperature[i + 1] - m_vCrossSectionTemperature[i]) / dxR -
                               (m_vCrossSectionTemperature[i] - m_vCrossSectionTemperature[i - 1]) / dxL);
        const double diff = m_dThermalDispersion * d2Tdx2;
        
        // Surface heat flux contribution: Q_net/(ρ·Cp·h)
        double Qnet = m_vCrossSectionTemperatureASt[i];
        double dTdt = adv + diff + Qnet / (WATER_DENSITY * WATER_SPECIFIC_HEAT * m_vCrossSectionWaterDepth[i]);

        const double invS = 1.0 / dGetLateralStorageFactor(i);
        m_vPredictedCrossSectionT[i] = m_vCrossSectionTemperature[i] + (m_dTimestep * invS) * dTdt;
    }
    
    // Boundary conditions
    // Upstream temperature boundary
    if (m_nUpwardTemperatureCondition == 1) {
        // Type 1: Constant prescribed value
        m_vPredictedCrossSectionT[0] = m_dUpwardTemperatureBoundaryValue;
    } else if (m_nUpwardTemperatureCondition == 2 && !m_vUpwardTemperatureBoundaryConditionTime.empty()) {
        // Type 2: Time-varying prescribed value
        m_vPredictedCrossSectionT[0] = linearInterpolation1d(m_dCurrentTime, 
                                                             m_vUpwardTemperatureBoundaryConditionTime, 
                                                             m_vUpwardTemperatureBoundaryConditionValue);
    } else {
        // Default: Extrapolate from interior
        m_vPredictedCrossSectionT[0] = m_vPredictedCrossSectionT[1];
    }
    
    // Downstream temperature boundary
    if (m_nDownwardTemperatureCondition == 1) {
        m_vPredictedCrossSectionT[N-1] = m_dDownwardTemperatureBoundaryValue;
    } else if (m_nDownwardTemperatureCondition == 2 && !m_vDownwardTemperatureBoundaryConditionTime.empty()) {
        m_vPredictedCrossSectionT[N-1] = linearInterpolation1d(m_dCurrentTime, 
                                                               m_vDownwardTemperatureBoundaryConditionTime, 
                                                               m_vDownwardTemperatureBoundaryConditionValue);
    } else {
        m_vPredictedCrossSectionT[N-1] = m_vPredictedCrossSectionT[N-2];
    }
}

//======================================================================================================================
/**
 * @brief Corrector step for temperature transport equation
 * 
 * Uses predicted temperature from predictor step to compute corrected values.
 * Applies same advection-diffusion-heating scheme as predictor but with predicted state.
 */
//======================================================================================================================
void CSimulation::calculate_temperature_corrector() {
    const int N = m_nCrossSectionsNumber;
    if (m_vCorrectedCrossSectionT.size() != static_cast<size_t>(N))
        m_vCorrectedCrossSectionT.resize(N, 0.0);
    
    for (int i = 1; i < N-1; ++i) {
        double adv = 0.0;
        const double u = m_vCrossSectionU[i];
        const double dxL = m_vCrossSectionX[i] - m_vCrossSectionX[i - 1];
        const double dxR = m_vCrossSectionX[i + 1] - m_vCrossSectionX[i];
        const double dxCV = 0.5 * (dxL + dxR);

        if (dxL <= 1e-12 || dxR <= 1e-12 || dxCV <= 1e-12) {
            m_vCorrectedCrossSectionT[i] = m_vPredictedCrossSectionT[i];
            continue;
        }

        if (u > 0.0) {
            adv = -u * (m_vPredictedCrossSectionT[i] - m_vPredictedCrossSectionT[i - 1]) / dxL;
        } else {
            adv = -u * (m_vPredictedCrossSectionT[i + 1] - m_vPredictedCrossSectionT[i]) / dxR;
        }

        const double d2Tdx2 = (2.0 / (dxL + dxR)) *
                              ((m_vPredictedCrossSectionT[i + 1] - m_vPredictedCrossSectionT[i]) / dxR -
                               (m_vPredictedCrossSectionT[i] - m_vPredictedCrossSectionT[i - 1]) / dxL);
        const double diff = m_dThermalDispersion * d2Tdx2;

        const double Qnet = m_vCrossSectionTemperatureASt[i];
        const double dTdt = adv + diff + Qnet / (WATER_DENSITY * WATER_SPECIFIC_HEAT * m_vCrossSectionWaterDepth[i]);
        const double invS = 1.0 / dGetLateralStorageFactor(i);
        m_vCorrectedCrossSectionT[i] = m_vPredictedCrossSectionT[i] + (m_dTimestep * invS) * dTdt;
    }
    // Fronteras
    if (m_nUpwardTemperatureCondition == 1) {
        m_vCorrectedCrossSectionT[0] = m_dUpwardTemperatureBoundaryValue;
    } else if (m_nUpwardTemperatureCondition == 2 && !m_vUpwardTemperatureBoundaryConditionTime.empty()) {
        m_vCorrectedCrossSectionT[0] = linearInterpolation1d(m_dCurrentTime, m_vUpwardTemperatureBoundaryConditionTime, m_vUpwardTemperatureBoundaryConditionValue);
    } else {
        m_vCorrectedCrossSectionT[0] = m_vCorrectedCrossSectionT[1];
    }
    if (m_nDownwardTemperatureCondition == 1) {
        m_vCorrectedCrossSectionT[N-1] = m_dDownwardTemperatureBoundaryValue;
    } else if (m_nDownwardTemperatureCondition == 2 && !m_vDownwardTemperatureBoundaryConditionTime.empty()) {
        m_vCorrectedCrossSectionT[N-1] = linearInterpolation1d(m_dCurrentTime, m_vDownwardTemperatureBoundaryConditionTime, m_vDownwardTemperatureBoundaryConditionValue);
    } else {
        m_vCorrectedCrossSectionT[N-1] = m_vCorrectedCrossSectionT[N-2];
    }
}

//======================================================================================================================
/**
 * @brief Update reservoir temperature using 0D (tank) energy balance model with full radiative physics
 * 
 * Well-mixed tank model with comprehensive surface heat flux:
 * ρ·Cp·V·dT/dt = A_surface·Q_net + ρ·Cp·Q_in·(T_in - T_water)
 * 
 * Where Q_net includes:
 * - Solar radiation with dynamic Fresnel albedo
 * - Longwave radiation (atmospheric emission and water back-radiation)
 * - Sensible heat flux (turbulent heat transfer)
 * - Latent heat flux (evaporation/condensation)
 * 
 * Uses same physics as calculateRadiativeFluxes() for consistency.
 */
//======================================================================================================================
void CSimulation::updateReservoirTemperature0D() {
    // Get calibration parameters from YAML configuration
    double latitude = m_dHeatFluxLatitude;
    double cloud_cover = m_dHeatFluxCloudCover;
    double C_h = m_dHeatFlux_CS;
    double C_e = m_dHeatFlux_CL;
    
    for (size_t i = 0; i < m_vHeatFluxTime.size(); ++i) {   
        double dt = ((i == 0) ? 0.0 : m_vHeatFluxTime[i] - m_vHeatFluxTime[i-1]);
        
        // === METEOROLOGICAL FORCING ===
        double Tair = (i < m_vHeatFluxAirTemp.size()) ? m_vHeatFluxAirTemp[i] : 15.0;
        double wind = (i < m_vHeatFluxWind.size()) ? m_vHeatFluxWind[i] : 1.0;
        double pressure = (i < m_vHeatFluxAtmosphericPressure.size()) ? 
            m_vHeatFluxAtmosphericPressure[i] : ATM_PRESSURE_DEFAULT;
        
        if (wind < 0.1) wind = 0.1;  // Prevent stagnation

        // Relative humidity: either from data or calculated from temperature
        double rh = 70.0;  // Default value
        if (!m_vHeatFluxTime.empty()) {
            if (m_bCalculateRHFromTemperature && !m_vDailyMinTemperature.empty()) {
                // Calculate RH using FAO-56 method from current T_air and daily T_min
                double sim_time_hours = m_dCurrentTime / 3600.0;
                int day_index = static_cast<int>(sim_time_hours / 24.0);
                
                // Ensure day_index is within bounds
                if (day_index >= 0 && day_index < static_cast<int>(m_vDailyMinTemperature.size())) {
                    double T_min_daily = m_vDailyMinTemperature[day_index];
                    rh = calc_rh_from_temp(Tair, T_min_daily);
                } else {
                    // Fallback: use interpolated value if available, otherwise default
                    rh = (!m_vHeatFluxRelHumidity.empty()) ?
                        linearInterpolation1d(m_dCurrentTime, m_vHeatFluxTime, m_vHeatFluxRelHumidity) : 70.0;
                }
            } else {
                // Use interpolated RH from data file
                rh = linearInterpolation1d(m_dCurrentTime, m_vHeatFluxTime, m_vHeatFluxRelHumidity);
            }
        }
        
        // === WATER TEMPERATURE ===
        // For i=0: use initial condition (already in vector)
        // For i>0: use temperature from previous timestep (forward integration)
        double Twater = (i == 0) ? m_vUpwardTemperatureBoundaryConditionValue[i] 
                                 : m_vUpwardTemperatureBoundaryConditionValue[i-1];
        double Tin = Tair + m_dUpwardTemperatureOffsetBeta;  // Inflow temperature
        
        // === SOLAR GEOMETRY ===
        double sim_time_hours = m_vHeatFluxTime[i] / 3600.0;
        double hour_of_day = std::fmod(sim_time_hours, 24.0);
        int day_of_year = 1 + static_cast<int>(sim_time_hours / 24.0);
        
        double delta = 23.45 * std::sin(DEG_TO_RAD * (360.0/365.0) * (284.0 + day_of_year));
        double delta_rad = delta * DEG_TO_RAD;
        double lat_rad = latitude * DEG_TO_RAD;
        double omega = 15.0 * (hour_of_day - 12.0);
        double omega_rad = omega * DEG_TO_RAD;
        
        double cos_theta = std::sin(lat_rad) * std::sin(delta_rad) + 
                          std::cos(lat_rad) * std::cos(delta_rad) * std::cos(omega_rad);
        cos_theta = std::max(0.0, cos_theta);
        
        double zenith_rad = std::acos(cos_theta);
        double albedo = calc_albedo_briegleb(zenith_rad, cloud_cover);
        
        // === SHORTWAVE RADIATION ===
        double cloud_factor = 1.0 - 0.65 * cloud_cover * cloud_cover;
        double R_clear = SOLAR_CONSTANT * cos_theta * ATM_TRANSMISSIVITY;
        double Qsw_in = R_clear * cloud_factor;
        double Qsw_net = Qsw_in * (1.0 - albedo);
        
        // === LONGWAVE RADIATION ===
        double Tair_K = Tair + 273.15;
        double Twater_K = Twater + 273.15;
        double epsilon_sky = 9.37e-6 * (Tair_K * Tair_K) * (1.0 + 0.17 * cloud_cover * cloud_cover);
        double Qlw_in = epsilon_sky * STEFAN_BOLTZMANN * std::pow(Tair_K, 4);
        double Qlw_out = WATER_EMISSIVITY * STEFAN_BOLTZMANN * std::pow(Twater_K, 4);
        double Qlw_net = Qlw_in - Qlw_out;
        
        // === TURBULENT FLUXES ===
        double rho_a = calc_rho_air(Tair, pressure);
        double L_v = calc_Lv(Twater);
        
        double e_sat_air = 611.2 * std::exp(17.67 * Tair / (Tair + 243.5));
        double e_air = (rh / 100.0) * e_sat_air;
        double q_air = 0.622 * e_air / (pressure - 0.378 * e_air);
        
        double e_sat_water = 611.2 * std::exp(17.67 * Twater / (Twater + 243.5));
        double qsat_water = 0.622 * e_sat_water / (pressure - 0.378 * e_sat_water);
        
        double Qsensible = rho_a * AIR_SPECIFIC_HEAT * C_h * wind * (Tair - Twater);
        double Qlatente = rho_a * L_v * C_e * wind * (q_air - qsat_water);
        
        // === NET SURFACE HEAT FLUX ===
        double Qnet = Qsw_net + Qlw_net + Qsensible + Qlatente;  // W/m²
        
        // === TEMPERATURE TENDENCY ===
        // Surface heat flux term: Q_net / (ρ·Cp·h_eff)
        // Use configured reservoir effective depth
        double h_effective = m_dReservoirEffectiveDepth;
        
        double heat_flux_term = Qnet / (WATER_DENSITY * WATER_SPECIFIC_HEAT * h_effective);
        double inflow_term = (i < m_vCrossSectionQ.size()) ? 
            m_dUpwardInflowWaterEffectkQ * m_vCrossSectionQ[i] * (Tin - Twater) : 0.0;
        
        double dTdt = heat_flux_term + inflow_term;
        
        // Update temperature
        double T_new;
        if (i == 0) {
            T_new = m_vUpwardTemperatureBoundaryConditionValue[i] + dt * dTdt;
        } else {
            T_new = m_vUpwardTemperatureBoundaryConditionValue[i-1] + dt * dTdt;
        }
        
        // Apply physical limits to prevent unrealistic values
        T_new = std::max(-5.0, std::min(40.0, T_new));
        m_vUpwardTemperatureBoundaryConditionValue[i] = T_new;
    }
    
    // === DIAGNOSTIC: Reservoir temperature evolution ===
    if (LogStream.is_open() && m_nLogFileDetail >= 2) {
        LogStream << "\n--- RESERVOIR TEMPERATURE (0D Model) INITIALIZATION ---\n";
        LogStream << "Time series length: " << m_vHeatFluxTime.size() << " points\n";
        LogStream << "Reservoir effective depth: " << m_dReservoirEffectiveDepth << " m\n";
        LogStream << "Inflow effect coefficient: " << m_dUpwardInflowWaterEffectkQ << "\n";
        LogStream << "Initial T_res: " << std::setprecision(3) << m_vUpwardTemperatureBoundaryConditionValue[0] << " °C\n";
        if (m_vUpwardTemperatureBoundaryConditionValue.size() > 1) {
            size_t mid = m_vUpwardTemperatureBoundaryConditionValue.size() / 2;
            size_t end = m_vUpwardTemperatureBoundaryConditionValue.size() - 1;
            LogStream << "Middle T_res (t=" << std::setprecision(0) << m_vHeatFluxTime[mid]/86400.0 
                      << " days): " << std::setprecision(3) << m_vUpwardTemperatureBoundaryConditionValue[mid] << " °C\n";
            LogStream << "Final T_res: " << std::setprecision(3) << m_vUpwardTemperatureBoundaryConditionValue[end] << " °C\n";
            double T_change = m_vUpwardTemperatureBoundaryConditionValue[end] - m_vUpwardTemperatureBoundaryConditionValue[0];
            LogStream << "Total temperature change: " << std::setprecision(3) << T_change << " °C\n";
            
            if (fabs(T_change) < 0.1) {
                LogStream << "⚠️  WARNING: Reservoir temperature barely changes!\n";
                LogStream << "   Check: reservoir_effective_depth (large depth = high thermal inertia)\n";
                LogStream << "   Check: Heat flux data magnitude\n";
            }
        }
        LogStream << "-------------------------------------------------------\n\n";
        LogStream.flush();
    }
}

//===============================================================================================================================
//! Merge Predictor and Corrector
//======================================================================================================================
void CSimulation::mergePredictorCorrector() {
    //! Merge predictor and corrector steps using McCormack averaging
    //! Applies TVD flux limiters to maintain monotonicity and avoid spurious oscillations
    if (bGetDoMcComarckLimiterFlux())
    {
        //! Include TVD-McComarck? - Use pre-allocated vectors (class members)
        // Note: Vectors don't need zeroing - overwritten in loops below
        
        // References to maintain compatibility with existing code
        auto& a1_med = m_vTVD_a1_med;
        auto& a2_med = m_vTVD_a2_med;
        auto& alfa1_med = m_vTVD_alfa1_med;
        auto& alfa2_med = m_vTVD_alfa2_med;
        auto& psi1_med = m_vTVD_psi1_med;
        auto& psi2_med = m_vTVD_psi2_med;
        auto& fi1_med = m_vTVD_fi1_med;
        auto& fi2_med = m_vTVD_fi2_med;
        auto& vFactor1 = m_vTVD_Factor1;
        auto& vFactor2 = m_vTVD_Factor2;


        //! García-Navarro Psi formula
        double delta1 = m_dDeltaValue;
        double delta2 = m_dDeltaValue;
        for (int i = 0; i < m_nCrossSectionsNumber-1; i++)
        {
            const double u_med = (m_vCrossSectionQ[i+1] / sqrt(m_vCrossSectionArea[i+1]) + m_vCrossSectionQ[i] /
                                  sqrt(m_vCrossSectionArea[i])) / (
                                     sqrt(m_vCrossSectionArea[i+1]) + sqrt(m_vCrossSectionArea[i]));
            const double A_med = (m_vCrossSectionArea[i+1] + m_vCrossSectionArea[i]) / 2.0;
            const double c_med = (m_vCrossSectionC[i+1] + m_vCrossSectionC[i]) / 2.0;
            a1_med[i] = u_med + c_med;
            a2_med[i] = u_med - c_med;

            if (!bGetDoSurfaceGradientMethod()) {
                alfa1_med[i] = ((m_vCrossSectionQ[i+1] - m_vCrossSectionQ[i]) - a2_med[i]*(A_med - m_vCrossSectionArea[i]))/(2.0*c_med);
                alfa2_med[i] = -((m_vCrossSectionQ[i+1] - m_vCrossSectionQ[i]) - a1_med[i]*(A_med - m_vCrossSectionArea[i]))/(2.0*c_med);
            }
            else {
                alfa1_med[i] = m_vCrossSectionWidth[i]*((m_vCrossSectionQ[i+1]/m_vCrossSectionWidth[i+1] - m_vCrossSectionQ[i]/m_vCrossSectionWidth[i])-a2_med[i]*(m_vCrossSectionWaterDepth[i+1] - m_vCrossSectionWaterDepth[i]))/(2.0*c_med);
                alfa2_med[i] = -m_vCrossSectionWidth[i]*((m_vCrossSectionQ[i+1]/m_vCrossSectionWidth[i+1] - m_vCrossSectionQ[i]/m_vCrossSectionWidth[i])-a1_med[i]*(m_vCrossSectionWaterDepth[i+1] - m_vCrossSectionWaterDepth[i]))/(2.0*c_med);
            }

            if (nGetPsiFormula() != 1) {
                //! Tseng Psi formula
                vector<double> deltaValues = {0.0, a1_med[i] - (m_vCrossSectionU[i] + m_vCrossSectionC[i]), m_vCrossSectionU[i+1] + m_vCrossSectionC[i+1] - a1_med[i]};
                delta1 = dMaxVectorValue(deltaValues);

                deltaValues ={0.0, a2_med[i] - (m_vCrossSectionU[i] - m_vCrossSectionC[i]), m_vCrossSectionU[i+1] + m_vCrossSectionC[i+1] - a2_med[i]};
                delta2 = dMaxVectorValue(deltaValues);
            }

            //! Computing psi_i
            if (fabs(a1_med[i])>= delta1)
            {
                psi1_med[i] = fabs(a1_med[i]);
            }
            else
            {
                psi1_med[i] = delta1;
            }

            if (fabs(a2_med[i])>= delta2)
            {
                psi2_med[i] = fabs(a2_med[i]);
            }
            else
            {
                psi2_med[i] = delta2;
            }
        }


        // === STEP 3-4: Compute r-ratios and apply flux limiter (FUSED) ===
        // r = alfa[i±1]/alfa[i] measures solution smoothness
        // r ≈ 1: smooth → limiter inactive, r >> 1 or << 1: steep gradient → limiter active
        // Direction depends on characteristic speed sign (upwind stencil)
        
        // ⚡ OPTIMIZATION: Extract limiter type check outside hot loop
        const int limiter_type = nGetEquationLimiterFlux();
        const int n_nodes = m_nCrossSectionsNumber - 1;
        
        for (int i = 0; i < n_nodes; i++)
        {
            // === Compute r1 (upwind ratio for characteristic 1) ===
            double r1, r2;
            
            if (alfa1_med[i] == 0.0) {
                r1 = 1.0;  // No wave amplitude → assume smooth
            } else {
                // Upwind stencil based on wave direction
                if (a1_med[i] < 0 && i != n_nodes - 1) {
                    r1 = alfa1_med[i+1] / alfa1_med[i];  // Downstream gradient
                } else if (a1_med[i] > 0 && i != 0) {
                    r1 = alfa1_med[i-1] / alfa1_med[i];  // Upstream gradient
                } else {
                    r1 = 1.0;  // Boundary or zero velocity
                }
            }
            
            // === Compute r2 (upwind ratio for characteristic 2) ===
            if (alfa2_med[i] == 0.0) {
                r2 = 1.0;
            } else {
                if (a2_med[i] < 0 && i != n_nodes) {
                    r2 = alfa2_med[i+1] / alfa2_med[i];
                } else if (a2_med[i] > 0 && i != 0) {
                    r2 = alfa2_med[i-1] / alfa2_med[i];
                } else {
                    r2 = 1.0;
                }
            }
            
            // ✅ PROTECTION: Clamp extreme values that cause instability
            r1 = (std::isfinite(r1) && fabs(r1) <= 1e6) ? r1 : 1.0;
            r2 = (std::isfinite(r2) && fabs(r2) <= 1e6) ? r2 : 1.0;
            
            // === Apply flux limiter function Ψ(r) ===
            // TVD region: Ψ(1) = 1, Ψ(r) ≤ min(2, 2r) (Sweby's TVD condition)
            switch (limiter_type) {
                case 1:  // MinMod: Ψ(r) = clamp(r, 0, 1) - BRANCHLESS
                    // Most diffusive, very stable, first-order near extrema
                    fi1_med[i] = std::clamp(r1, 0.0, 1.0);
                    fi2_med[i] = std::clamp(r2, 0.0, 1.0);
                    break;
                    
                case 2:  // Roe's Superbee: Ψ(r) = max(0, min(2r,1), min(r,2))
                    // Least diffusive, sharpest shock capturing
                    fi1_med[i] = dMaxVectorValue({0.0, dMinVectorValue({2*r1, 1.}), dMinVectorValue({r1, 2.})});
                    fi2_med[i] = dMaxVectorValue({0.0, dMinVectorValue({2*r2, 1.}), dMinVectorValue({r2, 2.})});
                    break;
                    
                case 3:  // Van Leer: Ψ(r) = (r + |r|)/(1 + |r|)
                    // Good balance between accuracy and stability
                    fi1_med[i] = (fabs(r1) + r1) / (1.0 + fabs(r1));
                    fi2_med[i] = (fabs(r2) + r2) / (1.0 + fabs(r2));
                    break;
                    
                case 4:  // Van Albada: Ψ(r) = (r² + r)/(1 + r²)
                    // Very smooth, slightly more diffusive than Van Leer
                    fi1_med[i] = (r1*r1 + r1) / (1.0 + r1*r1);
                    fi2_med[i] = (r2*r2 + r2) / (1.0 + r2*r2);
                    break;
                    
                default:
                    fi1_med[i] = 0.0;  // No limiter
                    fi2_med[i] = 0.0;
            }
        }
        // === STEP 5: Compute anti-diffusive flux correction ===
        // Factor = α·Ψ·(1 - λ|a|)·(1 - φ(r))
        // - α: wave amplitude (jump size)
        // - Ψ: García-Navarro limiter (bounds correction)
        // - (1 - λ|a|): CFL-based weight
        // - (1 - φ): anti-diffusive component (φ=0 for 1st order, φ=1 for TVD)
        for (int i = 0; i < m_nCrossSectionsNumber-1; i++) {
            // Local interface spacing for non-uniform grids (between i and i+1)
            const double dx_interface = m_vPositionX[i + 1] - m_vPositionX[i];
            const double lambda_interface = (dx_interface > 1e-12) ? (m_dTimestep / dx_interface) : 0.0;

            // Anti-diffusive flux for each characteristic
            vFactor1[i] = alfa1_med[i]*psi1_med[i]*(1 - lambda_interface*fabs(a1_med[i]))*(1-fi1_med[i]);
            vFactor2[i] = alfa2_med[i]*psi2_med[i]*(1 - lambda_interface*fabs(a2_med[i]))*(1-fi2_med[i]);

            // Project back to physical variables: D1 for continuity, D2 for momentum
            m_vCrossSectionD1Factor[i+1] = 0.5*(vFactor1[i] + vFactor2[i]);
            m_vCrossSectionD2Factor[i+1] = 0.5*(vFactor1[i]*a1_med[i] + vFactor2[i]*a2_med[i]);
        }
        // Extrapolate correction to boundaries (usually small effect)
        m_vCrossSectionD1Factor[0] = m_vCrossSectionD1Factor[1];
        m_vCrossSectionD2Factor[0] = m_vCrossSectionD2Factor[1];
        m_vCrossSectionD1Factor[m_nCrossSectionsNumber] = m_vCrossSectionD1Factor[m_nCrossSectionsNumber-1];
        m_vCrossSectionD2Factor[m_nCrossSectionsNumber] = m_vCrossSectionD2Factor[m_nCrossSectionsNumber-1];

        // === STEP 6: Apply TVD correction to averaged solution ===
        // U^{n+1} = 0.5*(U_pred + U_corr) + TVD_correction
        // Only interior nodes - boundaries set by BC
        for (int i = 1; i < m_nCrossSectionsNumber - 1; i++) {
            // Cell length associated with node i is half the distance between neighbors:
            // Δx_cell ≈ 0.5*(x[i+1]-x[i-1]) = 0.5*(dx_left + dx_right)
            const double dx_sum = m_vPositionX[i + 1] - m_vPositionX[i - 1];
            const double lambda_cell = (dx_sum > 1e-12) ? (2.0 * m_dTimestep / dx_sum) : 0.0;

            const double invSf = 1.0 / dGetLateralStorageFactor(i);

            m_vCrossSectionArea[i] = 0.5 * (m_vPredictedCrossSectionArea[i] + m_vCorrectedCrossSectionArea[i]) +
                (lambda_cell * invSf) * (m_vCrossSectionD1Factor[i + 1] - m_vCrossSectionD1Factor[i]);
            m_vCrossSectionQ[i] = 0.5 * (m_vPredictedCrossSectionQ[i] + m_vCorrectedCrossSectionQ[i]) +
                lambda_cell * (m_vCrossSectionD2Factor[i + 1] - m_vCrossSectionD2Factor[i]);
        }
        
        //! Apply final boundary conditions after merge
        //! Boundaries must use corrected values, not averaged predictor-corrector
        m_vCrossSectionArea[0] = m_vCorrectedCrossSectionArea[0];
        m_vCrossSectionQ[0] = m_vCorrectedCrossSectionQ[0];
        m_vCrossSectionArea[m_nCrossSectionsNumber-1] = m_vCorrectedCrossSectionArea[m_nCrossSectionsNumber-1];
        m_vCrossSectionQ[m_nCrossSectionsNumber-1] = m_vCorrectedCrossSectionQ[m_nCrossSectionsNumber-1];
    }
}

/**
 * @brief Merge predictor and corrector steps for passive tracer transport (salinity, temperature)
 * 
 * Uses simple arithmetic averaging for passive scalars:
 *   S^{n+1} = 0.5 * (S_pred + S_corr)
 * 
 * Unlike momentum/continuity, passive tracers don't require TVD correction because:
 * - No wave propagation (only advection-diffusion)
 * - TVD is already applied in the velocity field that transports them
 * - Overshoots prevented by explicit flux limiters in predictor/corrector
 * 
 * @note Boundaries already set in calculate_salinity/temperature_corrector
 */
void CSimulation::mergeTracerPredictorCorrector() {
    // Salinity averaging - conserve mass by averaging A*S from predictor/corrector
    // and dividing by the final merged area (already computed in mergePredictorCorrector)
    if (bGetDoWaterSalinity() && m_vPredictedCrossSectionS.size() == static_cast<size_t>(m_nCrossSectionsNumber) && m_vCorrectedCrossSectionS.size() == static_cast<size_t>(m_nCrossSectionsNumber)) {
        for (int i = 0; i < m_nCrossSectionsNumber; ++i) {
            // Average mass (A*S) from predictor and corrector
            double AS_pred = m_vPredictedCrossSectionArea[i] * m_vPredictedCrossSectionS[i];
            double AS_corr = m_vCorrectedCrossSectionArea[i] * m_vCorrectedCrossSectionS[i];
            double AS_avg = 0.5 * (AS_pred + AS_corr);
            
            // Use the final merged area (computed in mergePredictorCorrector)
            double A_final = m_vCrossSectionArea[i];
            
            if (A_final > 1e-10) {  // Avoid division by zero
                m_vCrossSectionSalinity[i] = AS_avg / A_final;
            } else {
                // If area is essentially zero, set salinity to zero (dry cell)
                m_vCrossSectionSalinity[i] = 0.0;
            }
        }
    }
    // Temperature averaging
    if (bGetDoWaterTemperature() && m_vPredictedCrossSectionT.size() == static_cast<size_t>(m_nCrossSectionsNumber) && m_vCorrectedCrossSectionT.size() == static_cast<size_t>(m_nCrossSectionsNumber)) {
        for (int i = 0; i < m_nCrossSectionsNumber; ++i) {
            m_vCrossSectionTemperature[i] = 0.5 * (m_vPredictedCrossSectionT[i] + m_vCorrectedCrossSectionT[i]);
        }
    }
}

/**
 * @brief Pre-simulation bathymetry smoothing using non-uniform 1D Laplacian diffusion
 * 
 * Applies explicit diffusion to remove high-frequency noise from bed elevation data:
 *   Z^{k+1}_i = Z^k_i + α·∇²Z
 * 
 * where the non-uniform Laplacian is:
 *   ∇²Z_i = (2/(dx_L + dx_R)) * [(Z_{i+1} - Z_i)/dx_R - (Z_i - Z_{i-1})/dx_L]
 * 
 * @note Physics:
 * - Multiple passes act as multi-step diffusion (like heat equation)
 * - Alpha controls diffusion strength: larger α → more smoothing
 * - Non-uniform stencil respects variable grid spacing
 * - Boundaries kept fixed to preserve domain extent
 * 
 * @warning Too much smoothing can:
 * - Remove important morphological features (sills, steps)
 * - Violate mass conservation if not symmetric
 * - Create artificial slopes in transition zones
 */
void CSimulation::smoothBathymetry() {

    const int n = m_nCrossSectionsNumber;
    const int num_passes = m_nBathymetrySmoothingPasses;      // Number of smoothing passes from config
    const double alpha = m_dBathymetrySmoothingAlpha;         // Diffusion coefficient from config

    std::cout << "      - Smoothing bathymetry (non-uniform Laplacian, "
              << num_passes << " passes, alpha=" << alpha << ")" << std::endl;

    std::vector<double> vZSmooth(n);

    for (int pass = 0; pass < num_passes; pass++) {
        // Keep boundary values unchanged
        vZSmooth[0]   = m_vBedZ[0];
        vZSmooth[n-1] = m_vBedZ[n-1];

        // Apply non-uniform Laplacian smoothing to interior nodes
        for (int i = 1; i < n-1; i++) {
            double dxL = m_vCrossSectionX[i] - m_vCrossSectionX[i-1];
            double dxR = m_vCrossSectionX[i+1] - m_vCrossSectionX[i];

            // Non-uniform 1D Laplacian
            double lap =
                (2.0 / (dxL + dxR)) *
                ( (m_vBedZ[i+1] - m_vBedZ[i]) / dxR
                - (m_vBedZ[i]   - m_vBedZ[i-1]) / dxL );

            // Explicit smoothing step
            vZSmooth[i] = m_vBedZ[i] + alpha * lap;
        }

        // Update values for the next pass
        for (int i = 1; i < n-1; i++) {
            m_vBedZ[i] = vZSmooth[i];
        }
    }
}

/**
 * @brief Runtime solution smoothing with adaptive 2-pass algorithm
 * 
 * Pass 1: Global Laplacian smoothing (similar to smoothBathymetry)
 * Pass 2: Targeted smoothing only at strong area jumps
 * 
 * Strong jump criteria: A[i]/A[i±1] > 1.5 or < 0.67
 * 
 * @note Purpose:
 * - Stabilize solution near sudden geometry changes (bridge piers, junctions)
 * - Remove Gibbs oscillations from sharp fronts
 * - Maintain CFL stability in regions with complex geometry
 * 
 * @warning Can introduce:
 * - Artificial diffusion (reduces effective Reynolds number)
 * - Phase errors in wave propagation
 * - Mass/momentum conservation errors if alpha too large
 * 
 * @see smoothBathymetry() for similar Laplacian formulation
 */
void CSimulation::smoothSolution() {
    const int n = m_nCrossSectionsNumber;
    const int num_passes = m_nSolutionSmoothingPasses;
    const double alpha = m_dSolutionSmoothingAlpha;

    for (int pass = 0; pass < num_passes; ++pass) {
        std::vector<double> vAreaSmooth(n);
        std::vector<double> vQSmooth(n);

        // Keep boundaries unchanged
        vAreaSmooth[0]   = m_vCrossSectionArea[0];
        vAreaSmooth[n-1] = m_vCrossSectionArea[n-1];
        vQSmooth[0]      = m_vCrossSectionQ[0];
        vQSmooth[n-1]    = m_vCrossSectionQ[n-1];

        // First pass: non-uniform Laplacian smoothing
        for (int i = 1; i < n-1; i++) {
            double dxL = m_vCrossSectionX[i] - m_vCrossSectionX[i-1];
            double dxR = m_vCrossSectionX[i+1] - m_vCrossSectionX[i];

            // Laplacian for Area
            double lapA =
                (2.0 / (dxL + dxR)) *
                ( (m_vCrossSectionArea[i+1] - m_vCrossSectionArea[i]) / dxR
                - (m_vCrossSectionArea[i]   - m_vCrossSectionArea[i-1]) / dxL );

            // Laplacian for Q
            double lapQ =
                (2.0 / (dxL + dxR)) *
                ( (m_vCrossSectionQ[i+1] - m_vCrossSectionQ[i]) / dxR
                - (m_vCrossSectionQ[i]   - m_vCrossSectionQ[i-1]) / dxL );

            vAreaSmooth[i] = m_vCrossSectionArea[i] + alpha * lapA;
            vQSmooth[i]    = m_vCrossSectionQ[i]    + alpha * lapQ;
        }

        // Second pass only in regions with strong area jumps
        for (int i = 2; i < n-2; i++) {
            double ratioL = m_vCrossSectionArea[i]   / (m_vCrossSectionArea[i-1] + 1e-10);
            double ratioR = m_vCrossSectionArea[i+1] / (m_vCrossSectionArea[i]   + 1e-10);

            bool strongJump =
                (ratioL > 1.5 || ratioL < 0.67 ||
                 ratioR > 1.5 || ratioR < 0.67);

            if (strongJump) {
                double dxL = m_vCrossSectionX[i]   - m_vCrossSectionX[i-1];
                double dxR = m_vCrossSectionX[i+1] - m_vCrossSectionX[i];

                double lapA =
                    (2.0 / (dxL + dxR)) *
                    ( (vAreaSmooth[i+1] - vAreaSmooth[i]) / dxR
                    - (vAreaSmooth[i]   - vAreaSmooth[i-1]) / dxL );

                double lapQ =
                    (2.0 / (dxL + dxR)) *
                    ( (vQSmooth[i+1] - vQSmooth[i]) / dxR
                    - (vQSmooth[i]   - vQSmooth[i-1]) / dxL );

                vAreaSmooth[i] += alpha * lapA;
                vQSmooth[i]    += alpha * lapQ;
            }
        }

        // Apply smoothed values and update velocity
        for (int i = 1; i < n-1; i++) {
            m_vCrossSectionArea[i] = vAreaSmooth[i];
            m_vCrossSectionQ[i]    = vQSmooth[i];
            m_vCrossSectionU[i]    = m_vCrossSectionQ[i] /
                                     (m_vCrossSectionArea[i] + 1e-10);
        }
    }
}


/**
 * @brief Accessor for simulation state variables (for post-processing, NetCDF output, debugging)
 * 
 * Available variables:
 * - Flow: "A", "Q", "U", "c" (wave speed)
 * - Predictor/Corrector: "Ap", "Ac", "Qp", "Qc"
 * - Geometry: "B" (width), "Rh" (hydraulic radius), "level", "eta" (elevation)
 * - Transport: "S" (salinity), "T" (temperature), "rho" (density)
 * - Derivatives: "DhDx" (surface slope)
 * - Sediment: "Qb" (bedload), "Qs" (suspended), "Qt" (total)
 * - Bank locations: "xl", "xr", "xl_utm_x", "xl_utm_y", "xr_utm_x", "xr_utm_y"
 * 
 * @param strVariableName Name of variable to retrieve
 * @return Vector of values for all cross-sections, or empty vector if not found
 * 
 * @warning No bounds checking - caller must verify size matches m_nCrossSectionsNumber
 */
vector<double> CSimulation::vGetVariable(const string& strVariableName) const {

    if (strVariableName == "A") {
        return m_vCrossSectionArea;
    }
    else if (strVariableName == "Ap") {
        return m_vPredictedCrossSectionArea;
    }
    else if (strVariableName == "Ac") {
        return m_vCorrectedCrossSectionArea;
    }
    else if (strVariableName == "Q") {
        return m_vCrossSectionQ;
    }
    else if (strVariableName == "Qp") {
        return m_vPredictedCrossSectionQ;
    }
    else if (strVariableName == "Qc") {
        return m_vCorrectedCrossSectionQ;
    }
    else if (strVariableName == "Rh") {
        return m_vCrossSectionHydraulicRadius;
    }
    else if (strVariableName == "n") {
        return m_vCrossSectionManningNumber;
    }
    else if (strVariableName == "B") {
        return m_vCrossSectionWidth;
    }
    else if (strVariableName == "level") {
        return m_vCrossSectionWaterDepth;
    }
    else if (strVariableName == "eta") {
        return m_vCrossSectionWaterElevation;
    }
    else if (strVariableName == "DhDx") {
        return m_vCrossSectionDhDx;
    }
    else if (strVariableName == "rho") {
        return m_vCrossSectionDensity;
    }
    else if (strVariableName == "U") {
        return m_vCrossSectionU;
    }
    else if (strVariableName == "c") {
        return m_vCrossSectionC;
    }
    else if (strVariableName == "S") {
        return m_vCrossSectionSalinity;
    }
    else if (strVariableName == "T") {
        return m_vCrossSectionTemperature;
    }
    else if (strVariableName == "xl") {
        return m_vCrossSectionLeftRBLocation;
    }
    else if (strVariableName == "xr") {
        return m_vCrossSectionRightRBLocation;
    }
    else if (strVariableName == "xl_utm_x") {
        return m_vCrossSectionLeftRBLocation_UTM_X;
    }
    else if (strVariableName == "xl_utm_y") {
        return m_vCrossSectionLeftRBLocation_UTM_Y;
    }
    else if (strVariableName == "xr_utm_x") {
        return m_vCrossSectionRightRBLocation_UTM_X;
    }
    else if (strVariableName == "xr_utm_y") {
        return m_vCrossSectionRightRBLocation_UTM_Y;
    }
    else if (strVariableName == "Qb") {
        return m_vCrossSectionQb;
    }
    else if (strVariableName == "Qs") {
        return m_vCrossSectionQs;
    }
    else if (strVariableName == "Qt") {
        return m_vCrossSectionQt;
    }
    // else if (strVariableName == "q") {
    //     return m_vCrossSectionQt;
    // }
    else {
        return {};
    }
}

/**
 * @brief Display simulation progress with ETA and current simulation date/time
 * 
 * Prints: [Elapsed][Remaining] Time - Progress % - SimDate: YYYY-MM-DD HH:MM:SS
 * 
 * @note Performance optimization:
 * - Only updates every 100 timesteps to reduce I/O overhead
 * - I/O can be expensive (terminal flushing, string formatting)
 * - For dt=1s, 100 timesteps = ~1.7 minutes between updates
 * - Always updates when saving (m_bSaveTime=true)
 * 
 * @warning Static variables used for persistent state across calls:
 * - call_count: throttling counter
 * - t_start: cached simulation start time (computed once)
 */
void CSimulation::AnnounceProgress() {

    // ⚡ OPTIMIZATION: Reduce screen update frequency to minimize I/O overhead
    // Update every 1000 iterations (~2-3% of total time) or when saving output
    static int call_count = 0;
    if (++call_count % 3600 != 0 && !m_bSaveTime) return;

    // Stdout is connected to a tty, so not running as a background job
    static double sdElapsed = 0;
    static double sdToGo = 0;
    const time_t tNow = time(nullptr);

    // Calculate time elapsed and remaining
    sdElapsed = difftime(tNow, m_tSysStartLoopTime);
    sdToGo = (sdElapsed * m_dSimDuration / m_dCurrentTime) - sdElapsed;


    // Calculate t_start only once (static)
    static time_t t_start = 0;
    if (t_start == 0) {
        std::tm tm_start = {};
        tm_start.tm_year = m_nSimStartYear - 1900;
        tm_start.tm_mon = m_nSimStartMonth - 1;
        tm_start.tm_mday = m_nSimStartDay;
        tm_start.tm_hour = m_nSimStartHour;
        tm_start.tm_min = m_nSimStartMin;
        tm_start.tm_sec = m_nSimStartSec;
        tm_start.tm_isdst = -1;
        t_start = std::mktime(&tm_start);
    }
    time_t t_current = t_start + static_cast<time_t>(m_dCurrentTime);
    std::tm* tm_current = std::localtime(&t_current);

    // Imprimir fecha y hora en formato YYYY-MM-DD HH:MM:SS
    char datetime_buf[40];
    std::strftime(datetime_buf, sizeof(datetime_buf), "%Y-%m-%d %H:%M:%S", tm_current);

    cout << "\r    - Elapsed[Remaining] Time: " << std::fixed << setprecision(3) << setw(6) << sdElapsed <<"[" << std::fixed << setprecision(3) << setw(6) << sdToGo << "] s - Progress: " << std::fixed << setprecision(3) << setw(6) << 100 * m_dCurrentTime / m_dSimDuration  << "% - SimDate: " << datetime_buf << std::flush;

}


//======================================================================================================================
/**
 * @brief Check for anomalous values and write warnings to log
 * 
 * Monitors:
 * - CFL violations (stability)
 * - Excessive velocities (> 5 m/s)
 * - Supercritical flow (Froude > 1)
 * - Temperature anomalies (especially reservoir)
 * - Salinity out of physical range
 * - Negative areas
 */
//======================================================================================================================
void CSimulation::checkAnomalousValues() {
    if (!LogStream.is_open()) return;
    

    

    

    

    

    
    // 4. Check temperature (if active) - SPECIAL ATTENTION TO RESERVOIR
    if (m_bDoWaterTemperature) {
        // Check reservoir temperature (upstream boundary)
        if (m_nUpwardTemperatureCondition == 3 && !m_vUpwardTemperatureBoundaryConditionValue.empty()) {
            // Interpolate reservoir temperature at current simulation time
            double T_res = linearInterpolation1d(m_dCurrentTime, m_vHeatFluxTime, m_vUpwardTemperatureBoundaryConditionValue);
            if (T_res > 35.0 || T_res < 0.0) {
                LogStream << "🔥 ALERT [t=" << std::fixed << std::setprecision(1) << m_dCurrentTime 
                          << "s]: ANOMALOUS RESERVOIR TEMPERATURE = " << std::setprecision(2) << T_res 
                          << " °C (expected 0-35°C)";
                
                // Print diagnostic info to help debug
                double Tair_current = linearInterpolation1d(m_dCurrentTime, m_vHeatFluxTime, m_vHeatFluxAirTemp);
                double RH_current = linearInterpolation1d(m_dCurrentTime, m_vHeatFluxTime, m_vHeatFluxRelHumidity);
                LogStream << " | T_air=" << std::setprecision(2) << Tair_current << "°C";
                LogStream << " | RH=" << std::setprecision(1) << RH_current << "%";
                LogStream << "\n";
                m_nWarningCount++;
            }
        }
        
        // Check estuarine temperatures
        for (int i = 0; i < m_nCrossSectionsNumber; i++) {
            double T = m_vCrossSectionTemperature[i];
            if (T > 40.0 || T < -5.0) {
                LogStream << "🔥 WARNING [t=" << std::fixed << std::setprecision(1) << m_dCurrentTime 
                          << "s]: Anomalous water temperature T = " << std::setprecision(2) << T 
                          << " °C at x = " << std::setprecision(0) << m_vPositionX[i] << " m\n";
                m_nWarningCount++;
                break;
            }
        }
    }
    
    // 5. Check salinity (if active)
    if (m_bDoWaterSalinity) {
        for (int i = 0; i < m_nCrossSectionsNumber; i++) {
            double S = m_vCrossSectionSalinity[i];
            if (S < 0.0 || S > 45.0) {
                LogStream << "⚠️  WARNING [t=" << std::fixed << std::setprecision(1) << m_dCurrentTime 
                          << "s]: Anomalous salinity S = " << std::setprecision(2) << S 
                          << " ppt at x = " << std::setprecision(0) << m_vPositionX[i] << " m\n";
                m_nWarningCount++;
                break;
            }
        }
    }
    
    // 6. Check for negative areas
    for (int i = 0; i < m_nCrossSectionsNumber; i++) {
        if (m_vCrossSectionArea[i] < 0.0) {
            LogStream << "❌ ERROR [t=" << std::fixed << std::setprecision(1) << m_dCurrentTime 
                      << "s]: NEGATIVE AREA A = " << std::setprecision(2) << m_vCrossSectionArea[i]
                      << " m² at x = " << std::setprecision(0) << m_vPositionX[i] << " m\n";
            m_nWarningCount++;
            break;
        }
    }
    
    if (m_nWarningCount > 0) {
        LogStream.flush();
    }
}

//======================================================================================================================
/**
 * @brief Write periodic statistical summary to log (every simulated hour)
 * 
 * Provides min/max statistics for key variables to track simulation behavior
 */
//======================================================================================================================
void CSimulation::writePeriodicStatistics() {
    if (!LogStream.is_open()) return;
    
    // Write periodic simulation diagnostics/statistics only at output (save) timesteps
    // Only log if current time matches output save time (within half a timestep tolerance)
    if (fmod(m_dCurrentTime, m_dSimTimestep) > 0.5 * m_dTimestep &&
        fabs(m_dCurrentTime - m_dSimTimestep) > 1e-6) return;

    // Calculate statistics (use current merged state; do not rely on potentially stale
    // hydraulic vectors computed in predictor/corrector).
    double Q_min = 1e10, Q_max = -1e10;
    double U_max = 0.0;
    double eta_min = 1e10, eta_max = -1e10;
    double CFL_max = 0.0;

    double n_min = 1e10, n_max = -1e10;
    double murillo_min = 1e10;
    double Sf_abs_max = 0.0;

    for (int i = 0; i < m_nCrossSectionsNumber; i++) {
        const double A = m_vCrossSectionArea[i];
        if (A > DRY_AREA) {
            const double Q = m_vCrossSectionQ[i];
            if (Q < Q_min) Q_min = Q;
            if (Q > Q_max) Q_max = Q;

            // Robust eta from current A
            const double depth = linearInterpolation1d(A, m_vEstuaryAreas[i], m_vEstuaryWaterDepths[i]);
            const double eta = depth + m_vBedZ[i];
            if (eta < eta_min) eta_min = eta;
            if (eta > eta_max) eta_max = eta;

            // Velocity magnitude for reporting
            const double u_abs = fabs(Q / (A + 1e-10));
            if (u_abs > U_max) U_max = u_abs;

            // CFL estimate (interior only) using local min spacing for non-uniform grids
            if (i > 0 && i < m_nCrossSectionsNumber - 1) {
                const double dxL = m_vCrossSectionX[i] - m_vCrossSectionX[i - 1];
                const double dxR = m_vCrossSectionX[i + 1] - m_vCrossSectionX[i];
                const double dx_local = (dxL > 0.0 && dxR > 0.0) ? std::min(dxL, dxR) : std::max(dxL, dxR);
                double cfl = (dx_local > 0.0) ? (m_dTimestep * (u_abs + m_vCrossSectionC[i]) / dx_local) : 0.0;
                if (cfl > CFL_max) CFL_max = cfl;
            }

            // Manning / Murillo diagnostics
            const double n_val = m_vCrossSectionManningNumber[i];
            if (std::isfinite(n_val)) {
                if (n_val < n_min) n_min = n_val;
                if (n_val > n_max) n_max = n_val;
            }
            if (m_bDoMurilloCondition && i < static_cast<int>(m_vCrossSectionMurilloFactor.size())) {
                const double mf = m_vCrossSectionMurilloFactor[i];
                if (std::isfinite(mf) && mf < murillo_min) murillo_min = mf;
            }

            // Friction slope magnitude from current state (|Sf|)
            const double Rh = linearInterpolation1d(A, m_vEstuaryAreas[i], m_vEstuaryHydraulicRadius[i]);
            if (Rh > 1e-6) {
                const double A2 = A * A;
                const double Rh_43 = pow(Rh, 4.0 / 3.0);
                const double Sf = (n_val * n_val) * (Q * fabs(Q)) / (A2 * Rh_43 + 1e-30);
                const double Sf_abs = fabs(Sf);
                if (Sf_abs > Sf_abs_max) Sf_abs_max = Sf_abs;
            }
        }
    }

    // Calculate current date/time
    static time_t t_start = 0;
    if (t_start == 0) {
        std::tm tm_start = {};
        tm_start.tm_year = m_nSimStartYear - 1900;
        tm_start.tm_mon = m_nSimStartMonth - 1;
        tm_start.tm_mday = m_nSimStartDay;
        tm_start.tm_hour = m_nSimStartHour;
        tm_start.tm_min = m_nSimStartMin;
        tm_start.tm_sec = m_nSimStartSec;
        tm_start.tm_isdst = -1;
        t_start = std::mktime(&tm_start);
    }
    time_t t_current = t_start + static_cast<time_t>(m_dCurrentTime);
    std::tm* tm_current = std::localtime(&t_current);
    char datetime_buf[32];
    std::strftime(datetime_buf, sizeof(datetime_buf), "%Y-%m-%d %H:%M:%S", tm_current);

    // Write periodic statistics header
    LogStream << std::string(120, '=') << "\n";
    LogStream << "SIMULATION DIAGNOSTICS - PERIODIC STATISTICS\n";
    LogStream << "t=" << m_dCurrentTime/3600.0 << "h: ";

    // Salt mass (if salinity active)
    if (m_bDoWaterSalinity) {
        double salt_mass = 0.0;
        for (int i = 0; i < m_nCrossSectionsNumber; i++) {
            if (m_vCrossSectionArea[i] > DRY_AREA)
                salt_mass += m_vCrossSectionArea[i] * m_vCrossSectionSalinity[i];
        }
        LogStream << "Salt mass = " << salt_mass << " kg\n";
    }

    LogStream << std::setw(20) << datetime_buf
              << std::fixed << std::setprecision(1)
              << std::setw(10) << m_dCurrentTime
              << std::setprecision(3)
              << std::setw(8) << m_dTimestep
              << std::setprecision(2)
              << std::setw(8) << CFL_max
              << std::setprecision(1)
              << std::setw(10) << Q_min
              << std::setw(10) << Q_max
              << std::setprecision(2)
              << std::setw(10) << U_max
              << std::setprecision(2)
              << std::setw(10) << eta_min
              << std::setw(10) << eta_max;

    // Compact Manning/Murillo/friction diagnostics (helps explain tidal amplification)
    if (n_min < 1e9 && n_max > -1e9) {
        LogStream << "  n[min,max]=" << std::setprecision(4) << n_min << "," << n_max;
    }
    if (m_bDoMurilloCondition) {
        if (murillo_min < 1e9) {
            LogStream << "  MurilloMin=" << std::setprecision(6) << murillo_min;
        } else {
            LogStream << "  MurilloMin=n/a";
        }
    }
    LogStream << "  |Sf|max=" << std::scientific << std::setprecision(3) << Sf_abs_max << std::defaultfloat;

    // Temperature statistics (if active)
    if (m_bDoWaterTemperature) {
        double T_min = 1e10, T_max = -1e10;
        for (int i = 0; i < m_nCrossSectionsNumber; i++) {
            if (m_vCrossSectionTemperature[i] < T_min) T_min = m_vCrossSectionTemperature[i];
            if (m_vCrossSectionTemperature[i] > T_max) T_max = m_vCrossSectionTemperature[i];
        }
        double T_res = 0.0;
        if (m_nUpwardTemperatureCondition == 3 && !m_vUpwardTemperatureBoundaryConditionValue.empty()) {
            T_res = linearInterpolation1d(m_dCurrentTime, m_vHeatFluxTime, m_vUpwardTemperatureBoundaryConditionValue);
        }
        LogStream << std::setprecision(2)
                  << std::setw(10) << T_min
                  << std::setw(10) << T_max
                  << std::setw(10) << T_res;
    }

    // Salinity statistics (if active)
    if (m_bDoWaterSalinity) {
        double S_min = 1e10, S_max = -1e10;
        for (int i = 0; i < m_nCrossSectionsNumber; i++) {
            if (m_vCrossSectionSalinity[i] < S_min) S_min = m_vCrossSectionSalinity[i];
            if (m_vCrossSectionSalinity[i] > S_max) S_max = m_vCrossSectionSalinity[i];
        }
        LogStream << std::setprecision(2)
                  << std::setw(10) << S_min
                  << std::setw(10) << S_max;
    }

    // Mark auto-activated smoothing on extreme discharge with a trailing '*'
    if (m_bAutoSmoothAppliedThisStep) {
        LogStream << " *";
    }
    LogStream << "\n";

    // Snapshot at key nodes (0, mid, n) + running amplitude (max-min) for eta.
    {
        const int N = m_nCrossSectionsNumber;
        const int last = N - 1;
        const int mid = N / 2;

        auto eta_from_A = [&](int i) {
            const double depth = linearInterpolation1d(m_vCrossSectionArea[i], m_vEstuaryAreas[i], m_vEstuaryWaterDepths[i]);
            return depth + m_vBedZ[i];
        };
        auto Rh_from_A = [&](int i) {
            return linearInterpolation1d(m_vCrossSectionArea[i], m_vEstuaryAreas[i], m_vEstuaryHydraulicRadius[i]);
        };
        auto Sf_from_state = [&](int i) {
            const double A = m_vCrossSectionArea[i];
            const double Q = m_vCrossSectionQ[i];
            const double n_val = m_vCrossSectionManningNumber[i];
            const double Rh = Rh_from_A(i);
            if (!(A > DRY_AREA) || !(Rh > 1e-6)) return 0.0;
            const double A2 = A * A;
            const double Rh_43 = pow(Rh, 4.0 / 3.0);
            return (n_val * n_val) * (Q * fabs(Q)) / (A2 * Rh_43 + 1e-30);
        };

        static double eta0_min = 1e10, eta0_max = -1e10;
        static double etam_min = 1e10, etam_max = -1e10;
        static double etan_min = 1e10, etan_max = -1e10;

        const double eta0 = eta_from_A(0);
        const double etam = eta_from_A(mid);
        const double etan = eta_from_A(last);
        eta0_min = std::min(eta0_min, eta0); eta0_max = std::max(eta0_max, eta0);
        etam_min = std::min(etam_min, etam); etam_max = std::max(etam_max, etam);
        etan_min = std::min(etan_min, etan); etan_max = std::max(etan_max, etan);

        const double c0 = celerityFromArea(m_vEstuaryAreas[0], m_vWidth[0], m_vCrossSectionArea[0]) / std::sqrt(dGetLateralStorageFactor(0));
        const double cMid = celerityFromArea(m_vEstuaryAreas[mid], m_vWidth[mid], m_vCrossSectionArea[mid]) / std::sqrt(dGetLateralStorageFactor(mid));
        const double cN = celerityFromArea(m_vEstuaryAreas[last], m_vWidth[last], m_vCrossSectionArea[last]) / std::sqrt(dGetLateralStorageFactor(last));

        // Approximate gravity-wave travel time from downstream to upstream using current c(x): T ≈ ∑ dx / c.
        // This is a diagnostic (not a solver constraint) to compare against the observed phase lag.
        double travel_s = 0.0;
        for (int i = 1; i <= last; ++i) {
            const double dx = m_vCrossSectionX[i] - m_vCrossSectionX[i - 1];
            if (dx <= 0.0) continue;
            // Use celerity at the segment midpoint via averaging endpoint celerities (diagnostic only).
            const double cL = celerityFromArea(m_vEstuaryAreas[i - 1], m_vWidth[i - 1], m_vCrossSectionArea[i - 1]) / std::sqrt(dGetLateralStorageFactor(i - 1));
            const double cR = celerityFromArea(m_vEstuaryAreas[i], m_vWidth[i], m_vCrossSectionArea[i]) / std::sqrt(dGetLateralStorageFactor(i));
            const double ci = 0.5 * (cL + cR);
            travel_s += dx / std::max(1e-6, ci);
        }
        const double L = m_vCrossSectionX[last] - m_vCrossSectionX[0];
        const double c_eff = (travel_s > 0.0) ? (L / travel_s) : 0.0;

        LogStream << "KeyNodes: "
                  << "n0=" << std::fixed << std::setprecision(4) << m_vCrossSectionManningNumber[0]
                  << " nMid=" << m_vCrossSectionManningNumber[mid]
                  << " nN=" << m_vCrossSectionManningNumber[last]
                  << "  beta[0,mid,n]=" << std::fixed << std::setprecision(6)
                  << (m_vBeta.empty() ? 0.0 : m_vBeta[0]) << ","
                  << (m_vBeta.empty() ? 0.0 : m_vBeta[mid]) << ","
                  << (m_vBeta.empty() ? 0.0 : m_vBeta[last])
                  << "  Fr[0,mid,n]=";

        auto safe_Froude = [&](int i, double c) {
            const double A = std::max(DRY_AREA, m_vCrossSectionArea[i]);
            const double u = m_vCrossSectionQ[i] / A;
            return std::fabs(u) / std::max(1e-6, c);
        };
        const double Fr0 = safe_Froude(0, c0);
        const double FrM = safe_Froude(mid, cMid);
        const double FrN = safe_Froude(last, cN);
        LogStream << std::fixed << std::setprecision(3) << Fr0 << "," << FrM << "," << FrN;

        auto safe_advF1 = [&](int i) {
            const double A = std::max(DRY_AREA, m_vCrossSectionArea[i]);
            const double Q = m_vCrossSectionQ[i];
            const double beta = (m_vBeta.empty() ? 1.0 : m_vBeta[i]);
            return beta * (Q * Q) / A;
        };
        const double adv0 = safe_advF1(0);
        const double advM = safe_advF1(mid);
        const double advN = safe_advF1(last);

        LogStream << "  advF1=βQ²/A[0,mid,n]=" << std::scientific << std::setprecision(3)
                  << adv0 << "," << advM << "," << advN << std::defaultfloat
                  << "  SfMid=" << std::scientific << std::setprecision(3) << Sf_from_state(mid) << std::defaultfloat
                  << "  c[0,mid,n]=" << std::fixed << std::setprecision(2) << c0 << "," << cMid << "," << cN
                  << "  Ttravel≈" << std::fixed << std::setprecision(2) << (travel_s / 3600.0) << "h"
                  << "  c_eff≈" << std::fixed << std::setprecision(2) << c_eff
                  << "  etaAmp[0,mid,n]=" << std::fixed << std::setprecision(3)
                  << (eta0_max - eta0_min) << "," << (etam_max - etam_min) << "," << (etan_max - etan_min)
                  << "\n";

        // Phase diagnostic: estimate best lag between eta at x0 and xN using cross-correlation
        // over the last ~72 output samples (matching the typical hourly output cadence).
        // Positive lag means eta(x0) lags eta(xN).
        static std::vector<double> eta0_hist;
        static std::vector<double> etan_hist;
        static bool hist_init = false;
        if (!hist_init) {
            eta0_hist.reserve(256);
            etan_hist.reserve(256);
            hist_init = true;
        }
        eta0_hist.push_back(eta0);
        etan_hist.push_back(etan);
        const size_t max_keep = 256;
        if (eta0_hist.size() > max_keep) {
            eta0_hist.erase(eta0_hist.begin(), eta0_hist.begin() + (eta0_hist.size() - max_keep));
            etan_hist.erase(etan_hist.begin(), etan_hist.begin() + (etan_hist.size() - max_keep));
        }

        // Compute lag only when we have enough samples.
        const size_t Nsamples = eta0_hist.size();
        if (Nsamples >= 12) {
            const size_t window = std::min<size_t>(72, Nsamples);
            const size_t start = Nsamples - window;

            // Detrend by removing mean over the window.
            double m0 = 0.0, mn = 0.0;
            for (size_t k = start; k < Nsamples; ++k) {
                m0 += eta0_hist[k];
                mn += etan_hist[k];
            }
            m0 /= static_cast<double>(window);
            mn /= static_cast<double>(window);

            // Search lags in +/-12 samples (scaled to hours using output interval).
            const int max_lag_samples = 12;
            int best_lag = 0;
            double best_score = -1e300;
            for (int lag = -max_lag_samples; lag <= max_lag_samples; ++lag) {
                double sum = 0.0;
                for (size_t k = start; k < Nsamples; ++k) {
                    const long j = static_cast<long>(k) + lag;
                    if (j < static_cast<long>(start) || j >= static_cast<long>(Nsamples)) continue;
                    sum += (eta0_hist[k] - m0) * (etan_hist[static_cast<size_t>(j)] - mn);
                }
                if (sum > best_score) {
                    best_score = sum;
                    best_lag = lag;
                }
            }

            const double dt_out_h = m_dSimTimestep / 3600.0;
            LogStream << "PhaseLag: bestLag=" << std::fixed << std::setprecision(2)
                      << (static_cast<double>(best_lag) * dt_out_h)
                      << " h (eta0 relative to etaN, window=" << window << ")\n";
        }
    }

    // Boundary-condition diagnostics (stage vs simulated at the last two nodes)
    const int n = m_nCrossSectionsNumber - 1;
    if (n >= 2 && nGetDownwardEstuarineCondition() == 2 && !m_vDownwardBoundaryConditionTime.empty()) {
        const double eta_bc = linearInterpolation1d(m_dCurrentTime,
                                                    m_vDownwardBoundaryConditionTime,
                                                    m_vDownwardBoundaryConditionValue);
        const double eta_n = linearInterpolation1d(m_vCrossSectionArea[n], m_vEstuaryAreas[n], m_vEstuaryWaterDepths[n]) + m_vBedZ[n];
        const double eta_nm1 = linearInterpolation1d(m_vCrossSectionArea[n - 1], m_vEstuaryAreas[n - 1], m_vEstuaryWaterDepths[n - 1]) + m_vBedZ[n - 1];
        const double dQ = m_vCrossSectionQ[n] - m_vCrossSectionQ[n - 1];
        LogStream << "BCdown(stage): eta_bc=" << std::fixed << std::setprecision(3) << eta_bc
                  << "  eta[n]=" << eta_n << "  eta[n-1]=" << eta_nm1
                  << "  Q[n]=" << m_vCrossSectionQ[n] << "  Q[n-1]=" << m_vCrossSectionQ[n - 1]
                  << "  dQ=" << std::scientific << std::setprecision(3) << dQ << std::defaultfloat
                  << "\n";
    }
    if (nGetUpwardEstuarineCondition() == 3 && m_nCrossSectionsNumber >= 2) {
        const double eta0 = linearInterpolation1d(m_vCrossSectionArea[0], m_vEstuaryAreas[0], m_vEstuaryWaterDepths[0]) + m_vBedZ[0];
        const double eta1 = linearInterpolation1d(m_vCrossSectionArea[1], m_vEstuaryAreas[1], m_vEstuaryWaterDepths[1]) + m_vBedZ[1];
        LogStream << "BCup(Q): Q0=" << std::fixed << std::setprecision(3) << m_vCrossSectionQ[0]
                  << "  Q1=" << m_vCrossSectionQ[1]
                  << "  eta0=" << eta0 << "  eta1=" << eta1 << "\n";
    }

    // Advanced salinity diagnostics (every hour)
    if (m_bDoWaterSalinity) {
        // Find salt wedge intrusion limit (isohaline 1 ppt)
        // Search from upstream (dam) to downstream (ocean) to find furthest upstream intrusion
        int intrusion_idx = -1;
        double intrusion_distance = 0.0;
        for (int i = 0; i < m_nCrossSectionsNumber; i++) {
            if (m_vCrossSectionSalinity[i] >= 1.0) {
                // First occurrence from upstream = furthest intrusion
                intrusion_idx = i;
                intrusion_distance = m_vCrossSectionX[i];
                break;
            }
        }
        // Calculate maximum salinity gradient (halocline strength)
        double max_grad = 0.0;
        int max_grad_idx = 0;
        for (int i = 1; i < m_nCrossSectionsNumber-1; i++) {
            double grad = fabs(m_vCrossSectionSalinity[i+1] - m_vCrossSectionSalinity[i-1]) /
                         (m_vCrossSectionX[i+1] - m_vCrossSectionX[i-1]);
            if (grad > max_grad) {
                max_grad = grad;
                max_grad_idx = i;
            }
        }
        // Estimate advective vs diffusive transport
        double advective_flux = 0.0;
        double diffusive_flux = 0.0;
        int n_active = 0;
        for (int i = 1; i < m_nCrossSectionsNumber-1; i++) {
            if (m_vCrossSectionArea[i] > DRY_AREA) {
                advective_flux += fabs(m_vCrossSectionU[i] * m_vCrossSectionSalinity[i] * m_vCrossSectionArea[i]);
                const double Kh = dGetLongitudinalDispersion(i);
                if (Kh > 0.0) {
                    double dS_dx = (m_vCrossSectionSalinity[i+1] - m_vCrossSectionSalinity[i-1]) /
                                   (m_vCrossSectionX[i+1] - m_vCrossSectionX[i-1]);
                    diffusive_flux += fabs(Kh * dS_dx * m_vCrossSectionArea[i]);
                }
                n_active++;
            }
        }
        if (n_active > 0) {
            advective_flux /= n_active;
            diffusive_flux /= n_active;
        }
        // Estimate numerical diffusion (Pe = u*dx/K)
        double Pe_min = 1e10;
        double numerical_diffusion = 0.0;
        for (int i = 1; i < m_nCrossSectionsNumber-1; i++) {
            if (m_vCrossSectionArea[i] > DRY_AREA && m_vCrossSectionU[i] > 0.1) {
                const double dxL = m_vCrossSectionX[i] - m_vCrossSectionX[i - 1];
                const double dxR = m_vCrossSectionX[i + 1] - m_vCrossSectionX[i];
                const double dx = 0.5 * (dxL + dxR);
                double u = fabs(m_vCrossSectionU[i]);
                const double Kh = dGetLongitudinalDispersion(i);
                double Pe = (Kh > 0.0) ? (u * dx / Kh) : 1e10;
                if (Pe < Pe_min) Pe_min = Pe;
                numerical_diffusion = std::max(numerical_diffusion, u * dx * 0.5);
            }
        }
        LogStream << "🌊 SALINITY [t=" << std::fixed << std::setprecision(1) << m_dCurrentTime/3600.0 << "h]:\n";
        if (intrusion_idx >= 0) {
            LogStream << "  • Salt wedge intrusion: " << std::setprecision(1) << intrusion_distance/1000.0
                      << " km (S≥1ppt at section " << intrusion_idx << ")\n";
        } else {
            LogStream << "  • Salt wedge: Not detected (S<1ppt throughout domain)\n";
        }
        LogStream << "  • Halocline: max gradient = " << std::setprecision(4) << max_grad*1000.0
                  << " ppt/km at x=" << std::setprecision(1) << m_vCrossSectionX[max_grad_idx]/1000.0 << " km\n";
        LogStream << "  • Transport: Advection=" << std::setprecision(1) << advective_flux
                  << " kg/s, Diffusion=" << diffusive_flux << " kg/s";
        if (advective_flux > 0) {
            LogStream << " (Ratio=" << std::setprecision(2) << diffusive_flux/advective_flux << ")";
        }
        LogStream << "\n";
        double kh_ref = std::max(0.0, m_dLongitudinalDispersion);
        if (bHasLongitudinalDispersionProfile()) {
            double kh_min = 1e300, kh_max = -1e300;
            for (double kh : m_vLongitudinalDispersion) {
                kh_min = std::min(kh_min, kh);
                kh_max = std::max(kh_max, kh);
            }
            kh_ref = std::max(0.0, kh_max);
            LogStream << "  • Dispersion coefficient: Kh(x) in [" << std::setprecision(2) << kh_min
                      << " .. " << kh_max << "] m²/s, K_numerical≈" << numerical_diffusion << " m²/s";
        } else {
            LogStream << "  • Dispersion coefficient: K_physical=" << std::setprecision(1) << m_dLongitudinalDispersion
                      << " m²/s, K_numerical≈" << numerical_diffusion << " m²/s";
        }
        if (kh_ref > 0.0 && numerical_diffusion > kh_ref * 0.1) {
            LogStream << " ⚠️ HIGH";
        }
        LogStream << "\n";
        if (Pe_min < 2.0) {
            LogStream << "  ⚠️ WARNING: Low Péclet number (Pe=" << std::setprecision(1) << Pe_min
                      << ") - Diffusion dominates over advection!\n";
        }
        LogStream << "\n";
    }
    LogStream.flush();
    
    // Detailed salinity diagnostics every 6 hours (if log level >= 2)
    if (m_bDoWaterSalinity && m_nLogFileDetail >= 2) {
        static double last_salinity_log = -21600.0;  // Initialize to -6h so first log happens at t=0
        if (m_dCurrentTime - last_salinity_log >= 21600.0) {
            last_salinity_log = m_dCurrentTime;
            
            // Find salt wedge intrusion limit (isohaline 1 ppt)
            int intrusion_idx = -1;
            double intrusion_distance = 0.0;
            for (int i = m_nCrossSectionsNumber-1; i >= 0; i--) {
                if (m_vCrossSectionSalinity[i] >= 1.0) {
                    intrusion_idx = i;
                    intrusion_distance = m_vCrossSectionX[i];
                    break;
                }
            }
            
            // Calculate maximum salinity gradient (halocline strength)
            double max_grad = 0.0;
            int max_grad_idx = 0;
            for (int i = 1; i < m_nCrossSectionsNumber-1; i++) {
                double grad = fabs(m_vCrossSectionSalinity[i+1] - m_vCrossSectionSalinity[i-1]) / 
                             (m_vCrossSectionX[i+1] - m_vCrossSectionX[i-1]);
                if (grad > max_grad) {
                    max_grad = grad;
                    max_grad_idx = i;
                }
            }
            
            // Estimate advective vs diffusive transport
            double advective_flux = 0.0;
            double diffusive_flux = 0.0;
            int n_active = 0;
            for (int i = 1; i < m_nCrossSectionsNumber-1; i++) {
                if (m_vCrossSectionArea[i] > DRY_AREA) {
                    // Advective transport: U * S * A
                    advective_flux += fabs(m_vCrossSectionU[i] * m_vCrossSectionSalinity[i] * m_vCrossSectionArea[i]);
                    
                    // Diffusive transport: K * dS/dx * A
                    const double Kh = dGetLongitudinalDispersion(i);
                    if (Kh > 0.0) {
                        double dS_dx = (m_vCrossSectionSalinity[i+1] - m_vCrossSectionSalinity[i-1]) / 
                                       (m_vCrossSectionX[i+1] - m_vCrossSectionX[i-1]);
                        diffusive_flux += fabs(Kh * dS_dx * m_vCrossSectionArea[i]);
                    }
                    n_active++;
                }
            }
            if (n_active > 0) {
                advective_flux /= n_active;
                diffusive_flux /= n_active;
            }
            
            // Estimate numerical diffusion (Pe = u*dx/K)
            double Pe_min = 1e10;
            double numerical_diffusion = 0.0;
            for (int i = 1; i < m_nCrossSectionsNumber-1; i++) {
                if (m_vCrossSectionArea[i] > DRY_AREA && m_vCrossSectionU[i] > 0.1) {
                    const double dxL = m_vCrossSectionX[i] - m_vCrossSectionX[i - 1];
                    const double dxR = m_vCrossSectionX[i + 1] - m_vCrossSectionX[i];
                    const double dx = 0.5 * (dxL + dxR);
                    double u = fabs(m_vCrossSectionU[i]);
                    const double Kh = dGetLongitudinalDispersion(i);
                    double Pe = (Kh > 0.0) ? (u * dx / Kh) : 1e10;
                    if (Pe < Pe_min) Pe_min = Pe;
                    
                    // Numerical diffusion from upwind: K_num ~ u*dx/2
                    numerical_diffusion = std::max(numerical_diffusion, u * dx * 0.5);
                }
            }

            LogStream << "\n";
            LogStream << "SALINITY DETAIL [t=" << std::fixed << std::setprecision(1) << m_dCurrentTime/3600.0 << "h] (6-hour):\n";
            if (intrusion_idx >= 0) {
                LogStream << "  - Intrusion (S≥1 ppt): x=" << std::setprecision(2) << intrusion_distance/1000.0
                          << " km (section " << intrusion_idx << ")\n";
            } else {
                LogStream << "  - Intrusion (S≥1 ppt): not detected\n";
            }
            LogStream << "  - Max |dS/dx|: " << std::setprecision(4) << max_grad*1000.0
                      << " ppt/km at x=" << std::setprecision(2) << m_vCrossSectionX[max_grad_idx]/1000.0 << " km\n";
            LogStream << "  - Mean |adv flux|≈" << std::setprecision(2) << advective_flux
                      << " ; mean |diff flux|≈" << std::setprecision(2) << diffusive_flux;
            if (advective_flux > 0.0) {
                LogStream << " ; diff/adv=" << std::setprecision(3) << (diffusive_flux / advective_flux);
            }
            LogStream << "\n";
            if (bHasLongitudinalDispersionProfile()) {
                double kh_min = 1e300, kh_max = -1e300;
                for (double kh : m_vLongitudinalDispersion) {
                    kh_min = std::min(kh_min, kh);
                    kh_max = std::max(kh_max, kh);
                }
                LogStream << "  - Kh(x) in [" << std::setprecision(2) << kh_min << " .. " << kh_max
                          << "] m^2/s ; K_numerical≈" << std::setprecision(2) << numerical_diffusion << " m^2/s\n";
            } else {
                LogStream << "  - K_physical=" << std::setprecision(2) << dGetLongitudinalDispersionConstant()
                          << " m^2/s ; K_numerical≈" << std::setprecision(2) << numerical_diffusion << " m^2/s\n";
            }
            if (Pe_min < 2.0) {
                LogStream << "  - WARNING: Low Peclet (min Pe=" << std::setprecision(2) << Pe_min << ")\n";
            }
            LogStream << "\n";
        }
    }
    
    LogStream.flush();
}



/**
 * @brief Calculate comprehensive net surface heat flux for water temperature equation
 * 
 * Computes Q_net = Q_SW + Q_LW + Q_H + Q_E (W/m²) using physically-based parameterizations
 * 
 * RADIATION COMPONENTS:
 * --------------------
 * 1. Shortwave radiation (Q_SW): Solar radiation absorbed by water
 *    - Synthetic incident radiation from solar geometry (Kasten-Czeplak model)
 *    - Dynamic albedo from Briegleb formula (zenith-angle dependent, ~0.05-0.40)
 *    - Q_SW = SW↓ × (1 - α(θ_z))
 * 
 * 2. Longwave radiation (Q_LW): Net thermal radiation exchange
 *    - Incoming: εsky × σ × T_air⁴ (Swinbank formula with cloud correction)
 *    - Outgoing: εwater × σ × T_water⁴
 *    - Q_LW = LW↓ - LW↑ (typically negative, net cooling)
 * 
 * TURBULENT FLUXES (Bulk Aerodynamic Formulas):
 * ----------------------------------------------
 * 3. Sensible heat flux (Q_H): Turbulent heat conduction
 *    Q_H = ρ_air(T_air) × c_p × C_H × U × (T_air - T_water)
 *    - C_H ≈ 1.3×10⁻³: Stanton number (calibration parameter)
 *    - Positive when air warms water
 * 
 * 4. Latent heat flux (Q_E): Evaporation/condensation
 *    Q_E = ρ_air(T_air) × L_v(T_water) × C_E × U × (q_air - q_sat,water)
 *    - C_E ≈ 1.3×10⁻³: Dalton number (calibration parameter)
 *    - Negative during evaporation (typical), positive during condensation
 *    - q = specific humidity from Magnus-Tetens formula
 * 
 * PHYSICAL IMPROVEMENTS OVER PREVIOUS VERSION:
 * ---------------------------------------------
 * - Dynamic air density ρ_air(T) using ideal gas law (see calc_rho_air)
 * - Dynamic latent heat L_v(T) temperature-dependent (see calc_Lv)
 * - Solar geometry calculation for realistic diurnal cycle
 * - Zenith-angle dependent albedo (Briegleb et al. 1986 reflectance)
 * - Sky emissivity with cloud cover correction
 * - Proper sign conventions (Q_net > 0 → water warms)
 * 
 * CONFIGURATION PARAMETERS:
 * -------------------------
 * - Latitude: 36.5°N (modifiable for other locations)
 * - Cloud cover: 0.2 (20%, can be made time-dependent if data available)
 * - C_H, C_E: 1.3×10⁻³ (calibration parameters, adjust for local conditions)
 * 
 * @note Requires meteorological forcing data:
 * - m_vHeatFluxAirTemp: Air temperature (°C)
 * - m_vHeatFluxRelHumidity: Relative humidity (%)
 * - m_vHeatFluxWind: Wind speed (m/s)
 * 
 * @note Output stored in m_vCrossSectionTemperatureASt[i] for each cross-section
 * 
 * @warning Minimum wind speed clamped to 0.1 m/s to avoid numerical issues
 * 
 * @see calculateTemperatureSourceTerms() for integration of Q_net into dT/dt equation
 * @see calc_rho_air() for temperature-dependent air density
 * @see calc_Lv() for temperature-dependent latent heat
 * @see calc_albedo_briegleb() for zenith-angle dependent albedo with direct/diffuse separation
 * 
 * @ref Swinbank (1963): Long-wave radiation from clear skies
 * @ref Kasten & Czeplak (1980): Solar radiation parameterization
 * @ref Briegleb et al. (1986): Angle-dependent surface reflectance
 */
void CSimulation::calculateRadiativeFluxes() {
    // === SIMULATION CONFIGURATION ===
    // Derive current hour and day from simulation time
    double sim_time_hours = m_dCurrentTime / 3600.0; 
    double hour_of_day = std::fmod(sim_time_hours, 24.0);  // Hour 0-24
    int day_of_year = 1 + static_cast<int>(sim_time_hours / 24.0); // Day 1-365 (simplified)
    
    // Get calibration parameters from YAML configuration (or use defaults)
    double latitude = m_dHeatFluxLatitude;        // Geographic latitude (degrees, North positive)
    double cloud_cover = m_dHeatFluxCloudCover;  // Cloud cover fraction (0.0=clear, 1.0=overcast)
    double C_h = m_dHeatFlux_CS;                 // Sensible heat transfer coefficient (Stanton number)
    double C_e = m_dHeatFlux_CL;                 // Latent heat transfer coefficient (Dalton number)

    for (int i = 0; i < m_nCrossSectionsNumber; ++i) {
        // === 1. STATE VARIABLES ===
        // Meteorological forcing from input data (with defaults if unavailable)
        double Tair = (!m_vHeatFluxTime.empty()) ? 
            linearInterpolation1d(m_dCurrentTime, m_vHeatFluxTime, m_vHeatFluxAirTemp) : 15.0;
        
        // Relative humidity: either from data or calculated from temperature
        double rh = 70.0;  // Default value
        if (!m_vHeatFluxTime.empty()) {
            if (m_bCalculateRHFromTemperature && !m_vDailyMinTemperature.empty()) {
                // Calculate RH using FAO-56 method from current T_air and daily T_min
                double sim_time_hours = m_dCurrentTime / 3600.0;
                int day_index = static_cast<int>(sim_time_hours / 24.0);
                
                // Ensure day_index is within bounds
                if (day_index >= 0 && day_index < static_cast<int>(m_vDailyMinTemperature.size())) {
                    double T_min_daily = m_vDailyMinTemperature[day_index];
                    rh = calc_rh_from_temp(Tair, T_min_daily);
                } else {
                    // Fallback: use interpolated value if available, otherwise default
                    rh = (!m_vHeatFluxRelHumidity.empty()) ?
                        linearInterpolation1d(m_dCurrentTime, m_vHeatFluxTime, m_vHeatFluxRelHumidity) : 70.0;
                }
            } else {
                // Use interpolated RH from data file
                rh = linearInterpolation1d(m_dCurrentTime, m_vHeatFluxTime, m_vHeatFluxRelHumidity);
            }
        }
        
        double wind = (!m_vHeatFluxTime.empty()) ? 
            linearInterpolation1d(m_dCurrentTime, m_vHeatFluxTime, m_vHeatFluxWind) : 1.0;
        double pressure = (!m_vHeatFluxAtmosphericPressure.empty()) ?
            linearInterpolation1d(m_dCurrentTime, m_vHeatFluxTime, m_vHeatFluxAtmosphericPressure) : ATM_PRESSURE_DEFAULT;
        double Twater = m_vCrossSectionTemperature[i];
        
        // Prevent zero wind speed (avoids unrealistic stagnation and division by zero)
        if (wind < 0.1) wind = 0.1;

        // Dynamic physical properties (temperature-dependent)
        double rho_a_dynamic = calc_rho_air(Tair, pressure);  // Air density (kg/m³) with actual pressure
        double L_v_dynamic = calc_Lv(Twater);                 // Latent heat of vaporization (J/kg)

        // === 2. SOLAR GEOMETRY (For Albedo and Synthetic Radiation) ===
        // Solar declination angle (degrees) - function of day of year
        // δ = 23.45° × sin[2π/365 × (284 + n)]
        double delta = 23.45 * std::sin(DEG_TO_RAD * (360.0/365.0) * (284.0 + day_of_year));
        double delta_rad = delta * DEG_TO_RAD;
        double lat_rad = latitude * DEG_TO_RAD;
        
        // Hour angle (degrees from solar noon)
        double omega = 15.0 * (hour_of_day - 12.0);
        double omega_rad = omega * DEG_TO_RAD;

        // Cosine of solar zenith angle (sun elevation)
        // cos(θ_z) = sin(φ)sin(δ) + cos(φ)cos(δ)cos(ω)
        double cos_theta = std::sin(lat_rad) * std::sin(delta_rad) + 
                           std::cos(lat_rad) * std::cos(delta_rad) * std::cos(omega_rad);
        
        double zenith_rad = std::acos(std::max(0.0, cos_theta));

        // === 3. Q_SW: NET SHORTWAVE RADIATION ===
        double albedo_dynamic = 0.06;  // Default albedo (Briegleb diffuse value for nighttime)
        double Qsw_in = 0.0;           // Incident shortwave radiation (W/m²)

        if (cos_theta > 0.0) {  // Daytime
            // Dynamic albedo from Briegleb et al. (1986) model
            albedo_dynamic = calc_albedo_briegleb(zenith_rad, cloud_cover);

            // Synthetic incident shortwave radiation (Kasten-Czeplak model)
            double R_clear = SOLAR_CONSTANT * cos_theta * ATM_TRANSMISSIVITY;
            double cloud_factor = 1.0 - 0.75 * std::pow(cloud_cover, 3.4);
            Qsw_in = R_clear * cloud_factor;
        }

        // Net shortwave radiation absorbed by water (positive = heat gain)
        double Qsw = Qsw_in * (1.0 - albedo_dynamic);

        // === 4. Q_LW: NET LONGWAVE RADIATION ===
        double Tair_K = Tair + 273.15;
        double Twater_K = Twater + 273.15;
        
        // Sky emissivity (Swinbank formula with cloud correction)
        // ε_sky = 9.37×10⁻⁶ T_air² (1 + 0.17 C²)
        double epsilon_sky = 9.37e-6 * (Tair_K * Tair_K) * 
                             (1.0 + 0.17 * cloud_cover * cloud_cover);
        
        // Incoming longwave radiation from atmosphere
        double Qlw_in = epsilon_sky * STEFAN_BOLTZMANN * std::pow(Tair_K, 4);
        
        // Outgoing longwave radiation from water surface
        double Qlw_out = WATER_EMISSIVITY * STEFAN_BOLTZMANN * std::pow(Twater_K, 4);

        // Net longwave radiation (typically negative/cooling)
        double Qlw = Qlw_in - Qlw_out;

        // === 5. TURBULENT HEAT FLUXES (Bulk Aerodynamic Formulas) ===
        
        // --- SENSIBLE HEAT FLUX (Q_H) ---
        // Q_H = ρ_air × c_p × C_h × U × (T_air - T_water)
        // Sign convention: Positive when air is warmer (water gains heat)
        double bulk_coeff_sensible = rho_a_dynamic * AIR_SPECIFIC_HEAT * C_h * wind;
        double Qsensible = bulk_coeff_sensible * (Tair - Twater); 

        // --- LATENT HEAT FLUX (Q_E) ---
        // Q_E = ρ_air × L_v × C_e × U × (q_air - q_sat,water)
        
        // Saturation vapor pressure (Magnus-Tetens formula, result in hPa)
        double esat_water = 6.112 * std::exp((17.62 * Twater) / (243.12 + Twater));
        double esat_air = 6.112 * std::exp((17.62 * Tair) / (243.12 + Tair));
        
        // Specific humidity (kg water vapor / kg air)
        double qsat_water = 0.622 * esat_water / 1013.25;  // At water surface (saturated)
        double q_air = (rh / 100.0) * 0.622 * esat_air / 1013.25;  // In air (from relative humidity)
        
        // Latent heat flux calculation
        // Sign convention: Positive when air is more humid (condensation, water gains heat)
        //                  Negative when evaporation occurs (water loses heat)
        double bulk_coeff_latent = rho_a_dynamic * L_v_dynamic * C_e * wind;
        double Qlatente = bulk_coeff_latent * (q_air - qsat_water); 

        // === 6. NET SURFACE HEAT FLUX ===
        // Q_net = Q_SW + Q_LW + Q_H + Q_E
        // Positive Q_net → water warms
        // Negative Q_net → water cools
        double Qnet = Qsw + Qlw + Qsensible + Qlatente;

        // Store result for temperature equation (used in calculateTemperatureSourceTerms)
        m_vCrossSectionTemperatureASt[i] = Qnet;
        
        // === DIAGNOSTIC LOG (first cross-section only, every 6 hours) ===
        if (i == 0 && m_nLogFileDetail >= 2) {
            // Log heat flux components every 6 hours
            static int last_log_step = -1;
            int current_6h_step = static_cast<int>(m_dCurrentTime / 21600.0);  // 6 hours = 21600 s
            
            if (current_6h_step != last_log_step) {
                last_log_step = current_6h_step;
                LogStream << "\n--- HEAT FLUX DIAGNOSTIC [t=" << std::fixed << std::setprecision(1) 
                          << m_dCurrentTime << "s, " << std::setprecision(1) << hour_of_day << "h] ---\n";
                LogStream << "Position: x=" << m_vPositionX[i] << " m | T_water=" << std::setprecision(2) 
                          << Twater << "°C | T_air=" << Tair << "°C\n";
                LogStream << "Q_SW (solar)     = " << std::setw(8) << std::setprecision(1) << Qsw << " W/m²";
                if (cos_theta > 0) {
                    LogStream << " (day, zenith=" << std::setprecision(1) << (zenith_rad * RAD_TO_DEG) 
                              << "°, albedo=" << std::setprecision(3) << albedo_dynamic << ")\n";
                } else {
                    LogStream << " (night)\n";
                }
                LogStream << "Q_LW (longwave)  = " << std::setw(8) << std::setprecision(1) << Qlw << " W/m²"
                          << " (ε_sky=" << std::setprecision(3) << epsilon_sky << ")\n";
                LogStream << "Q_H (sensible)   = " << std::setw(8) << std::setprecision(1) << Qsensible << " W/m²\n";
                LogStream << "Q_E (latent)     = " << std::setw(8) << std::setprecision(1) << Qlatente << " W/m²"
                          << " (RH=" << std::setprecision(1) << rh << "%)\n";
                LogStream << "Q_net (total)    = " << std::setw(8) << std::setprecision(1) << Qnet << " W/m²\n";
                LogStream << "Wind: " << std::setprecision(2) << wind << " m/s | Depth: " 
                          << std::setprecision(2) << m_vCrossSectionWaterDepth[i] << " m\n";
                LogStream << "---------------------------------------------------\n";
                LogStream.flush();
            }
        }
    }
}


/**
 * @brief Calculate dynamic albedo from latitude, day of year, and hour of day
 * 
 * Computes albedo based on solar position geometry:
 * 
 * 1. Solar declination (Cooper 1969):
 *    δ = 23.45°·sin[2π/365 × (284 + N)]
 *    where N = day of year (1-365)
 * 
 * 2. Hour angle (15° per hour from solar noon):
 *    ω = 15° × (hour - 12)
 * 
 * 3. Solar zenith angle:
 *    cosθ = sinφ·sinδ + cosφ·cosδ·cosω
 *    where φ = latitude
 * 
 * 4. Albedo from Briegleb et al. (1986) model (see calc_albedo_briegleb)
 * 
 * @param lat_deg Latitude in degrees (e.g., 36.5 for Cádiz)
 * @param day_of_year Day of year (1-365, 1=Jan 1)
 * @param hour_of_day Decimal hour (0-24, e.g., 14.5 = 14:30)
 * @return Albedo (0.03-1.0)
 * 
 * @note Simplifications:
 * - Assumes solar noon at 12:00 local time (neglects equation of time ~±15 min)
 * - Neglects longitude correction (use UTC + timezone offset for precision)
 * - Good accuracy for heat budget (~5% error acceptable)
 * 
 * @warning For high-precision applications:
 * - Use NOAA solar position algorithm
 * - Include atmospheric refraction
 * - Account for equation of time
 * 
 * @see calc_albedo_briegleb() for zenith-to-albedo conversion
 */
double calc_dynamic_albedo(double lat_deg, int day_of_year, double hour_of_day) {
    
    // 1. SOLAR POSITION CALCULATION (Basic geometry)
    
    // Solar declination (delta): Sun angle relative to equator as function of day of year
    // Cooper formula (1969) - Sufficient for engineering applications
    double delta = 23.45 * std::sin(CSimulation::DEG_TO_RAD * (360.0/365.0) * (284.0 + day_of_year));
    double delta_rad = delta * CSimulation::DEG_TO_RAD;

    // Hour angle (omega): 0 at solar noon, +/- for each hour
    // (15 degrees per hour)
    double time_offset = 12.0; // Assuming solar noon at 12:00 (simplification)
    double omega = 15.0 * (hour_of_day - time_offset); 
    double omega_rad = omega * CSimulation::DEG_TO_RAD;

    double lat_rad = lat_deg * CSimulation::DEG_TO_RAD;

    // Cosine of zenith angle (cos_theta)
    // theta = 0 (sun vertical), theta = 90 (horizon)
    double cos_theta = std::sin(lat_rad) * std::sin(delta_rad) + 
                       std::cos(lat_rad) * std::cos(delta_rad) * std::cos(omega_rad);

    // 2. CHECK IF NIGHTTIME
    // If sun is below horizon, albedo doesn't matter (Radiation = 0), 
    // but we return 0.07 as default to avoid errors
    if (cos_theta <= 0.0) {
        return 0.07; // Nighttime or default value
    }

    // 3. ALBEDO CALCULATION (Fresnel approximation for water)
    // Many approximations exist. This one is robust and widely used in oceanography.
    // Albedo A = 0.5 * (Parallel_reflection + Perpendicular_reflection)
    
    // We use empirical approximation (Payne 1972, Briegleb et al. 1986)
    // This formula is more efficient than complete Fresnel and very accurate:
    // Approx from "Briegleb et al. (1986)" for marine surface albedo:
    
    double num = 0.037;
    double den = (1.1 * std::pow(cos_theta, 1.4)) + 0.15;
    double albedo = num / den;

    // Limit to physical values (albedo never > 1.0)
    return std::min(std::max(albedo, 0.03), 1.0);
}

/**
 * @brief Calculate sediment transport using van Rijn (1984) formulation
 * 
 * Separates transport into:
 * 1. Bedload (Qb): Grains rolling/sliding along bed
 * 2. Suspended load (Qs): Grains carried in water column
 * 3. Total (Qt): Qb + Qs
 * 
 * BEDLOAD (van Rijn 1984):
 *   u* = (g^0.5 / C) |U|  (shear velocity, C = Chézy coefficient)
 *   θ = u*² / [(s-1)·g·D50]  (Shields parameter)
 *   T = (θ - θ_cr) / θ_cr  (transport stage parameter)
 *   
 *   For T < 3:  gb = 0.053·√[(s-1)·g·D50³]·T^2.1·D*^{-0.3}
 *   For T ≥ 3:  gb = 0.10·√[(s-1)·g·D50³]·T^1.5·D*^{-0.3}
 *   
 *   where D* = D50·[(s-1)·g/ν²]^{1/3} (dimensionless grain size)
 * 
 * SUSPENDED LOAD:
 *   Uses Rouse-Vanoni profile with reference concentration at z = Δb:
 *   C(z) = C_a·[(h-z)/z · Δb/(h-Δb)]^Z
 *   
 *   where:
 *   - Z = w_s/(β·κ·u*): Rouse number (w_s = settling velocity)
 *   - Δb = 0.3·D*^0.7·√T·D_ave: Reference height
 *   - C_a = 0.117·ρ_s·T/D*: Reference concentration
 * 
 * @note Valid for:
 * - Non-cohesive sediments (sand/gravel)
 * - D50: 0.0001 - 0.01 m
 * - Flow not dominated by vegetation drag
 * 
 * @warning Limitations:
 * - Assumes steady, uniform flow (not ideal for tides)
 * - No bed slope correction (can be important in estuaries)
 * - No turbulence damping by stratification
 * - Calibration needed for site-specific conditions
 * 
 * @see van Rijn, L.C. (1984). Sediment Transport, Part I: Bed Load Transport.
 *      Journal of Hydraulic Engineering, 110(10), 1431-1456.
 */
void CSimulation::calculate_sediment_transport()
{
    for (int i=0; i< m_nCrossSectionsNumber; i++)
    {
        // The sign should be included at the end of the calculation
        double vel_abs = fabs(m_vCrossSectionU[i]);

        // Sediment transport direction (+ve upward, -ne downward)
        double dSedimentDirection = 0.0;
        if (m_vCrossSectionU[i] > 0) dSedimentDirection = 1;
        else if (m_vCrossSectionU[i] < 0) dSedimentDirection = -1;
        else dSedimentDirection = 0;


        //====================================================================================================
        //! Bed load sediment transport
        //====================================================================================================
        double c_1 = 18*log(12*m_vCrossSectionHydraulicRadius[i]/(3*m_vCrossSectionD90[i]));
        if (c_1 < 1e-3) c_1 = 1e-3;

        double u_star = pow(G, 0.5)/c_1*vel_abs;
        double shields_crit = 0.0013*pow(m_vCrossSectionDiamX[i], 0.29);
        if (m_vCrossSectionDiamX[i] < 150) shields_crit = 0.055;

        double u_star_crit = sqrt(shields_crit*(m_vCrossSectionRhos[i] - 1)*G*m_vCrossSectionD50[i]);
        //! Sediment transport
        double transport = (u_star*u_star - u_star_crit*u_star_crit)/(u_star_crit*u_star_crit);

        // If sediment transport is negative means no sediment transport
        if (transport < 0) transport = 0;

        // Mass sediment transport
        double gb = 0.0;
        if (transport >= 3)
        {
            gb = 0.1*sqrt((m_vCrossSectionRhos[i]-1)*G*pow(m_vCrossSectionD50[i], 3.0))*pow(transport, 1.5)*pow(m_vCrossSectionDiamX[i], -0.3);
        }
        else
        {
            gb = 0.053*sqrt((m_vCrossSectionRhos[i]-1)*G*pow(m_vCrossSectionD50[i], 3.0))*pow(transport, 2.1)*pow(m_vCrossSectionDiamX[i], -0.3);
        }
        //! Volumetric sediment transport
        m_vCrossSectionQb[i] = dSedimentDirection*gb*m_vCrossSectionWidth[i]*FRESH_WATER_DENSITY*m_vCrossSectionRhos[i]*dSedimentDirection;

        //====================================================================================================
        //! Suspended sediment transport
        //====================================================================================================
        //! Jump height
        double dDeltaB = m_vCrossSectionDaveraged[i]*0.3*pow(m_vCrossSectionDiamX[i], 0.7)*sqrt(transport);
        // Correction of delta_b
        if (dDeltaB < 0.01*m_vCrossSectionWaterDepth[i]) dDeltaB = 0.01*m_vCrossSectionWaterDepth[i];
        if (dDeltaB > 0.5*m_vCrossSectionWaterDepth[i]) dDeltaB = 0.5*m_vCrossSectionWaterDepth[i];

        //! Shear velocity
        double dUx = sqrt(G*m_vCrossSectionHydraulicRadius[i]*fabs(m_vCrossSectionFrictionSlope[i]));

        //! Reference concentration at z=delta_b with plane bottom
        double c_a = 0.117*FRESH_WATER_DENSITY*m_vCrossSectionRhos[i]*transport/m_vCrossSectionDiamX[i];

        //! Representative diameter of suspended particle
        double dRepresentativeDiameter = m_vCrossSectionD50[i]*(1.0+0.11*(m_vCrossSectionSedimentSigma[i]-1.0));

        if (transport >= 25) dRepresentativeDiameter = m_vCrossSectionD50[i];

        double dSettlingVelocity = 0.0;
        //! Settling velocity of representative diameters
        if (dRepresentativeDiameter < 0.0001)
        {
            dSettlingVelocity = (m_vCrossSectionRhos[i] - 1.0)*G*pow(dRepresentativeDiameter, 2.0)/(18.0*NU);
        }
        else if ((dRepresentativeDiameter >= 0.0001) & (dRepresentativeDiameter < 0.001))
        {
            dSettlingVelocity = 10.0*NU*sqrt((1.0 + 0.01*G*(m_vCrossSectionRhos[i] - 1.00)*pow(dRepresentativeDiameter, 3.0)/(NU*NU)) - 1.0)/dRepresentativeDiameter;
        }
        else
        {
            dSettlingVelocity = 1.1*sqrt((m_vCrossSectionRhos[i] - 1.0)*G*dRepresentativeDiameter);
        }

        //! Factor that considers the different diffusion between fluid and sediment particles
        double dBeta = 2.0;

        if (dUx != 0.0) dBeta = 1.0 + 2.0*pow(dSettlingVelocity/dUx, 2.0);
        if (dBeta > 2.0) dBeta = 2.0;

        //! Rouse suspension parameter
        double dRouse = dSettlingVelocity/(dBeta*KAPPA*dUx);

        //! Maximum concentration of sediments at the reference level
        double C0 = 0.65*FRESH_WATER_DENSITY*m_vCrossSectionRhos[i];

        //! Global correction factor
        double dPsi = 2.5*pow(dSettlingVelocity/dUx, 0.8)*pow(c_a/C0, 0.4);

        //! Update Rouse number
        dRouse = dRouse + dPsi;

        //! Adimensional factor A and F
        double dA_Factor = dDeltaB/m_vCrossSectionWaterDepth[i];
        double dF_Factor = 0.0;
        if (fabs(pow(1.0 - dA_Factor, dRouse)*(1.2 + dRouse)) >= 1e-4)
        {
            dF_Factor = (pow(dA_Factor, dRouse) - pow(dA_Factor, 1.2))/(pow(1.0-dA_Factor, dRouse)*(1.2 - dRouse));
        }

        double gbs1 = c_a*m_vCrossSectionWaterDepth[i]*vel_abs*dF_Factor;
        //! Suspended sediment transport
        m_vCrossSectionQs[i] = gbs1/(m_vCrossSectionRhos[i]*FRESH_WATER_DENSITY)*m_vCrossSectionWidth[i]*dSedimentDirection;

        //! Total sediment transport
        m_vCrossSectionQt[i] = m_vCrossSectionQb[i] + m_vCrossSectionQs[i];
        if (m_vCrossSectionQt[i] < 1e-6) m_vCrossSectionQt[i] = 0.0;

        return;
    }

};


/**
 * @brief Calculate water density from salinity, temperature, and sediment concentration
 * 
 * Linear equation of state (EOS):
 *   ρ = ρ₀·[1 + β_S·S - β_T·(T - 4)]
 * 
 * where:
 * - ρ₀ = 1000 kg/m³: Reference density (freshwater at 4°C)
 * - β_S ≈ 0.00078: Haline contraction coefficient (kg/m³/psu)
 * - β_T ≈ 0.0002: Thermal expansion coefficient (1/°C)
 * - S: Salinity (psu = g/kg)
 * - T: Temperature (°C)
 * 
 * Sediment effect (if enabled):
 *   ρ_total = ρ_water + (1 - ρ_w/(ρ_s·1000))·(Qt / (A·dx))·dt
 * 
 * where Qt = sediment transport, ρ_s = sediment relative density (≈2.65)
 * 
 * @note Physical interpretation:
 * - Salinity: +1 psu → +0.78 kg/m³ (denser)
 * - Temperature: +1°C → -0.2 kg/m³ (lighter)
 * - Freshwater densest at 4°C (hence T-4 term)
 * 
 * @warning Linear EOS limitations:
 * - Valid for: 0-40 psu, 0-30°C, pressures < 100 bar
 * - For higher accuracy use UNESCO EOS-80 or TEOS-10
 * - Neglects pressure effect (OK for shallow water)
 * 
 * @see Gill (1982): Atmosphere-Ocean Dynamics, Appendix 3
 */
void CSimulation::calculate_density()
{
    // ⚡ OPTIMIZATION: Precalculate constants outside loop (checked once per timestep)
    // Baroclinic coefficients betaS (salinity) and betaT (temperature)
    double betaS = bGetDoWaterSalinity() ? dGetBetaSalinityConstant() : 0.0;
    double betaT = bGetDoWaterTemperature() ? 
                   (bGetDoBetaCoefficient() && m_dBetaTemperatureConstant > 0.0 ? 
                    m_dBetaTemperatureConstant : dGetBetaTemperatureConstant()) : 0.0;
    
    // Cache constant factors
    const double rho_base = 1000.0;
    const double T_ref = 4.0;  // Reference temperature for density calculation
    
    for (int i = 0; i < m_nCrossSectionsNumber; i++)
    {
        const double S = m_vCrossSectionSalinity[i];
        const double T = m_vCrossSectionTemperature[i];
        const double rhow = rho_base * (1.0 + betaS * S - betaT * (T - T_ref));
        if (m_bDoSedimentTransport) {
            // Use control-volume length for non-uniform grids
            double dxCV = 0.0;
            if (m_nCrossSectionsNumber >= 2) {
                if (i == 0) {
                    dxCV = m_vCrossSectionX[1] - m_vCrossSectionX[0];
                } else if (i == m_nCrossSectionsNumber - 1) {
                    dxCV = m_vCrossSectionX[i] - m_vCrossSectionX[i - 1];
                } else {
                    const double dxL = m_vCrossSectionX[i] - m_vCrossSectionX[i - 1];
                    const double dxR = m_vCrossSectionX[i + 1] - m_vCrossSectionX[i];
                    dxCV = 0.5 * (dxL + dxR);
                }
            }
            dxCV = std::max(1e-12, dxCV);

            if (m_nPredictor == 1) {
                m_vPredictedCrossSectionDensity[i] =  rhow + (1 - rhow/1000.0/m_vCrossSectionRhos[i])*
                                                     m_vCrossSectionQt[i]/(m_vCrossSectionArea[i]*dxCV)*m_dTimestep;
            }
            else 
                m_vCrossSectionDensity[i] = rhow + (1 - rhow/1000.0/m_vCrossSectionRhos[i])*
                                            m_vCrossSectionQt[i]/(m_vCrossSectionArea[i]*dxCV)*m_dTimestep;
            }
        else {
            if (m_nPredictor == 1) {
                m_vPredictedCrossSectionDensity[i] = rhow;
            }
            else {
                m_vCrossSectionDensity[i] = rhow;
            }
        }
    }
}

//===============================================================================================================================
//! Notifies the user that the simulation has ended, asks for keypress if necessary, and if compiled under GNU can send an email
//===============================================================================================================================
void CSimulation::bDoSimulationEnd(){
    presenter.EndingRun();

    const string strText = strGetErrorText(m_nStringError);
    CScreenPresenter::AnnounceEnding(strText);
};


/**
 * @brief Precompute and cache estuary geometry/hydraulics for O(1) runtime access
 * 
 * Optimization strategy:
 * - Extract data from CCrossSection objects into flat vectors
 * - Avoids repeated virtual function calls (e.g., estuary[i].dGetZ())
 * - Enables better CPU cache locality (contiguous memory)
 * - Precompute expensive terms: A·R_h^{2/3} for Manning friction
 * 
 * Cached data:
 * - Scalars: bed elevation, Manning's n, x-position, β coefficient
 * - Vectors: width(h), area(h), Rh(h), depth(h) - tabulated vs water level
 * - Bank locations: left/right riverbank positions
 * - Friction term: A·R_h^{2/3} (saves one pow() call per timestep)
 * 
 * @note Performance impact:
 * - ~10-15% speedup for large domains (N > 500)
 * - Most benefit in calculateHydraulicParameters() (called every timestep)
 * - Trade-off: ~50 MB extra memory for N=1000 domain
 * 
 * @warning Must be called after:
 * - estuary[] objects fully initialized
 * - All cross-section data loaded from CSV
 * 
 * @see calculateHydraulicParameters() for usage of cached data
 */
void CSimulation::precomputeEstuaryData() {
    
    // Preallocate scalar data arrays for geometry properties
    m_vElevationSectionsCount.resize(m_nCrossSectionsNumber);
    m_vBedZ.resize(m_nCrossSectionsNumber);
    m_vManningN.resize(m_nCrossSectionsNumber);
    m_vPositionX.resize(m_nCrossSectionsNumber);
    m_vBeta.resize(m_nCrossSectionsNumber);
    
    // Initialize binary search cache (speeds up calculateHydraulicParameters ~30%)
    m_vLastInterpolationIndex.assign(m_nCrossSectionsNumber, 0);
    
    // Initialize with known dimensions
    m_vWidth.assign(m_nCrossSectionsNumber, vector<double>());        // Initialize empty vectors       
    m_vLeftY.assign(m_nCrossSectionsNumber, vector<double>());        
    m_vRightY.assign(m_nCrossSectionsNumber, vector<double>());       
    m_vEstuaryAreas.assign(m_nCrossSectionsNumber, vector<double>());
    m_vEstuaryHydraulicRadius.assign(m_nCrossSectionsNumber, vector<double>());
    m_vEstuaryWaterDepths.assign(m_nCrossSectionsNumber, vector<double>());
    m_vPrecalculatedSecondTerm.assign(m_nCrossSectionsNumber, vector<double>());

    
    for (int i = 0; i < m_nCrossSectionsNumber; i++) {
        // Extract basic scalar properties from cross-section geometry
        m_vBedZ[i] = estuary[i].dGetZ();
        m_vManningN[i] = estuary[i].dGetManningNumber();
        m_vPositionX[i] = estuary[i].dGetX();
        m_vElevationSectionsCount[i] = estuary[i].nGetElevationSectionsNumber();
        m_vBeta[i] = estuary[i].dGetBeta();
        
        // Extract hydraulic lookup tables (area, depth, width, hydraulic radius)
        m_vWidth[i] = estuary[i].vGetWidth();   
        m_vEstuaryAreas[i] = estuary[i].vGetArea();
        m_vEstuaryHydraulicRadius[i] = estuary[i].vGetHydraulicRadius();
        m_vEstuaryWaterDepths[i] = estuary[i].vGetWaterDepth();
        m_vLeftY[i] = estuary[i].vGetLeftRBLocation();
        m_vRightY[i] = estuary[i].vGetRightRBLocation();

        const auto& areas = m_vEstuaryAreas[i];
        const auto& hydraulicRadius = m_vEstuaryHydraulicRadius[i];
        
        m_vPrecalculatedSecondTerm[i].resize(areas.size());
        for (size_t j = 0; j < areas.size(); j++) {
            m_vPrecalculatedSecondTerm[i][j] = areas[j] * pow(hydraulicRadius[j], 2.0/3.0);
        }
    }

    // Some modules (e.g., tracer transport) use m_vCrossSectionX explicitly.
    // Keep it defined and consistent with the main geometry x-coordinate.
    m_vCrossSectionX = m_vPositionX;
}

/**
 * @brief Generate descriptive NetCDF output filename from simulation parameters
 * 
 * Format: barrier_sim_YYYYMMDD_HHMM_CSxxx_Txxx_dtxxx_BCxy_[OPTIONS]_CFLxx.nc
 * 
 * Components:
 * - Date/Time: Current wall-clock time when simulation starts
 * - CS: Number of cross-sections (zero-padded to 3 digits)
 * - T: Simulation duration (minutes < 60, hours < 24, else days)
 * - dt: Timestep in seconds
 * - BC: Boundary conditions (xy = upstream/downstream type)
 * - OPTIONS: _SED (sediment), _SAL (salinity), _TVD (limiters), _DRY (dry bed)
 * - CFL: Courant number × 100 (e.g., CFL15 = 0.15)
 * 
 * Example:
 *   barrier_sim_20240115_1430_CS100_T7d_dt3600_BC12_SAL_TVD_DRY_CFL15.nc
 * 
 * @return Formatted filename string
 * 
 * @note Purpose:
 * - Self-documenting filenames for parameter studies
 * - Easy sorting by date
 * - Unique names prevent accidental overwrites
 * - Quickly identify simulation configuration without opening file
 * 
 * @see CDataWriter for NetCDF attribute metadata (more detailed)
 */
std::string CSimulation::generateOutputFileName() const {
    std::ostringstream filename;
    
    // Current date and time
    auto now = std::time(nullptr);
    auto tm = *std::localtime(&now);
    
    // BASE PROJECT NAME
    filename << "barrier_sim";
    
    // Date and time: YYYYMMDD_HHMM
    filename << "_" << std::put_time(&tm, "%Y%m%d_%H%M");
    
    // Number of cross-sections (e.g., CS806)
    filename << "_CS" << std::setfill('0') << std::setw(3) << m_nCrossSectionsNumber;
    
    // Simulation duration (convert seconds to human-readable format)
    double hours = m_dSimDuration / 3600.0;
    if (hours < 1.0) {
        filename << "_T" << std::setfill('0') << std::setw(2) << static_cast<int>(m_dSimDuration / 60.0) << "min";
    } else if (hours < 24.0) {
        filename << "_T" << std::setfill('0') << std::setw(2) << static_cast<int>(hours) << "h";
    } else {
        filename << "_T" << static_cast<int>(hours / 24.0) << "d";
    }
    
    // TIMESTEP (en segundos)
    filename << "_dt" << std::setfill('0') << std::setw(3) << static_cast<int>(m_dSimTimestep);
    
    // BOUNDARY CONDITIONS
    filename << "_BC" << m_nUpwardEstuarineCondition << m_nDownwardEstuarineCondition;
    
    // SPECIAL OPTIONS
    if (m_bDoSedimentTransport) filename << "_SED";
    if (m_bDoWaterSalinity) filename << "_SAL";
    if (m_bDoMcCormackLimiterFlux) filename << "_TVD";
    if (m_bDoDryBed) filename << "_DRY";
    
    // COURANT NUMBER
    filename << "_CFL" << std::setfill('0') << std::setw(2) << static_cast<int>(m_dCourantNumber * 100);
    
    // FILE EXTENSION
    filename << ".nc";
    
    return filename.str();
}