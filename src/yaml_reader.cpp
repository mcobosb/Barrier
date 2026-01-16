/*!\n * \file yaml_reader.cpp
 * \brief YAML configuration file parser for Barrier simulation
 * \details Reads YAML config files and populates CSimulation object with:
 *          - Run parameters (timestep, duration, output variables)
 *          - Geometry file paths
 *          - Hydrodynamic settings (BC, TVD, Courant number)
 *          - Transport modules (salinity, temperature, sediment)
 *          - Smoothing options (bathymetry, solution regularization)
 * \author Manuel Cobos Budia
 * \date 2026
 * \copyright GNU General Public License
 */

#include "yaml_reader.h"
#include "simulation.h"
#include <filesystem>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <cctype>

namespace fs = std::filesystem;

/**
 * @brief Construct CYAMLReader object
 * 
 * Initializes error message string to empty.
 */
CYAMLReader::CYAMLReader() : m_strErrorMessage("") {
}

/**
 * @brief Destructor (default implementation)
 */
CYAMLReader::~CYAMLReader() {
}

/**
 * @brief Load and parse YAML configuration file
 * 
 * Workflow:
 * 1. Load YAML file using yaml-cpp library
 * 2. Extract configuration directory path for relative file references
 * 3. Parse sections in sequence:
 *    - run: Simulation parameters, output settings
 *    - geometry: File paths for bathymetry and cross-sections
 *    - hydrodynamics: Initial/boundary conditions, numerical methods
 *    - transport: Salinity, temperature, sediment modules
 *    - smoothing: Bathymetry and solution regularization
 * 4. Populate CSimulation object via setter methods
 * 
 * @param filepath Absolute or relative path to YAML config file
 * @param simulation Pointer to CSimulation object to configure
 * @return true if successful, false if YAML parsing fails
 * 
 * @note Error handling:
 * - YAML::Exception: Syntax errors, missing required fields
 * - std::exception: File not found, I/O errors
 * - Error message stored in m_strErrorMessage
 * 
 * @see Example config: docs/config_example.yaml
 */
bool CYAMLReader::loadConfiguration(const std::string& filepath, CSimulation* simulation) {
    try {
        // Load YAML file
        YAML::Node config = YAML::LoadFile(filepath);
        
        // Get input path for relative file references
        fs::path configPath(filepath);
        m_strInputPath = configPath.parent_path().string();
        if (!m_strInputPath.empty() && m_strInputPath.back() != '/') {
            m_strInputPath += "/";
        }
        
        // Parse each section
        if (config["run"]) {
            parseRunSection(config["run"], simulation);
        }
        
        if (config["geometry"]) {
            parseGeometrySection(config["geometry"]);
        }
        
        if (config["hydrodynamics"]) {
            parseHydrodynamicsSection(config["hydrodynamics"], simulation);
        }
        
        if (config["transport"]) {
            parseTransportSection(config["transport"], simulation);
        }
        
        if (config["smoothing"]) {
            parseSmoothingSection(config["smoothing"], simulation);
        }
        
        return true;
        
    } catch (const YAML::Exception& e) {
        m_strErrorMessage = "YAML parsing error: " + std::string(e.what());
        std::cerr << m_strErrorMessage << std::endl;
        return false;
    } catch (const std::exception& e) {
        m_strErrorMessage = "Error loading configuration: " + std::string(e.what());
        std::cerr << m_strErrorMessage << std::endl;
        return false;
    }
}

/**
 * @brief Parse 'run' section of YAML config
 * 
 * Fields:
 * - name: Output file basename (generates .nc and .log files)
 * - log_level: Detail level (0=minimal, 1=normal, 2=debug/all timesteps)
 * - start_date: ISO 8601 format "YYYY-MM-DDTHH:MM:SS"
 * - duration: Simulation length in seconds (mutually exclusive with end_date)
 * - end_date: ISO 8601 format "YYYY-MM-DDTHH:MM:SS" (mutually exclusive with duration)
 * - timestep: dt (seconds)
 * - output_variables: List of variables to save, or "full" for all
 * - continue_simulation: Resume from previous NetCDF (optional)
 * - continue_netcdf_path: Path to restart file (if continue=true)
 * 
 * @param node YAML node containing 'run' section
 * @param m_pSimulation Pointer to simulation object
 * 
 * @note Either 'duration' or 'end_date' must be specified (not both)
 * @note output_variables="full" includes:
 *       A, Ap, Ac, Q, Qp, Qc, Rh, B, eta, level, rho, U, c, S, T,
 *       Qb, Qs, Qt, xl, xr, UTM coordinates
 */
void CYAMLReader::parseRunSection(const YAML::Node& node, CSimulation* m_pSimulation) {
    // Output file names
    if (node["name"]) {
        m_pSimulation->m_strOutFile = node["name"].as<std::string>() + "_";
        m_pSimulation->m_strLogFile = node["name"].as<std::string>() + ".log";
    }
    
    // Log level
    if (node["log_level"]) {
        m_pSimulation->m_nLogFileDetail = node["log_level"].as<int>();
    }
    
    // Start date/time parsing: "2007-02-09T00:00:00"
    std::tm tm_start = {};
    bool has_start_date = false;
    
    if (node["start_date"]) {
        std::string dateStr = node["start_date"].as<std::string>();
        int year, month, day, hour, min, sec;
        
        // Parse ISO 8601 format: YYYY-MM-DDTHH:MM:SS
        if (sscanf(dateStr.c_str(), "%d-%d-%dT%d:%d:%d", 
                   &year, &month, &day, &hour, &min, &sec) == 6) {
            m_pSimulation->setSimulationStartDateTime(year, month, day, hour, min, sec);
            
            // Store for potential end_date calculation
            tm_start.tm_year = year - 1900;
            tm_start.tm_mon = month - 1;
            tm_start.tm_mday = day;
            tm_start.tm_hour = hour;
            tm_start.tm_min = min;
            tm_start.tm_sec = sec;
            tm_start.tm_isdst = -1;  // Let mktime determine DST
            has_start_date = true;
        }
    }
    
    // Duration in seconds (option 1: direct specification)
    if (node["duration"]) {
        m_pSimulation->dSetSimulationDuration(node["duration"].as<double>());
        
        if (node["end_date"]) {
            std::cerr << "⚠️ Warning: Both 'duration' and 'end_date' specified. Using 'duration'." << std::endl;
        }
    }
    // End date (option 2: calculate duration from start_date to end_date)
    else if (node["end_date"]) {
        if (!has_start_date) {
            std::cerr << "❌ Error: 'end_date' requires 'start_date' to be specified" << std::endl;
            exit(EXIT_FAILURE);
        }
        
        std::string endDateStr = node["end_date"].as<std::string>();
        int year, month, day, hour, min, sec;
        
        if (sscanf(endDateStr.c_str(), "%d-%d-%dT%d:%d:%d", 
                   &year, &month, &day, &hour, &min, &sec) == 6) {
            
            std::tm tm_end = {};
            tm_end.tm_year = year - 1900;
            tm_end.tm_mon = month - 1;
            tm_end.tm_mday = day;
            tm_end.tm_hour = hour;
            tm_end.tm_min = min;
            tm_end.tm_sec = sec;
            tm_end.tm_isdst = -1;
            
            // Calculate duration in seconds: difftime(end, start)
            std::time_t t_start = std::mktime(&tm_start);
            std::time_t t_end = std::mktime(&tm_end);
            
            if (t_start == -1 || t_end == -1) {
                std::cerr << "❌ Error: Invalid date format in start_date or end_date" << std::endl;
                exit(EXIT_FAILURE);
            }
            
            double duration_seconds = std::difftime(t_end, t_start);
            
            if (duration_seconds <= 0) {
                std::cerr << "❌ Error: end_date must be after start_date" << std::endl;
                exit(EXIT_FAILURE);
            }
            
            m_pSimulation->dSetSimulationDuration(duration_seconds);
            std::cout << "        - Calculated duration: " << duration_seconds 
                      << " s (" << duration_seconds / 86400.0 << " days)" << std::endl;
        } else {
            std::cerr << "❌ Error: Invalid end_date format. Expected YYYY-MM-DDTHH:MM:SS" << std::endl;
            exit(EXIT_FAILURE);
        }
    } else {
        std::cerr << "❌ Error: Either 'duration' or 'end_date' must be specified in run section" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    // Timestep in seconds
    if (node["timestep"]) {
        m_pSimulation->dSetSimulationTimestep(node["timestep"].as<double>());
    }
    
    // Output variables
    if (node["output_variables"]) {
        if (node["output_variables"].IsScalar()) {
            std::string varStr = node["output_variables"].as<std::string>();
            if (varStr == "full") {
                m_pSimulation->m_vOutputVariables = {"A", "Ap", "Ac", "Q", "Qp", "Qc", "Rh", "B", 
                                           "eta", "level", "rho", "U", "c", "S", "T", 
                                           "Qb", "Qs", "Qt", "xl", "xr", "xl_utm_x", 
                                           "xl_utm_y", "xr_utm_x", "xr_utm_y"};
            }
        } else if (node["output_variables"].IsSequence()) {
            for (const auto& var : node["output_variables"]) {
                m_pSimulation->strAddOutputVariable(var.as<std::string>());
            }
        }
    }

    // Read simulation continuation options
    if (node["continue_simulation"]) {
        m_pSimulation->m_bContinueSimulation = node["continue_simulation"].as<bool>();
    } else {
        m_pSimulation->m_bContinueSimulation = false;
    }
    if (node["continue_netcdf_path"]) {
        m_pSimulation->m_strContinueNetcdfPath = node["continue_netcdf_path"].as<std::string>();
    } else {
        m_pSimulation->m_strContinueNetcdfPath = "";
    }
}

/**
 * @brief Parse 'geometry' section of YAML config
 * 
 * Fields:
 * - along_channel_file: CSV with x, Z, Manning's n (main channel axis data)
 * - cross_sections_file: CSV with hydraulic properties vs elevation
 * 
 * @param node YAML node containing 'geometry' section
 * 
 * @note Filenames are relative to config file directory
 * @see CDataReader for CSV parsing details
 */
void CYAMLReader::parseGeometrySection(const YAML::Node& node) {
    if (node["along_channel_file"]) {
        m_strAlongChannelDataFilename = m_strInputPath + node["along_channel_file"].as<std::string>();
    }
    
    if (node["cross_sections_file"]) {
        m_strCrossSectionGeometryFilename = m_strInputPath + node["cross_sections_file"].as<std::string>();
    }
}


/**
 * @brief Parse 'hydrodynamics' section of YAML config
 * 
 * Subsections:
 * 
 * 1. initial_conditions:
 *    - type: 0=calm, 1=flow, 2=elevation
 * 
 * 2. boundary_conditions:
 *    - upstream: {type: 0=free/1=reflective, file: Q(t) or h(t), value: constant}
 *    - downstream: {type: 0=free/1=reflective/2=elevation, file: data, value: constant}
 * 
 * 3. forcing:
 *    - lateral_inflows: CSV with tributary discharge time series
 * 
 * 4. Numerical parameters:
 *    - courant_number: CFL condition (typically 0.15-0.5)
 *    - tvd_limiter: {enabled, method: minmod/roe/vanleer/vanalbada, transport_method (optional), psi_formula, delta}
 *    - surface_gradient_method: Use improved gradient calculation
 *    - source_term_balance: Balance source terms with flux gradients
 *    - beta_coefficient: Use momentum correction factor β
 *    - dry_bed: Enable dry bed treatment
 *    - murillo_condition: Use Murillo wetting/drying condition
 *    - manning_eta: Manning's n depends on water level
 *    - manning_db: Threshold for width gradient detection
 * 
 * @param node YAML node containing 'hydrodynamics' section
 * @param m_pSimulation Pointer to simulation object
 * 
 * @note BC type can be specified as int or string for readability
 */
void CYAMLReader::parseHydrodynamicsSection(const YAML::Node& node, CSimulation* m_pSimulation) {
    
    const auto& initial_conditions = node["initial_conditions"];
    // Initial estuarine condition
    if (initial_conditions["type"]) {
        if (initial_conditions["type"].IsScalar()) {
            // Try to parse as integer first
            try {
                int type = initial_conditions["type"].as<int>();
                m_pSimulation->nSetInitialEstuarineCondition(type);
            } catch (...) {
                // Parse as string
                std::string type = initial_conditions["type"].as<std::string>();
                if (type == "calm") {
                    m_pSimulation->nSetInitialEstuarineCondition(0);
                } else if (type == "flow") {
                    m_pSimulation->nSetInitialEstuarineCondition(1);
                } else if (type == "elevation") {
                    m_pSimulation->nSetInitialEstuarineCondition(2);
                }
            }
        }
    }

    // Upstream BC
    const auto& boundary_conditions = node["boundary_conditions"];
    if (boundary_conditions["upstream"]) {
        const auto& upstream = boundary_conditions["upstream"];
        if (upstream["type"]) {
            if (upstream["type"].IsScalar()) {
                try {
                    int type = upstream["type"].as<int>();
                    m_pSimulation->nSetUpwardEstuarineCondition(type);
                } catch (...) {
                    std::string type = upstream["type"].as<std::string>();
                    if (type == "free") m_pSimulation->nSetUpwardEstuarineCondition(0);
                    else if (type == "reflective") m_pSimulation->nSetUpwardEstuarineCondition(1);
                    else if (type == "elevation") m_pSimulation->nSetUpwardEstuarineCondition(2);
                    else if (type == "discharge") m_pSimulation->nSetUpwardEstuarineCondition(3);
                }
            }
        }
        if (upstream["file"] && !upstream["file"].IsNull()) {
            std::string filename = upstream["file"].as<std::string>();
            if (!filename.empty()) {
                m_pSimulation->m_strUpwardBoundaryConditionFilename = m_strInputPath + filename;
            }
        }
        if (upstream["value"] && upstream["value"].IsScalar()) {
            try {
                m_pSimulation->m_dUpwardBoundaryValue = upstream["value"].as<double>();
            } catch (...) {
                // Si no es double, ignorar
            }
        }
    }
    // Downstream BC
    if (boundary_conditions["downstream"]) {
        const auto& downstream = boundary_conditions["downstream"];
        if (downstream["type"]) {
            if (downstream["type"].IsScalar()) {
                try {
                    int type = downstream["type"].as<int>();
                    m_pSimulation->nSetDownwardEstuarineCondition(type);
                } catch (...) {
                    std::string type = downstream["type"].as<std::string>();
                    if (type == "free") m_pSimulation->nSetDownwardEstuarineCondition(0);
                    else if (type == "reflective") m_pSimulation->nSetDownwardEstuarineCondition(1);
                    else if (type == "elevation") m_pSimulation->nSetDownwardEstuarineCondition(2);
                }
            }
        }
        if (downstream["file"] && !downstream["file"].IsNull()) {
            std::string filename = downstream["file"].as<std::string>();
            if (!filename.empty()) {
                m_pSimulation->m_strDownwardBoundaryConditionFilename = m_strInputPath + filename;
            }
        }
        if (downstream["value"] && downstream["value"].IsScalar()) {
            try {
                m_pSimulation->m_dDownwardBoundaryValue = downstream["value"].as<double>();
            } catch (...) {
                // Si no es double, ignorar
            }
        }
    }

    // Support both "lateral_inflows" and "tributaries_file" field names
    const auto& forcing = node["forcing"];
    std::string filename;
    if (forcing["lateral_inflows"]) {
        filename = forcing["lateral_inflows"].as<std::string>();
    } else if (forcing["tributaries_file"]) {
        filename = forcing["tributaries_file"].as<std::string>();
    }
    
    if (!filename.empty()) {
        m_strHydrographsFilename = m_strInputPath + filename;
        m_pSimulation->m_bHydroFile = true;
    }
    
    // Courant number
    if (node["courant_number"]) {
        m_pSimulation->dSetCourantNumber(node["courant_number"].as<double>());
    }

    // Lateral storage (optional)
    // Accept either:
    //   lateral_storage: true/false
    // or
    //   lateral_storage: { enabled: true/false }
    if (node["lateral_storage"]) {
        const auto& ls = node["lateral_storage"];
        try {
            if (ls.IsScalar()) {
                m_pSimulation->bSetDoLateralStorage(ls.as<bool>());
            } else if (ls["enabled"]) {
                m_pSimulation->bSetDoLateralStorage(ls["enabled"].as<bool>());
            }
        } catch (...) {
            // Ignore malformed lateral_storage
        }
    }
    
    // TVD limiter
    if (node["tvd_limiter"]) {
        const auto& tvd = node["tvd_limiter"];
        
        if (tvd["enabled"]) {
            m_pSimulation->bSetDoMcComarckLimiterFlux(tvd["enabled"].as<bool>());
        }
        
        if (tvd["method"]) {
            std::string method = tvd["method"].as<std::string>();
            if (method == "minmod") m_pSimulation->nSetEquationLimiterFlux(1);
            else if (method == "roe") m_pSimulation->nSetEquationLimiterFlux(2);
            else if (method == "vanleer") m_pSimulation->nSetEquationLimiterFlux(3);
            else if (method == "vanalbada") m_pSimulation->nSetEquationLimiterFlux(4);
            std::cout << "        - Hydrodynamics limiter: " << method << std::endl;
        }
        
        // Optional: separate limiter for transport (salinity/temperature)
        // If not specified, transport uses same limiter as hydrodynamics
        if (tvd["transport_method"]) {
            std::string method = tvd["transport_method"].as<std::string>();
            int transport_limiter = 0;
            bool use_independent = true;
            
            if (method == "none" || method == "upwind") {
                // No TVD limiter: simple first-order upwind (very stable, diffusive)
                transport_limiter = 0;
                use_independent = true;
            }
            else if (method == "minmod") transport_limiter = 1;
            else if (method == "roe") transport_limiter = 2;
            else if (method == "vanleer") transport_limiter = 3;
            else if (method == "vanalbada") transport_limiter = 4;
            
            m_pSimulation->nSetTransportLimiterFlux(transport_limiter, use_independent);
            
            if (transport_limiter == 0) {
                std::cout << "        - Transport limiter: NONE (1st-order upwind)" << std::endl;
            } else {
                std::cout << "        - Transport limiter: " << method << std::endl;
            }
        } else {
            // Use same limiter for transport as for hydrodynamics
            m_pSimulation->nSetTransportLimiterFlux(0, false);
            std::cout << "        - Transport limiter: same as hydrodynamics" << std::endl;
        }
        
        if (tvd["psi_formula"]) {
            std::string psi = tvd["psi_formula"].as<std::string>();
            if (psi == "garcia_navarro") m_pSimulation->nSetPsiFormula(1);
            else if (psi == "tseng") m_pSimulation->nSetPsiFormula(2);
        }
        
        if (tvd["delta"]) {
            m_pSimulation->dSetDeltaValue(tvd["delta"].as<double>());
        }
    }
    
    // Other numeric methods
    if (node["surface_gradient_method"]) {
        bool use_sgm = node["surface_gradient_method"].as<bool>();
        m_pSimulation->bSetDoSurfaceGradientMethod(use_sgm);
        if (use_sgm) {
            std::cout << "      ⚠ WARNING: Surface gradient method enabled" << std::endl;
            std::cout << "         This method may produce spurious oscillations in channels" << std::endl;
            std::cout << "         with abrupt geometry transitions. Use with caution." << std::endl;
        }
    }
    
    if (node["source_term_balance"]) {
        m_pSimulation->bSetDoSurfaceTermBalance(node["source_term_balance"].as<bool>());
    }
    
    if (node["beta_coefficient"]) {
        m_pSimulation->bSetDoBetaCoefficient(node["beta_coefficient"].as<bool>());
    }

    if (node["dry_bed"]) {
        m_pSimulation->bSetDoDryBed(node["dry_bed"].as<bool>());
    }

    if (node["murillo_condition"]) {
        m_pSimulation->bSetDoMurilloCondition(node["murillo_condition"].as<bool>());
    }

    // Manning dependence on water level
    if (node["manning_eta"]) {
        m_pSimulation->bSetManningDependsOnLevel(node["manning_eta"].as<bool>());
    }
}

/**
 * @brief Parse 'transport' section of YAML config
 * 
 * Subsections:
 * 
 * 1. temperature:
 *    - enabled: Activate temperature transport
 *    - initial_file: Initial T(x) distribution CSV
 *    - dispersion_kh: Thermal diffusivity (m²/s)
 *    - beta: Thermal expansion coefficient (1/°C, ~0.0002)
 *    - boundary_conditions:
 *      * upstream: {type: 0=free/1=constant/2=timeseries/3=0D_model, value/file}
 *      * downstream: {type, value/file}
 *    - heat_flux: {heat_flux_file: meteorological forcing, cs/cl: bulk transfer coefficients}
 * 
 * 2. salinity:
 *    - enabled: Activate salinity transport
 *    - initial_file: Initial S(x) distribution CSV
 *    - dispersion_kh: Salt diffusivity (m²/s)
 *    - beta: Haline contraction coefficient (1/psu, ~0.00078)
 *    - boundary_conditions:
 *      * upstream: {type: 0=free/1=null/2=ocean, value/file}
 *      * downstream: {type, value/file}
 * 
 * 3. sediment:
 *    - enabled: Activate sediment transport
 *    - properties_file: D50, rho_s, sigma per cross-section
 *    - equation: Transport formula (0=van Rijn)
 * 
 * 4. density:
 *    - enabled: Include baroclinic pressure gradient
 * 
 * @param node YAML node containing 'transport' section
 * @param m_pSimulation Pointer to simulation object
 * 
 * @note 0D temperature model (type=3 upstream) uses reservoir energy balance:
 *       dT/dt = kA·(T_air - T) + kQ·(T_inflow - T)
 */
void CYAMLReader::parseTransportSection(const YAML::Node& node, CSimulation* m_pSimulation) {

        // Temperature transport
        if (node["temperature"]) {
            const auto& temperature = node["temperature"];
            if (temperature["enabled"]) 
            {
                const YAML::Node enabledNode = temperature["enabled"];

                if (enabledNode.IsScalar())
                {
                    std::string enabledStr = enabledNode.as<std::string>();
                    std::transform(enabledStr.begin(), enabledStr.end(), enabledStr.begin(),
                                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

                    if (enabledStr == "given")
                    {
                        m_pSimulation->m_bDoWaterTemperature = true;
                        m_pSimulation->m_eTemperatureMode = CSimulation::ETemperatureMode::Given;
                    }
                    else
                    {
                        // Backwards compatible: allow true/false
                        try
                        {
                            const bool enabled = enabledNode.as<bool>();
                            m_pSimulation->m_bDoWaterTemperature = enabled;
                            m_pSimulation->m_eTemperatureMode = enabled ? CSimulation::ETemperatureMode::Transport
                                                                        : CSimulation::ETemperatureMode::Off;
                        }
                        catch (...)
                        {
                            m_pSimulation->m_bDoWaterTemperature = false;
                            m_pSimulation->m_eTemperatureMode = CSimulation::ETemperatureMode::Off;
                        }
                    }
                }
                else
                {
                    m_pSimulation->m_bDoWaterTemperature = enabledNode.as<bool>();
                    m_pSimulation->m_eTemperatureMode = m_pSimulation->m_bDoWaterTemperature ? CSimulation::ETemperatureMode::Transport
                                                                                               : CSimulation::ETemperatureMode::Off;
                }
            }
            if (temperature["initial_file"]) {
                std::string filename = temperature["initial_file"].as<std::string>();
                if (!filename.empty()) {
                    m_pSimulation->m_strInitialTemperatureConditionFilename = m_strInputPath + filename;
                }
            }
            if (temperature["dispersion_kh"]) {
                m_pSimulation->m_dThermalDispersion = temperature["dispersion_kh"].as<double>();
            }
            // Read beta_temperature or betaT if present
            if (temperature["beta"]) {
                m_pSimulation->dSetBetaTemperatureConstant(temperature["beta"].as<double>());
            }
            if (temperature["given_file"]) {
                std::string filename = temperature["given_file"].as<std::string>();
                if (!filename.empty()) {
                    m_pSimulation->m_strGivenTemperatureFilename = m_strInputPath + filename;
                }
            }
            else if (temperature["file"]) {
                std::string filename = temperature["file"].as<std::string>();
                if (!filename.empty()) {
                    m_pSimulation->m_strGivenTemperatureFilename = m_strInputPath + filename;
                }
            }

            if (temperature["boundary_conditions"]) {
                const auto& bc = temperature["boundary_conditions"];
                // UPSTREAM
                if (bc["upstream"]) {
                    const auto& up = bc["upstream"];
                    if (up["type"]) {
                        try { m_pSimulation->m_nUpwardTemperatureCondition = up["type"].as<int>(); }
                        catch (...) {}
                    }
                    if (m_pSimulation->m_nUpwardTemperatureCondition == 1) {
                        if (up["value"] && up["value"].IsScalar()) {
                            try { m_pSimulation->m_dUpwardTemperatureBoundaryValue = up["value"].as<double>(); } catch (...) {}
                        }
                    }
                    if (m_pSimulation->m_nUpwardTemperatureCondition == 2) {
                        // Series temporal
                        if (up["file"] && !up["file"].IsNull()) {
                            std::string filename = up["file"].as<std::string>();
                            if (!filename.empty()) m_pSimulation->m_strUpwardTemperatureBoundaryConditionFilename = m_strInputPath + filename;
                        }
                    }
                    if (m_pSimulation->m_nUpwardTemperatureCondition == 3) {
                        // Model 0D a reservoir
                        if (up["offset_beta"] && up["offset_beta"].IsScalar()) {
                            try { m_pSimulation->m_dUpwardTemperatureOffsetBeta = up["offset_beta"].as<double>(); } catch (...) {}
                        }
                        if (up["atmospheric_exchange_kA"] && up["atmospheric_exchange_kA"].IsScalar()) {
                            try { m_pSimulation->m_dUpwardAtmosphericExchangekA = up["atmospheric_exchange_kA"].as<double>(); } catch (...) {}
                        }
                        if (up["inflow_water_effect_kQ"] && up["inflow_water_effect_kQ"].IsScalar()) {
                            try { m_pSimulation->m_dUpwardInflowWaterEffectkQ = up["inflow_water_effect_kQ"].as<double>(); } catch (...) {}
                        }
                        // Series temporal
                        if (up["file"] && !up["file"].IsNull()) {
                            std::string filename = up["file"].as<std::string>();
                            if (!filename.empty()) m_pSimulation->m_strUpwardTemperatureBoundaryConditionFilename = m_strInputPath + filename;
                        }
                    }
                }
                // DOWNSTREAM
                if (bc["downstream"]) {
                    const auto& down = bc["downstream"];
                    if (down["type"]) {
                        try { m_pSimulation->m_nDownwardTemperatureCondition = down["type"].as<int>(); }
                        catch (...) {}
                    }
                    if (down["file"] && !down["file"].IsNull()) {
                        std::string filename = down["file"].as<std::string>();
                        if (!filename.empty()) m_pSimulation->m_strDownwardTemperatureBoundaryConditionFilename = m_strInputPath + filename;
                    }
                    if (down["value"] && down["value"].IsScalar()) {
                        try { m_pSimulation->m_dDownwardTemperatureBoundaryValue = down["value"].as<double>(); } catch (...) {}
                    }
                }
            }
            if (temperature["heat_flux"]) {
                const auto& hf = temperature["heat_flux"];
                if (hf["heat_flux_file"]) {
                    std::string filename = hf["heat_flux_file"].as<std::string>();
                    if (!filename.empty()) m_pSimulation->m_strHeatFluxFile = m_strInputPath + filename;
                }
                // Bulk transfer coefficients
                if (hf["cs"]) {
                    m_pSimulation->m_dHeatFlux_CS = hf["cs"].as<double>();
                }
                if (hf["cl"]) {
                    m_pSimulation->m_dHeatFlux_CL = hf["cl"].as<double>();
                }
                // Calibration parameters for radiative flux calculations
                if (hf["latitude"]) {
                    m_pSimulation->m_dHeatFluxLatitude = hf["latitude"].as<double>();
                }
                if (hf["cloud_cover"]) {
                    m_pSimulation->m_dHeatFluxCloudCover = hf["cloud_cover"].as<double>();
                }
                // Effective depth for 0D reservoir temperature model
                if (hf["reservoir_depth"]) {
                    m_pSimulation->m_dReservoirEffectiveDepth = hf["reservoir_depth"].as<double>();
                }
                // Option to calculate RH from temperature using FAO-56 method
                if (hf["calculate_rh_from_temperature"]) {
                    m_pSimulation->m_bCalculateRHFromTemperature = hf["calculate_rh_from_temperature"].as<bool>();
                }
            }
        }
    // Salinity transport
    if (node["salinity"]) {
        const auto& salinity = node["salinity"];
        
        if (salinity["enabled"]) {
            m_pSimulation->bSetDoWaterSalinity(salinity["enabled"].as<bool>());
        }
        
        if (salinity["initial_file"]) {
            std::string filename = salinity["initial_file"].as<std::string>();
            if (!filename.empty()) {
                m_strSalinityFilename = m_strInputPath + filename;
            }
        }
        
        if (salinity["boundary_conditions"]) {
            const auto& bc = salinity["boundary_conditions"];
            // UPSTREAM
            if (bc["upstream"]) {
                const auto& up = bc["upstream"];
                if (up["type"]) {
                    try { m_pSimulation->nSetUpwardSalinityCondition(up["type"].as<int>()); }
                    catch (...) {
                        std::string t = up["type"].as<std::string>();
                        if (t == "free") m_pSimulation->nSetUpwardSalinityCondition(0);
                        else if (t == "null") m_pSimulation->nSetUpwardSalinityCondition(1);
                        else if (t == "ocean") m_pSimulation->nSetUpwardSalinityCondition(2);
                    }
                }
                if (up["file"] && !up["file"].IsNull()) {
                    std::string filename = up["file"].as<std::string>();
                    if (!filename.empty()) m_pSimulation->m_strUpwardSalinityBoundaryConditionFilename = m_strInputPath + filename;
                }
                if (up["value"] && up["value"].IsScalar()) {
                    try { m_pSimulation->m_dUpwardSalinityBoundaryValue = up["value"].as<double>(); } catch (...) {}
                }
            }
            // DOWNSTREAM
            if (bc["downstream"]) {
                const auto& down = bc["downstream"];
                if (down["type"]) {
                    try { m_pSimulation->nSetDownwardSalinityCondition(down["type"].as<int>()); }
                    catch (...) {
                        std::string t = down["type"].as<std::string>();
                        if (t == "free") m_pSimulation->nSetDownwardSalinityCondition(0);
                        else if (t == "null") m_pSimulation->nSetDownwardSalinityCondition(1);
                        else if (t == "ocean") m_pSimulation->nSetDownwardSalinityCondition(2);
                    }
                }
                if (down["file"] && !down["file"].IsNull()) {
                    std::string filename = down["file"].as<std::string>();
                    if (!filename.empty()) m_pSimulation->m_strDownwardSalinityBoundaryConditionFilename = m_strInputPath + filename;
                }
                if (down["value"] && down["value"].IsScalar()) {
                    try { m_pSimulation->m_dDownwardSalinityBoundaryValue = down["value"].as<double>(); } catch (...) {}
                }
            }
        }                
        
        if (salinity["beta"]) {
            m_pSimulation->dSetBetaSalinityConstant(salinity["beta"].as<double>());
        }

        if (salinity["dispersion_kh"]) {
            m_pSimulation->dSetLongitudinalDispersionConstant(salinity["dispersion_kh"].as<double>());
        }
    }

    // Sediment transport
    if (node["sediment"]) {
        const auto& sediment = node["sediment"];
        if (sediment["enabled"]) m_pSimulation->bSetDoSedimentTransport(sediment["enabled"].as<bool>());
        if (sediment["properties_file"]) m_strSedimentPropertiesFilename = m_strInputPath + sediment["properties_file"].as<std::string>();
        if (sediment["equation"]) {
            int eq = 0;
            try { eq = sediment["equation"].as<int>(); } catch (...) {}
            m_pSimulation->nSetEquationSedimentTransport(eq);
        }
    }

    // Water density
    if (node["density"]) {
        const auto& density = node["density"];
        if (density["enabled"]) m_pSimulation->bSetDoWaterDensity(density["enabled"].as<bool>());
    }
}

/**
 * @brief Parse 'smoothing' section of YAML config
 * 
 * Subsections:
 * 
 * 1. bathymetry:
 *    - enabled: Apply pre-simulation smoothing to bed elevation
 *    - num_passes: Number of diffusion iterations (typically 1-5)
 *    - alpha: Diffusion coefficient (0.1-0.5, larger = more smoothing)
 * 
 * 2. regularization:
 *    - enabled: Apply runtime solution smoothing
 *    - num_passes: Passes per timestep (typically 1-2)
 *    - alpha: Smoothing strength (0.05-0.2)
 * 
 * Purpose:
 * - Bathymetry smoothing: Remove survey noise, prevent numerical oscillations
 * - Solution regularization: Stabilize at sharp geometry transitions
 * 
 * @param node YAML node containing 'smoothing' section
 * @param m_pSimulation Pointer to simulation object
 * 
 * @warning Excessive smoothing can:
 * - Remove important morphological features
 * - Introduce artificial diffusion
 * - Violate conservation laws
 * 
 * @note Default: Both disabled (alpha=0, num_passes=0)
 */
void CYAMLReader::parseSmoothingSection(const YAML::Node& node, CSimulation* m_pSimulation) {
    if (node["bathymetry"]) {
        const auto& bathy = node["bathymetry"];
        if (bathy["enabled"]) {
            m_pSimulation->bSetDoSmoothBathymetry(bathy["enabled"].as<bool>());
        }
        if (bathy["num_passes"]) {
            m_pSimulation->m_nBathymetrySmoothingPasses = bathy["num_passes"].as<int>();
        }
        if (bathy["alpha"]) {
            m_pSimulation->m_dBathymetrySmoothingAlpha = bathy["alpha"].as<double>();
        }
    }
    if (node["regularization"]) {
        const auto& reg = node["regularization"];
        if (reg["enabled"]) {
            m_pSimulation->bSetDoSmoothSolution(reg["enabled"].as<bool>());
        }
        if (reg["num_passes"]) {
            m_pSimulation->m_nSolutionSmoothingPasses = reg["num_passes"].as<int>();
        }
        if (reg["alpha"]) {
            m_pSimulation->m_dSolutionSmoothingAlpha = reg["alpha"].as<double>();
        }
    }
}