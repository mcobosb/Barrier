/*!
 * \file yaml_reader.cpp
 * \brief Implementation of CYAMLReader class
 */

#include "yaml_reader.h"
#include "simulation.h"
#include <iostream>
#include <filesystem>
#include <sstream>

namespace fs = std::filesystem;

CYAMLReader::CYAMLReader() : m_strErrorMessage("") {
}

CYAMLReader::~CYAMLReader() {
}

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
        
        if (config["initial_conditions"]) {
            parseInitialConditionsSection(config["initial_conditions"], simulation);
        }
        
        if (config["boundary_conditions"]) {
            parseBoundaryConditionsSection(config["boundary_conditions"], simulation);
        }
        
        if (config["forcing"]) {
            parseForcingSection(config["forcing"], simulation);
        }
        
        if (config["numerics"]) {
            parseNumericsSection(config["numerics"], simulation);
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
    if (node["start_date"]) {
        std::string dateStr = node["start_date"].as<std::string>();
        int year, month, day, hour, min, sec;
        
        // Parse ISO 8601 format: YYYY-MM-DDTHH:MM:SS
        if (sscanf(dateStr.c_str(), "%d-%d-%dT%d:%d:%d", 
                   &year, &month, &day, &hour, &min, &sec) == 6) {
            m_pSimulation->setSimulationStartDateTime(year, month, day, hour, min, sec);
        }
    }
    
    // Duration in seconds
    if (node["duration"]) {
        m_pSimulation->dSetSimulationDuration(node["duration"].as<double>());
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
                                           "eta", "level", "rho", "U", "c", "S", 
                                           "Qb", "Qs", "Qt", "xl", "xr", "xl_utm_x", 
                                           "xl_utm_y", "xr_utm_x", "xr_utm_y"};
            }
        } else if (node["output_variables"].IsSequence()) {
            for (const auto& var : node["output_variables"]) {
                m_pSimulation->strAddOutputVariable(var.as<std::string>());
            }
        }
    }

    // Leer opciones de continuación de simulación
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

void CYAMLReader::parseGeometrySection(const YAML::Node& node) {
    if (node["along_channel_file"]) {
        m_strAlongChannelDataFilename = m_strInputPath + node["along_channel_file"].as<std::string>();
    }
    
    if (node["cross_sections_file"]) {
        m_strCrossSectionGeometryFilename = m_strInputPath + node["cross_sections_file"].as<std::string>();
    }
}

void CYAMLReader::parseInitialConditionsSection(const YAML::Node& node, CSimulation* m_pSimulation) {
    if (node["type"]) {
        if (node["type"].IsScalar()) {
            // Try to parse as integer first
            try {
                int type = node["type"].as<int>();
                m_pSimulation->nSetInitialEstuarineCondition(type);
            } catch (...) {
                // Parse as string
                std::string type = node["type"].as<std::string>();
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
}

void CYAMLReader::parseBoundaryConditionsSection(const YAML::Node& node, CSimulation* m_pSimulation) {
    // Upstream BC
    // Upstream BC
    if (node["upstream"]) {
        const auto& upstream = node["upstream"];
        if (upstream["type"]) {
            if (upstream["type"].IsScalar()) {
                try {
                    int type = upstream["type"].as<int>();
                    m_pSimulation->nSetUpwardEstuarineCondition(type);
                } catch (...) {
                    std::string type = upstream["type"].as<std::string>();
                    if (type == "reflective") m_pSimulation->nSetUpwardEstuarineCondition(1);
                    else if (type == "free") m_pSimulation->nSetUpwardEstuarineCondition(0);
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
    if (node["downstream"]) {
        const auto& downstream = node["downstream"];
        if (downstream["type"]) {
            if (downstream["type"].IsScalar()) {
                try {
                    int type = downstream["type"].as<int>();
                    m_pSimulation->nSetDownwardEstuarineCondition(type);
                } catch (...) {
                    std::string type = downstream["type"].as<std::string>();
                    if (type == "elevation") m_pSimulation->nSetDownwardEstuarineCondition(2);
                    else if (type == "reflective") m_pSimulation->nSetDownwardEstuarineCondition(1);
                    else if (type == "free") m_pSimulation->nSetDownwardEstuarineCondition(0);
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
}


void CYAMLReader::parseForcingSection(const YAML::Node& node, CSimulation* m_pSimulation) {
    // Support both "lateral_inflows" and "tributaries_file" field names
    std::string filename;
    if (node["lateral_inflows"]) {
        filename = node["lateral_inflows"].as<std::string>();
    } else if (node["tributaries_file"]) {
        filename = node["tributaries_file"].as<std::string>();
    }
    
    if (!filename.empty()) {
        m_strHydrographsFilename = m_strInputPath + filename;
        m_pSimulation->m_bHydroFile = true;
    }
}

void CYAMLReader::parseNumericsSection(const YAML::Node& node, CSimulation* m_pSimulation) {
    // Courant number
    if (node["courant_number"]) {
        m_pSimulation->dSetCourantNumber(node["courant_number"].as<double>());
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
        m_pSimulation->bSetDoSurfaceGradientMethod(node["surface_gradient_method"].as<bool>());
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
}

void CYAMLReader::parseTransportSection(const YAML::Node& node, CSimulation* m_pSimulation) {

        // Temperature transport
        if (node["temperature"]) {
            const auto& temperature = node["temperature"];
            if (temperature["enabled"]) {
                m_pSimulation->m_bDoWaterTemperature = temperature["enabled"].as<bool>();
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
            // Leer beta_temperature o betaT si está presente
            if (temperature["beta"]) {
                m_pSimulation->dSetBetaTemperatureConstant(temperature["beta_temperature"].as<double>());
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
                    if (up["file"] && !up["file"].IsNull()) {
                        std::string filename = up["file"].as<std::string>();
                        if (!filename.empty()) m_pSimulation->m_strUpwardTemperatureBoundaryConditionFilename = m_strInputPath + filename;
                    }
                    if (up["value"] && up["value"].IsScalar()) {
                        try { m_pSimulation->m_dUpwardTemperatureBoundaryValue = up["value"].as<double>(); } catch (...) {}
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
                if (hf["cs"]) {
                    m_pSimulation->m_dHeatFlux_CS = hf["cs"].as<double>();
                }
                if (hf["cl"]) {
                    m_pSimulation->m_dHeatFlux_CL = hf["cl"].as<double>();
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