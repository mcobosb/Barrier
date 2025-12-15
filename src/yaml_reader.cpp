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
            parseGeometrySection(config["geometry"], simulation);
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

void CYAMLReader::parseRunSection(const YAML::Node& node, CSimulation* sim) {
    // Output file names
    if (node["name"]) {
        sim->m_strOutFile = node["name"].as<std::string>() + "_";
        sim->m_strLogFile = node["name"].as<std::string>() + ".log";
    }
    
    // Log level
    if (node["log_level"]) {
        sim->m_nLogFileDetail = node["log_level"].as<int>();
    }
    
    // Start date/time parsing: "2007-02-09T00:00:00"
    if (node["start_date"]) {
        std::string dateStr = node["start_date"].as<std::string>();
        int year, month, day, hour, min, sec;
        
        // Parse ISO 8601 format: YYYY-MM-DDTHH:MM:SS
        if (sscanf(dateStr.c_str(), "%d-%d-%dT%d:%d:%d", 
                   &year, &month, &day, &hour, &min, &sec) == 6) {
            sim->setSimulationStartDateTime(year, month, day, hour, min, sec);
        }
    }
    
    // Duration in seconds
    if (node["duration"]) {
        sim->dSetSimulationDuration(node["duration"].as<double>());
    }
    
    // Timestep in seconds
    if (node["timestep"]) {
        sim->dSetSimulationTimestep(node["timestep"].as<double>());
    }
    
    // Output variables
    if (node["output_variables"]) {
        if (node["output_variables"].IsScalar()) {
            std::string varStr = node["output_variables"].as<std::string>();
            if (varStr == "full") {
                sim->m_vOutputVariables = {"A", "Ap", "Ac", "Q", "Qp", "Qc", "Rh", "B", 
                                           "eta", "level", "beta", "rho", "U", "c", "S", 
                                           "Qb", "Qs", "Qt", "xl", "xr", "xl_utm_x", 
                                           "xl_utm_y", "xr_utm_x", "xr_utm_y"};
            }
        } else if (node["output_variables"].IsSequence()) {
            for (const auto& var : node["output_variables"]) {
                sim->strAddOutputVariable(var.as<std::string>());
            }
        }
    }

    // Leer opciones de continuación de simulación
    if (node["continue_simulation"]) {
        sim->m_bContinueSimulation = node["continue_simulation"].as<bool>();
    } else {
        sim->m_bContinueSimulation = false;
    }
    if (node["continue_netcdf_path"]) {
        sim->m_strContinueNetcdfPath = node["continue_netcdf_path"].as<std::string>();
    } else {
        sim->m_strContinueNetcdfPath = "";
    }
}

void CYAMLReader::parseGeometrySection(const YAML::Node& node, CSimulation* sim) {
    if (node["along_channel_file"]) {
        m_strAlongChannelDataFilename = m_strInputPath + node["along_channel_file"].as<std::string>() + ".csv";
    }
    
    if (node["cross_sections_file"]) {
        m_strCrossSectionGeometryFilename = m_strInputPath + node["cross_sections_file"].as<std::string>() + ".csv";
    }
    
    (void)sim;
}

void CYAMLReader::parseInitialConditionsSection(const YAML::Node& node, CSimulation* sim) {
    if (node["type"]) {
        if (node["type"].IsScalar()) {
            // Try to parse as integer first
            try {
                int type = node["type"].as<int>();
                sim->nSetInitialEstuarineCondition(type);
            } catch (...) {
                // Parse as string
                std::string type = node["type"].as<std::string>();
                if (type == "calm") {
                    sim->nSetInitialEstuarineCondition(0);
                } else if (type == "flow") {
                    sim->nSetInitialEstuarineCondition(1);
                } else if (type == "elevation") {
                    sim->nSetInitialEstuarineCondition(2);
                }
            }
        }
    }
}

void CYAMLReader::parseBoundaryConditionsSection(const YAML::Node& node, CSimulation* sim) {
    // Upstream BC
    if (node["upstream"]) {
        const auto& upstream = node["upstream"];
        if (upstream["type"]) {
            // Try to parse as integer first, then as string
            try {
                int type = upstream["type"].as<int>();
                sim->nSetUpwardEstuarineCondition(type);
            } catch (...) {
                std::string type = upstream["type"].as<std::string>();
                if (type == "open") {
                    sim->nSetUpwardEstuarineCondition(0);
                } else if (type == "reflective") {
                    sim->nSetUpwardEstuarineCondition(1);
                } else if (type == "elevation") {
                    sim->nSetUpwardEstuarineCondition(2);
                }
            }
        }
        if (upstream["file"] && !upstream["file"].IsNull()) {
            std::string filename = upstream["file"].as<std::string>();
            if (!filename.empty()) {
                m_strUpwardBoundaryConditionFilename = m_strInputPath + filename + ".csv";
            }
        }
    }
    
    // Downstream BC
    if (node["downstream"]) {
        const auto& downstream = node["downstream"];
        if (downstream["type"]) {
            // Try to parse as integer first, then as string
            try {
                int type = downstream["type"].as<int>();
                sim->nSetDownwardEstuarineCondition(type);
            } catch (...) {
                std::string type = downstream["type"].as<std::string>();
                if (type == "open") {
                    sim->nSetDownwardEstuarineCondition(0);
                } else if (type == "reflective") {
                    sim->nSetDownwardEstuarineCondition(1);
                } else if (type == "elevation") {
                    sim->nSetDownwardEstuarineCondition(2);
                }
            }
        }
        if (downstream["file"] && !downstream["file"].IsNull()) {
            std::string filename = downstream["file"].as<std::string>();
            if (!filename.empty()) {
                m_strDownwardBoundaryConditionFilename = m_strInputPath + filename + ".csv";
            }
        }
    }
}

void CYAMLReader::parseForcingSection(const YAML::Node& node, CSimulation* sim) {
    // Support both "lateral_inflows" and "tributaries_file" field names
    std::string filename;
    if (node["lateral_inflows"]) {
        filename = node["lateral_inflows"].as<std::string>();
    } else if (node["tributaries_file"]) {
        filename = node["tributaries_file"].as<std::string>();
    }
    
    if (!filename.empty()) {
        m_strHydrographsFilename = m_strInputPath + filename + ".csv";
        sim->m_bHydroFile = true;
    }
}

void CYAMLReader::parseNumericsSection(const YAML::Node& node, CSimulation* sim) {
    // Courant number
    if (node["courant_number"]) {
        sim->dSetCourantNumber(node["courant_number"].as<double>());
    }
    
    // TVD limiter
    if (node["tvd_limiter"]) {
        const auto& tvd = node["tvd_limiter"];
        
        if (tvd["enabled"]) {
            sim->bSetDoMcComarckLimiterFlux(tvd["enabled"].as<bool>());
        }
        
        if (tvd["method"]) {
            std::string method = tvd["method"].as<std::string>();
            if (method == "minmod") sim->nSetEquationLimiterFlux(1);
            else if (method == "roe") sim->nSetEquationLimiterFlux(2);
            else if (method == "vanleer") sim->nSetEquationLimiterFlux(3);
            else if (method == "vanalbada") sim->nSetEquationLimiterFlux(4);
        }
        
        if (tvd["psi_formula"]) {
            std::string psi = tvd["psi_formula"].as<std::string>();
            if (psi == "garcia_navarro") sim->nSetPsiFormula(1);
            else if (psi == "tseng") sim->nSetPsiFormula(2);
        }
        
        if (tvd["delta"]) {
            sim->dSetDeltaValue(tvd["delta"].as<double>());
        }
    }
    
    // Other numeric methods
    if (node["surface_gradient_method"]) {
        sim->bSetDoSurfaceGradientMethod(node["surface_gradient_method"].as<bool>());
    }
    
    if (node["source_term_balance"]) {
        sim->bSetDoSurfaceTermBalance(node["source_term_balance"].as<bool>());
    }
    
    if (node["beta_coefficient"]) {
        sim->bSetDoBetaCoefficient(node["beta_coefficient"].as<bool>());
    }
    
    if (node["dry_bed"]) {
        sim->bSetDoDryBed(node["dry_bed"].as<bool>());
    }
    
    if (node["murillo_condition"]) {
        sim->bSetDoMurilloCondition(node["murillo_condition"].as<bool>());
    }
}

void CYAMLReader::parseTransportSection(const YAML::Node& node, CSimulation* sim) {
    // Salinity transport
    if (node["salinity"]) {
        const auto& salinity = node["salinity"];
        
        if (salinity["enabled"]) {
            sim->bSetDoWaterSalinity(salinity["enabled"].as<bool>());
        }
        
        if (salinity["initial_file"]) {
            std::string filename = salinity["initial_file"].as<std::string>();
            if (!filename.empty()) {
                m_strSalinityFilename = m_strInputPath + filename + ".csv";
            }
        }
        
        // Support both "upstream_bc" and "upstream_condition"
        auto upstreamNode = salinity["upstream_bc"] ? salinity["upstream_bc"] : salinity["upstream_condition"];
        if (upstreamNode) {
            try {
                int bc = upstreamNode.as<int>();
                sim->nSetUpwardSalinityCondition(bc);
            } catch (...) {
                std::string bc = upstreamNode.as<std::string>();
                if (bc == "free") sim->nSetUpwardSalinityCondition(0);
                else if (bc == "null") sim->nSetUpwardSalinityCondition(1);
                else if (bc == "ocean") sim->nSetUpwardSalinityCondition(2);
            }
        }
        
        // Support both "downstream_bc" and "downstream_condition"
        auto downstreamNode = salinity["downstream_bc"] ? salinity["downstream_bc"] : salinity["downstream_condition"];
        if (downstreamNode) {
            try {
                int bc = downstreamNode.as<int>();
                sim->nSetDownwardSalinityCondition(bc);
            } catch (...) {
                std::string bc = downstreamNode.as<std::string>();
                if (bc == "free") sim->nSetDownwardSalinityCondition(0);
                else if (bc == "null") sim->nSetDownwardSalinityCondition(1);
                else if (bc == "ocean") sim->nSetDownwardSalinityCondition(2);
            }
        }
        
        if (salinity["beta"]) {
            sim->dSetBetaSalinityConstant(salinity["beta"].as<double>());
        }
        
        if (salinity["dispersion_kh"]) {
            sim->dSetLongitudinalDispersionConstant(salinity["dispersion_kh"].as<double>());
        }
    }
    
    // Sediment transport
    if (node["sediment"]) {
        const auto& sediment = node["sediment"];
        
        if (sediment["enabled"]) {
            sim->bSetDoSedimentTransport(sediment["enabled"].as<bool>());
        }
        
        if (sediment["properties_file"]) {
            m_strSedimentPropertiesFilename = m_strInputPath + sediment["properties_file"].as<std::string>() + ".csv";
        }
        
        if (sediment["method"]) {
            std::string method = sediment["method"].as<std::string>();
            if (method == "vanrijn") {
                sim->nSetEquationSedimentTransport(0);
            }
        }
    }
    
    // Water density
    if (node["density"]) {
        const auto& density = node["density"];
        
        if (density["enabled"]) {
            sim->bSetDoWaterDensity(density["enabled"].as<bool>());
        }
    }
}

void CYAMLReader::parseSmoothingSection(const YAML::Node& node, CSimulation* sim) {
    if (node["bathymetry"]) {
        sim->bSetDoSmoothBathymetry(node["bathymetry"].as<bool>());
    }
    
    if (node["solution"]) {
        sim->bSetDoSmoothSolution(node["solution"].as<bool>());
    }
}
