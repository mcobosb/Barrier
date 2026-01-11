/*!
 * \class CYAMLReader
 * \brief This class reads YAML configuration files
 * \details Replaces the old .conf format with modern YAML
 * \author Manuel Cobos Budia
 * \date 2025
 * \copyright GNU General Public License
 * 
 * \file yaml_reader.h
 * \brief Contains CYAMLReader definitions
 */

#ifndef YAML_READER_H
#define YAML_READER_H

#include <string>
#include <vector>
#include <yaml-cpp/yaml.h>

class CSimulation;

class CYAMLReader {
public:
    CYAMLReader();
    ~CYAMLReader();
    
    //! Load and parse YAML configuration file
    bool loadConfiguration(const std::string& filepath, CSimulation* simulation);
    
    //! Get error message if loading failed
    std::string getErrorMessage() const { return m_strErrorMessage; }
    
private:
    //! Parse run section
    void parseRunSection(const YAML::Node& node, CSimulation* sim);
    
    //! Parse geometry section  
    void parseGeometrySection(const YAML::Node& node);
    
    //! Parse hydrodynamics section
    void parseHydrodynamicsSection(const YAML::Node& node, CSimulation* sim);
    
    //! Parse transport section
    void parseTransportSection(const YAML::Node& node, CSimulation* sim);
    
    //! Parse smoothing section
    void parseSmoothingSection(const YAML::Node& node, CSimulation* sim);
    
    //! Error message storage
    std::string m_strErrorMessage;
    
    //! Input path for relative file references
    std::string m_strInputPath;
    
public:
    //! File paths extracted from YAML (to be used by data_reader)
    std::string m_strAlongChannelDataFilename;
    std::string m_strCrossSectionGeometryFilename;
    std::string m_strUpwardBoundaryConditionFilename;
    std::string m_strDownwardBoundaryConditionFilename;
    std::string m_strUpwardSalinityBoundaryConditionFilename;
    std::string m_strDownwardSalinityBoundaryConditionFilename;
    std::string m_strSedimentPropertiesFilename;
    std::string m_strHydrographsFilename;
    std::string m_strSalinityFilename;
};

#endif // YAML_READER_H
