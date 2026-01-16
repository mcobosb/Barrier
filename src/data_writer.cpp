/*!
*
 * \file data_writer.cpp
 * \brief Write output NetCDF file with simulation results
 * \details Handles creation and writing of NetCDF output files containing
 *          temporal evolution of hydraulic variables (area, flow, velocity, etc.)
 *          along the estuarine channel.
 * \author Manuel Cobos Budia

 * \date 2026
 * \copyright GNU General Public License
 *
 */

/*==============================================================================================================================

==============================================================================================================================*/
#include <cstdlib>

#include <sstream>
using std::stringstream;
#include <cstring>  // Para strlen

#include <iostream>
using std::cerr;
using std::cout;
using std::endl;
using std::ios;

#include <string>
using std::to_string;
using std::stoi;
using std::string;

#include <vector>
using std::vector;

#include <algorithm>
using std::find;

#include <ctime>
using std::time;
using std::localtime;

#include <iomanip>
using std::put_time;
using std::setfill;
using std::setw;

#include "simulation.h"
#include "data_writer.h"
#include "main.h"


/**
 * @brief Construct CDataWriter and initialize variable metadata dictionary
 * 
 * Creates metadata map (m_mVariableDefinitions) with CF-compliant attributes:
 * - description: Long description of variable
 * - longname: Short display name
 * - units: Physical units (m, m/s, m², psu, °C, kg/m³, etc.)
 * 
 * Variables defined:
 * - Hydraulics: A (area), Q (discharge), U (velocity), Rh (hydraulic radius)
 * - Water surface: level (depth), eta (elevation)
 * - Transport: S (salinity), T (temperature), rho (density)
 * - Sediment: Qb (bedload), Qs (suspended), Qt (total)
 * - Geometry: xl/xr (bank distances), UTM coordinates
 * - Coefficients: beta (momentum correction), I1/I2 (pressure integrals)
 * 
 * @note CF conventions: www.cfconventions.org
 * @see nDefineNetCDFFile() for writing metadata to file
 */
CDataWriter::CDataWriter() {

    // Creating the attributes of the main variables
    m_mVariableDefinitions["A"]["description"] = "Cross-sectional area";
    m_mVariableDefinitions["A"]["longname"] = "Area";
    m_mVariableDefinitions["A"]["units"] = "m2";

    m_mVariableDefinitions["Q"]["description"] = "Cross-sectional averaged water flow";
    m_mVariableDefinitions["Q"]["longname"] = "water flow";
    m_mVariableDefinitions["Q"]["units"] = "m3/s";

    m_mVariableDefinitions["q"]["description"] = "Fluvial contributions along the estuary";
    m_mVariableDefinitions["q"]["longname"] = "water flux contributions";
    m_mVariableDefinitions["q"]["units"] = "m3/s";

    m_mVariableDefinitions["Rh"]["description"] = "Hydraulic radius";
    m_mVariableDefinitions["Rh"]["longname"] = "hydraulic radius";
    m_mVariableDefinitions["Rh"]["units"] = "m";

    m_mVariableDefinitions["n"]["description"] = "Manning roughness coefficient";
    m_mVariableDefinitions["n"]["longname"] = "manning coefficient";
    m_mVariableDefinitions["n"]["units"] = "s/m^(1/3)";

    m_mVariableDefinitions["Sf"]["description"] = "Lateral storage factor S (dimensionless, >=1)";
    m_mVariableDefinitions["Sf"]["longname"] = "lateral storage factor";
    m_mVariableDefinitions["Sf"]["units"] = "";

    m_mVariableDefinitions["level"]["description"] = "Water depth";
    m_mVariableDefinitions["level"]["longname"] = "water depth";
    m_mVariableDefinitions["level"]["units"] = "m";

    m_mVariableDefinitions["eta"]["description"] = "Free surface elevation";
    m_mVariableDefinitions["eta"]["longname"] = "free surface elevation";
    m_mVariableDefinitions["eta"]["units"] = "m";

    m_mVariableDefinitions["beta"]["description"] = "Beta coefficient";
    m_mVariableDefinitions["beta"]["longname"] = "beta";
    m_mVariableDefinitions["beta"]["units"] = "";

    m_mVariableDefinitions["I2"]["description"] = "Term I2 of the balance";
    m_mVariableDefinitions["I2"]["longname"] = "I2";
    m_mVariableDefinitions["I2"]["units"] = "";

    m_mVariableDefinitions["U"]["description"] = "Cross-sectional averaged water velocity";
    m_mVariableDefinitions["U"]["longname"] = "mean water velocity";
    m_mVariableDefinitions["U"]["units"] = "m/s";

    m_mVariableDefinitions["S"]["description"] = "Cross-sectional averaged salinity";
    m_mVariableDefinitions["S"]["longname"] = "salinity";
    m_mVariableDefinitions["S"]["units"] = "psu";

    m_mVariableDefinitions["T"]["description"] = "Cross-sectional averaged water temperature";
    m_mVariableDefinitions["T"]["longname"] = "water temperature";
    m_mVariableDefinitions["T"]["units"] = "degC";

    m_mVariableDefinitions["Qb"]["description"] = "Bed-load sediment transport";
    m_mVariableDefinitions["Qb"]["longname"] = "bed-load sediment transport";
    m_mVariableDefinitions["Qb"]["units"] = "m3/s";

    m_mVariableDefinitions["Qs"]["description"] = "Suspended sediment transport";
    m_mVariableDefinitions["Qs"]["longname"] = "suspended sediment transport";
    m_mVariableDefinitions["Qs"]["units"] = "m3/s";

    m_mVariableDefinitions["Qt"]["description"] = "Total sediment transport";
    m_mVariableDefinitions["Qt"]["longname"] = "total sediment transport";
    m_mVariableDefinitions["Qt"]["units"] = "m3/s";

    m_mVariableDefinitions["rho"]["description"] = "Fluid density (salinity, bed-load and suspended sediment)";
    m_mVariableDefinitions["rho"]["longname"] = "fluid density";
    m_mVariableDefinitions["rho"]["units"] = "kg/m3";

    m_mVariableDefinitions["xl"]["description"] = "Distance from the thalweg to the left riverbank";
    m_mVariableDefinitions["xlb"]["longname"] = "x-left";
    m_mVariableDefinitions["xl"]["units"] = "m";

    m_mVariableDefinitions["xr"]["description"] = "Distance from the thalweg to the right riverbank";
    m_mVariableDefinitions["xr"]["longname"] = "x-right";
    m_mVariableDefinitions["xr"]["units"] = "m";
    
    m_mVariableDefinitions["xl_utm_x"]["description"] = "UTM X coordinate of left riverbank";
    m_mVariableDefinitions["xl_utm_x"]["longname"] = "left bank UTM X";
    m_mVariableDefinitions["xl_utm_x"]["units"] = "m";
    
    m_mVariableDefinitions["xl_utm_y"]["description"] = "UTM Y coordinate of left riverbank";
    m_mVariableDefinitions["xl_utm_y"]["longname"] = "left bank UTM Y";
    m_mVariableDefinitions["xl_utm_y"]["units"] = "m";
    
    m_mVariableDefinitions["xr_utm_x"]["description"] = "UTM X coordinate of right riverbank";
    m_mVariableDefinitions["xr_utm_x"]["longname"] = "right bank UTM X";
    m_mVariableDefinitions["xr_utm_x"]["units"] = "m";
    
    m_mVariableDefinitions["xr_utm_y"]["description"] = "UTM Y coordinate of right riverbank";
    m_mVariableDefinitions["xr_utm_y"]["longname"] = "right bank UTM Y";
    m_mVariableDefinitions["xr_utm_y"]["units"] = "m";
}

/**
 * @brief Destructor (default implementation)
 */
CDataWriter::~CDataWriter() {}


/**
 * @brief Create and define NetCDF output file structure
 * 
 * Execution flow:
 * 1. Create NetCDF4 file with NC_CLOBBER (overwrite if exists)
 * 2. Define dimensions:
 *    - x: Along-channel coordinate (size = nCrossSections)
 *    - time: Temporal dimension (unlimited or fixed)
 * 3. Define variables: 2D arrays [time, x]
 * 4. Add CF-compliant attributes (long_name, units, description)
 * 5. Add global attributes (program name, author, version, run timestamp)
 * 6. End definition mode (nc_enddef)
 * 7. Write coordinate vectors (x, time)
 * 
 * NetCDF4 features used:
 * - Compression: Enabled by default (reduces file size 3-10x)
 * - Chunking: Automatic for unlimited dimension
 * - Parallel I/O: Not used (serial writes)
 * 
 * @param m_pSimulation Pointer to simulation (provides dimensions, metadata)
 * 
 * @note Time dimension behavior:
 * - Unlimited (NC_UNLIMITED): If detail==2 or saveAllTimesteps==true
 * - Fixed size: Otherwise, uses m_vOutputTimes.size()
 * 
 * @warning Error handling: Prints to cerr but continues execution
 * @see nSetOutputData() for writing actual data
 */
void CDataWriter::nDefineNetCDFFile(const CSimulation* m_pSimulation) {
    // Create a new NetCDF4 file (overwrite if exists)
    int status = nc_create(m_pSimulation->m_strOutFile.c_str(), NC_NETCDF4 | NC_CLOBBER, &m_ncId);
    if (status != NC_NOERR) {
        std::cerr << "Error creating NetCDF file: " << nc_strerror(status) << std::endl;
        return;
    }

    //! Define the dimensions: x and t
    //==================================================================================================================
    nc_def_dim(m_ncId, "x", static_cast<size_t>(m_pSimulation->m_nCrossSectionsNumber), &n_XId);
    nc_def_var(m_ncId, "x", NC_DOUBLE, 1, &n_XId, &n_VariableId);

    // Use strlen() to get actual string length
    const char* x_long_name = "along estuary coordinate";
    const char* x_units = "m";
    nc_put_att_text(m_ncId, n_VariableId, "long_name", strlen(x_long_name), x_long_name);
    nc_put_att_text(m_ncId, n_VariableId, "units", strlen(x_units), x_units);

    if (m_pSimulation->m_nLogFileDetail == 2 || m_pSimulation->bGetSaveAllTimesteps()) {
        size_t time_len = NC_UNLIMITED;  // Unlimited time
        nc_def_dim(m_ncId, "time", time_len, &n_TId);
        
        // Define time variable with NC_UNLIMITED dimension
        nc_def_var(m_ncId, "time", NC_DOUBLE, 1, &n_TId, &n_VariableId);
        const char* time_long_name = "time";
        const char* time_units = "s";
        nc_put_att_text(m_ncId, n_VariableId, "long_name", strlen(time_long_name), time_long_name);
        nc_put_att_text(m_ncId, n_VariableId, "units", strlen(time_units), time_units);
    }
    else {
        nc_def_dim(m_ncId, "time", static_cast<size_t>(m_pSimulation->m_vOutputTimes.size()), &n_TId);
        nc_def_var(m_ncId, "time", NC_DOUBLE, 1, &n_TId, &n_VariableId);
        
        // Use strlen() for time variable attributes
        const char* time_long_name = "time";
        const char* time_units = "s";
        nc_put_att_text(m_ncId, n_VariableId, "long_name", strlen(time_long_name), time_long_name);
        nc_put_att_text(m_ncId, n_VariableId, "units", strlen(time_units), time_units);
    }

    // Define the variable with two dimensions (time, x)
    //==================================================================================================================
    n_DimensionsIds[0] = n_TId;
    n_DimensionsIds[1] = n_XId;

    // Limpiar el mapa de IDs de variables antes de llenarlo
    m_mVariableIds.clear();

    for (const auto & m_vOutputVariable : m_pSimulation->m_vOutputVariables) {
        const char *outputVariableName = m_vOutputVariable.c_str();
        int varId;  // Variable local para cada variable

        // Special case: Manning n is constant in time → store as 1D n(x)
        if (m_vOutputVariable == "n") {
            const int x_dim = n_XId;
            status = nc_def_var(m_ncId, outputVariableName, NC_FLOAT, 1, &x_dim, &varId);
            if (status != NC_NOERR) {
                std::cerr << "Error defining variable '" << outputVariableName
                          << "': " << nc_strerror(status) << std::endl;
                continue;
            }

            // Compression for 1D variable is optional but harmless
            status = nc_def_var_deflate(m_ncId, varId, 1, 1, 5);
            if (status != NC_NOERR) {
                std::cerr << "Warning: Could not enable compression for '" << outputVariableName
                          << "': " << nc_strerror(status) << std::endl;
            }

            // Chunking for 1D variable
            size_t chunks[1] = {static_cast<size_t>(m_pSimulation->m_nCrossSectionsNumber)};
            status = nc_def_var_chunking(m_ncId, varId, NC_CHUNKED, chunks);
            if (status != NC_NOERR) {
                std::cerr << "Warning: Could not set chunking for '" << outputVariableName
                          << "': " << nc_strerror(status) << std::endl;
            }

            std::string longname = strGetVariableMetadata(outputVariableName, "longname");
            std::string units = strGetVariableMetadata(outputVariableName, "units");
            std::string description = strGetVariableMetadata(outputVariableName, "description");
            if (!longname.empty()) nc_put_att_text(m_ncId, varId, "long_name", longname.length(), longname.c_str());
            if (!units.empty()) nc_put_att_text(m_ncId, varId, "units", units.length(), units.c_str());
            if (!description.empty()) nc_put_att_text(m_ncId, varId, "description", description.length(), description.c_str());

            m_mVariableIds[outputVariableName] = varId;
            m_nVarIdManningX = varId;
            continue;
        }

        status = nc_def_var(m_ncId, outputVariableName, NC_FLOAT, 2, n_DimensionsIds, &varId);
        if (status != NC_NOERR) {
            std::cerr << "Error defining variable '" << outputVariableName 
                      << "': " << nc_strerror(status) << std::endl;
            continue;
        }

        // Enable compression (deflate level 4 = good balance speed/size)
        // Shuffle filter reorganizes bytes for better compression (especially for floats)
        status = nc_def_var_deflate(m_ncId, varId, 
                                     1,  // Shuffle enabled (improves compression ~30%)
                                     1,  // Deflate compression enabled
                                     5); // Compression level (1=fast, 9=max, 4=balanced)
        if (status != NC_NOERR) {
            std::cerr << "Warning: Could not enable compression for '" << outputVariableName 
                      << "': " << nc_strerror(status) << std::endl;
            // Continue anyway - compression is optional optimization
        }

        // Optimize chunk size for sequential time-series writes
        // Chunk = [1 timestep, all cross-sections] → best for time-series access pattern
        size_t chunks[2] = {1, static_cast<size_t>(m_pSimulation->m_nCrossSectionsNumber)};
        status = nc_def_var_chunking(m_ncId, varId, NC_CHUNKED, chunks);
        if (status != NC_NOERR) {
            std::cerr << "Warning: Could not set chunking for '" << outputVariableName 
                      << "': " << nc_strerror(status) << std::endl;
            // Continue anyway - default chunking will be used
        }

        // Get metadata and use actual length
        std::string longname = strGetVariableMetadata(outputVariableName, "longname");
        std::string units = strGetVariableMetadata(outputVariableName, "units");
        std::string description = strGetVariableMetadata(outputVariableName, "description");

        // Use actual length for each string
        if (!longname.empty()) {
            nc_put_att_text(m_ncId, varId, "long_name", longname.length(), longname.c_str());
        }
        
        if (!units.empty()) {
            nc_put_att_text(m_ncId, varId, "units", units.length(), units.c_str());
        }
        
        if (!description.empty()) {
            nc_put_att_text(m_ncId, varId, "description", description.length(), description.c_str());
        }

        //! Store the specific ID for each variable
        m_mVariableIds[outputVariableName] = varId;
    }

    // Also store static (x-only) profiles for variables that are constant in time.
    // - n(x): Manning roughness profile
    // - Sf(x): lateral storage factor profile
    // These are always defined if their names are not already used by 2D outputs.
    {
        const int x_dim = n_XId;

        if (m_mVariableIds.find("n") == m_mVariableIds.end()) {
            int varId = -1;
            status = nc_def_var(m_ncId, "n", NC_FLOAT, 1, &x_dim, &varId);
            if (status == NC_NOERR) {
                const std::string longname = strGetVariableMetadata("n", "longname");
                const std::string units = strGetVariableMetadata("n", "units");
                const std::string description = strGetVariableMetadata("n", "description");
                if (!longname.empty()) nc_put_att_text(m_ncId, varId, "long_name", longname.length(), longname.c_str());
                if (!units.empty()) nc_put_att_text(m_ncId, varId, "units", units.length(), units.c_str());
                if (!description.empty()) nc_put_att_text(m_ncId, varId, "description", description.length(), description.c_str());
                m_mVariableIds["n"] = varId;
                m_nVarIdManningX = varId;
            } else {
                std::cerr << "Warning: could not define static variable 'n': " << nc_strerror(status) << std::endl;
            }
        }

        if (m_mVariableIds.find("Sf") == m_mVariableIds.end()) {
            int varId = -1;
            status = nc_def_var(m_ncId, "Sf", NC_FLOAT, 1, &x_dim, &varId);
            if (status == NC_NOERR) {
                const std::string longname = strGetVariableMetadata("Sf", "longname");
                const std::string units = strGetVariableMetadata("Sf", "units");
                const std::string description = strGetVariableMetadata("Sf", "description");
                if (!longname.empty()) nc_put_att_text(m_ncId, varId, "long_name", longname.length(), longname.c_str());
                if (!units.empty()) nc_put_att_text(m_ncId, varId, "units", units.length(), units.c_str());
                if (!description.empty()) nc_put_att_text(m_ncId, varId, "description", description.length(), description.c_str());
                m_nVarIdStorageSfX = varId;
            } else {
                std::cerr << "Warning: could not define static variable 'Sf': " << nc_strerror(status) << std::endl;
            }
        }
    }

    //! Include global attributes - use actual length
    const char* program_name = NAME;
    const char* author_name = AUTHOR;
    const char* version_name = VERSION;
    
    nc_put_att_text(m_ncId, NC_GLOBAL, "Program", strlen(program_name), program_name);
    nc_put_att_text(m_ncId, NC_GLOBAL, "Author", strlen(author_name), author_name);
    nc_put_att_text(m_ncId, NC_GLOBAL, "Version", strlen(version_name), version_name);

    //! Include running time
    const time_t t = time(nullptr);
    const tm* m_tSysTime = localtime(&t);

    ostringstream oss;
    oss << put_time(m_tSysTime, "%Y-%m-%d %H:%M:%S");

    const string formatted_time = oss.str();
    nc_put_att_text(m_ncId, NC_GLOBAL, "Running on", formatted_time.length(), formatted_time.c_str());

    // Get start date from CDataReader (must be passed to CSimulation)
    ostringstream time_offset_oss;
    time_offset_oss << setfill('0') << setw(4) << m_pSimulation->nGetSimStartYear() << "-"
                    << setw(2) << m_pSimulation->nGetSimStartMonth() << "-"
                    << setw(2) << m_pSimulation->nGetSimStartDay() << " "
                    << setw(2) << m_pSimulation->nGetSimStartHour() << ":"
                    << setw(2) << m_pSimulation->nGetSimStartMin() << ":"
                    << setw(2) << m_pSimulation->nGetSimStartSec();

    const string time_offset = time_offset_oss.str();
    nc_put_att_text(m_ncId, NC_GLOBAL, "time_offset", time_offset.length(), time_offset.c_str());

    //! Time_offset description
    const char* time_offset_description = "Simulation start time used as time offset for all time variables";
    nc_put_att_text(m_ncId, NC_GLOBAL, "time_offset_description", strlen(time_offset_description), time_offset_description);

    //! Relative time units
    const char* time_units = "seconds since simulation start";
    nc_put_att_text(m_ncId, NC_GLOBAL, "time_units", strlen(time_units), time_units);
    
    //! YAML configuration as separate dimension + variable (preserves formatting better)
    if (!m_pSimulation->m_strYAMLConfigContent.empty()) {
        // Create dimension for config string length
        int config_len_dim;
        size_t config_length = m_pSimulation->m_strYAMLConfigContent.length();
        nc_def_dim(m_ncId, "config_length", config_length, &config_len_dim);
        
        // Create char variable to store config
        int config_var;
        nc_def_var(m_ncId, "yaml_configuration", NC_CHAR, 1, &config_len_dim, &config_var);
        nc_put_att_text(m_ncId, config_var, "long_name", 
                       strlen("YAML configuration file content"), 
                       "YAML configuration file content");
        
        // Store config in variable (done after nc_enddef below)
    }

    //! End the definition of NetCDF file
    status = nc_enddef(m_ncId);
    if (status != NC_NOERR) {
        std::cerr << "Error ending NetCDF definition: " << nc_strerror(status) << std::endl;
        return;
    }

    //! Write coordinates and times
    nc_put_var_double(m_ncId, n_XId, m_pSimulation->m_vCrossSectionX.data());
    //! Only write all times if NOT using unlimited mode (detail != 2 and not saveAllTimesteps)
    if (m_pSimulation->m_nLogFileDetail != 2 && !m_pSimulation->bGetSaveAllTimesteps()) {
        nc_put_var_double(m_ncId, n_TId, m_pSimulation->m_vOutputTimes.data());
    }
    
    //! Write YAML configuration content (after enddef)
    if (!m_pSimulation->m_strYAMLConfigContent.empty()) {
        int config_var;
        if (nc_inq_varid(m_ncId, "yaml_configuration", &config_var) == NC_NOERR) {
            nc_put_var_text(m_ncId, config_var, m_pSimulation->m_strYAMLConfigContent.c_str());
        }
    }
}


/**
 * @brief Retrieve metadata for a specific variable
 * 
 * Accessor for m_mVariableDefinitions map.
 * 
 * @param strVariable Variable name (e.g., "A", "Q", "S")
 * @param strField Metadata field ("longname", "units", "description")
 * @return Metadata string, or empty string if not found
 * 
 * @note No error checking - returns empty string for missing entries
 */
string CDataWriter::strGetVariableMetadata(const string& strVariable, const string& strField) {
    return m_mVariableDefinitions[strVariable][strField];
}
/**
 * @brief Write simulation data to NetCDF file at current timestep
 * 
 * Procedure:
 * 1. Determine time index based on save mode:
 *    - All timesteps (detail==2 or saveAllTimesteps): Use m_nTimeLogId, increment after write
 *    - Selected times: Use m_nTimeId (pre-computed index)
 * 2. Write current_time value if using unlimited dimension
 * 3. Loop over output variables:
 *    - Get data from simulation via vGetVariable()
 *    - Write 1D slice [time=current, x=0:N] to NetCDF
 * 4. Report errors to cerr
 * 
 * NetCDF write call:
 *   nc_put_vara_double(ncId, varId, start={time_idx, 0}, count={1, nCrossSections}, data)
 * 
 * @param m_pSimulation Pointer to simulation (provides data)
 * 
 * @note Performance:
 * - Each variable written separately (not efficient for many variables)
 * - No buffering (writes directly to disk)
 * - For large outputs, consider chunking/buffering strategies
 * 
 * @warning Does not sync/flush - relies on automatic NetCDF buffering
 * @see nCloseNetCDFFile() for final sync
 */
void CDataWriter::nSetOutputData(CSimulation *m_pSimulation) const {

    int status = 0;

    // Write static (x-only) variables once.
    if (!m_bWroteStaticX) {
        const int N = m_pSimulation->m_nCrossSectionsNumber;
        if (N > 0) {
            if (m_nVarIdManningX >= 0) {
                vector<float> n_x(static_cast<size_t>(N), 0.0f);
                // Prefer geometry-defined Manning profile if available.
                if (static_cast<int>(m_pSimulation->estuary.size()) >= N) {
                    for (int i = 0; i < N; ++i) {
                        n_x[static_cast<size_t>(i)] = static_cast<float>(m_pSimulation->estuary[i].dGetManningNumber());
                    }
                } else if (m_pSimulation->m_vCrossSectionManningNumber.size() == static_cast<size_t>(N)) {
                    for (int i = 0; i < N; ++i) {
                        n_x[static_cast<size_t>(i)] = static_cast<float>(m_pSimulation->m_vCrossSectionManningNumber[static_cast<size_t>(i)]);
                    }
                }
                status = nc_put_var_float(m_ncId, m_nVarIdManningX, n_x.data());
                if (status != NC_NOERR) {
                    std::cerr << "Error writing static variable 'n': " << nc_strerror(status) << std::endl;
                }
            }

            if (m_nVarIdStorageSfX >= 0) {
                vector<float> sf_x(static_cast<size_t>(N), 1.0f);
                for (int i = 0; i < N; ++i) {
                    sf_x[static_cast<size_t>(i)] = static_cast<float>(m_pSimulation->dGetLateralStorageFactor(i));
                }
                status = nc_put_var_float(m_ncId, m_nVarIdStorageSfX, sf_x.data());
                if (status != NC_NOERR) {
                    std::cerr << "Error writing static variable 'Sf': " << nc_strerror(status) << std::endl;
                }
            }
        }
        m_bWroteStaticX = true;
    }

    //! Choose to save all timestep if m_nLogFileDetail == 2 or only when indicated
    size_t start[2] = {static_cast<size_t>(m_pSimulation->m_nTimeId), 0};
    if (m_pSimulation->m_nLogFileDetail == 2 || m_pSimulation->bGetSaveAllTimesteps()) {
        start[0] = static_cast<size_t>(m_pSimulation->m_nTimeLogId);
        
        // Write current time value to unlimited time dimension
        size_t time_start = static_cast<size_t>(m_pSimulation->m_nTimeLogId);
        size_t time_count = 1;
        double current_time = m_pSimulation->m_dCurrentTime;
        status = nc_put_vara_double(m_ncId, n_TId, &time_start, &time_count, &current_time);
        if (status != NC_NOERR) {
            std::cerr << "Error writing time variable: " << nc_strerror(status) << std::endl;
        }
        
        m_pSimulation->m_nTimeLogId++;
    }

    //! Write variable data
    for (const auto& strOutputVariable : m_pSimulation->m_vOutputVariables) {
        // n is stored as 1D n(x) (static) and already written above.
        if (strOutputVariable == "n" && m_nVarIdManningX >= 0) {
            continue;
        }

        if (m_mVariableIds.find(strOutputVariable) == m_mVariableIds.end()) {
            std::cerr << "Variable ID for '" << strOutputVariable << "' not found!" << std::endl;
            continue;
        }

        vector<double> data = m_pSimulation->vGetVariable(strOutputVariable);

        if (data.empty()) {
            std::cerr << "No data available for variable '" << strOutputVariable << "'." << std::endl;
            continue;
        }

        // ✅ FIX: Use actual data size, not m_nCrossSectionsNumber
        // Some variables (e.g., D1Factor, D2Factor) have size n+1 for interfaces
        const size_t count[2] = {1, data.size()};

        // Convert double to float for 50% size reduction
        vector<float> data_float(data.size());
        for (size_t i = 0; i < data.size(); ++i) {
            data_float[i] = static_cast<float>(data[i]);
        }

        status = nc_put_vara_float(
            m_ncId,
            m_mVariableIds.at(strOutputVariable),
            start,
            count,
            data_float.data()
        );

        if (status != NC_NOERR) {
            std::cerr << "Error writing variable '" << strOutputVariable
                << "': " << nc_strerror(status) << std::endl;
            break;
        }
    }



    if (status != 0) {
        cerr << "Error while creating output variables" << endl;
    }
}



/**
 * @brief Close NetCDF output file and flush buffers
 * 
 * Operations:
 * 1. Check if file is open (m_ncId >= 0)
 * 2. Call nc_close() to:
 *    - Flush all pending writes
 *    - Write final metadata
 *    - Close file descriptor
 * 3. Set m_ncId = -1 to prevent double-close
 * 
 * @note nc_close() is idempotent-safe due to check
 * @warning File may be corrupted if not closed properly (power failure, kill -9)
 * 
 * Best practices:
 * - Always call in destructor or finally block
 * - Check disk space before long simulations
 * - Use nc_sync() periodically for crash recovery
 */
void CDataWriter::nCloseNetCDFFile()
{
    // Verify file is open before closing
    if (m_ncId < 0) {
        std::cout << "NetCDF file already closed or never opened" << std::endl;
        return;
    }
        
    int status = nc_close(m_ncId);
    if (status != NC_NOERR) {
        std::cerr << "Error closing NetCDF file: " << nc_strerror(status) << std::endl;
    }
    // Mark as closed regardless of result
    m_ncId = -1;
}
