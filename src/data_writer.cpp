/*!
*
 * \file data_writer.cpp
 * \brief Write output NetCDF file
 * \details TODO 001 A more detailed description of these routines.
 * \author Manuel Cobos Budia

 * \date 20240905
 * \copyright GNU General Public License
 *
 */

/*==============================================================================================================================

==============================================================================================================================*/
#include <cstdlib>

#include <sstream>
using std::stringstream;

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

#include "simulation.h"
#include "data_writer.h"
#include "main.h"


//======================================================================================================================
//! Class constructor
//======================================================================================================================
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

    m_mVariableDefinitions["I1"]["description"] = "Term I1 of the balance"; //!TODO 021: include a best description
    m_mVariableDefinitions["I1"]["longname"] = "I1";
    m_mVariableDefinitions["I1"]["units"] = "";

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

    m_mVariableDefinitions["Qb"]["description"] = "Distance from the thalweg to the right riverbank";
    m_mVariableDefinitions["Qb"]["longname"] = "x-right";
    m_mVariableDefinitions["Qb"]["units"] = "m";
}

CDataWriter::~CDataWriter() = default;



//======================================================================================================================
//! The definition of NetCDF file
//======================================================================================================================
void CDataWriter::nDefineNetCDFFile(const CSimulation* m_pSimulation) {
    //! Create a new NetCDF file
    int status = nc_create(m_pSimulation->m_strOutFile.c_str(), NC_NETCDF4 | NC_CLOBBER, &m_ncId);

    //! Define the dimensions: x and t
    //==================================================================================================================
    nc_def_dim(m_ncId, "x", static_cast<size_t>(m_pSimulation->m_nCrossSectionsNumber), &n_XId);
    nc_def_var(m_ncId, "x", NC_DOUBLE, 1, &n_XId, &n_VariableId);

    // Long name
    nc_put_att_text(m_ncId, n_VariableId, "long_name", 30, "along estuary coordinate");
    // Units
    nc_put_att_text(m_ncId, n_VariableId, "units", 15, "m");

    if (m_pSimulation->m_nLogFileDetail == 2) {
        size_t time_len = NC_UNLIMITED;  // Unlimited time
        nc_def_dim(m_ncId, "time", time_len, &n_TId);
    }
    else {
        nc_def_dim(m_ncId, "time", static_cast<size_t>(m_pSimulation->m_vOutputTimes.size()), &n_TId);
        nc_def_var(m_ncId, "time", NC_DOUBLE, 1,&n_TId, &n_VariableId);
    }
    // Long name
    nc_put_att_text(m_ncId, n_VariableId, "long_name", 30, "time");
    // Units
    nc_put_att_text(m_ncId, n_VariableId, "units", 15, "s");

    // Define the variable with two dimensions (time, x)
    //==================================================================================================================
    n_DimensionsIds[0] = n_TId;
    n_DimensionsIds[1] = n_XId;

    for (const auto & m_vOutputVariable : m_pSimulation->m_vOutputVariables) {
        const char *outputVariableName = m_vOutputVariable.c_str();
        nc_def_var(m_ncId, outputVariableName, NC_DOUBLE, 2, n_DimensionsIds, &n_VariableId);

        // Long name
        nc_put_att_text(m_ncId, n_VariableId, "long_name", 30, strGetVariableMetadata(outputVariableName, "longname").c_str());

        // Units
        nc_put_att_text(m_ncId, n_VariableId, "units", 15, strGetVariableMetadata(outputVariableName, "units").c_str());

        // Description
        nc_put_att_text(m_ncId, n_VariableId, "description", 50, strGetVariableMetadata(outputVariableName, "description").c_str());

        //! Include the variable Id to the dictionary
        m_mVariableIds[outputVariableName] = n_VariableId;
    }

    //! Include global attributes
    nc_put_att_text(m_ncId, NC_GLOBAL, "Program", 40, NAME);
    nc_put_att_text(m_ncId, NC_GLOBAL, "Author", 20, AUTHOR);
    nc_put_att_text(m_ncId, NC_GLOBAL, "Version", 17, VERSION);

    //! Include running time
    const time_t t = time(nullptr);           // Get the current time
    const tm* m_tSysTime = localtime(&t);

    // Capture the output of time
    ostringstream oss;
    oss << put_time(m_tSysTime, "%Y-%m-%d %H:%M:%S"); // Formatear la fecha y hora

    // Give format to time
    const string formatted_time = oss.str();
    nc_put_att_text(m_ncId, NC_GLOBAL, "Running on", 40, formatted_time.c_str());

    // End the definition of NetCDF file
    nc_enddef(m_ncId);

    //! Write coordenates and times. NetCDF only allow to write data outside the definition mode
    nc_put_var_double(m_ncId,  n_XId, m_pSimulation->m_vCrossSectionX.data());
    if (m_pSimulation->m_nLogFileDetail != 2) {
        nc_put_var_double(m_ncId,  n_TId, m_pSimulation->m_vOutputTimes.data());
    }
}


//======================================================================================================================
//! Get variable metadata
//======================================================================================================================
string CDataWriter::strGetVariableMetadata(const string& strVariable, const string& strField) {
    return m_mVariableDefinitions[strVariable][strField];
}
//======================================================================================================================
//! Set data of variables at a given timestep to the NetCDF file
//======================================================================================================================
void CDataWriter::nSetOutputData(CSimulation *m_pSimulation) const {

    int status = 0;
    //! Choose to save all timestep if m_nLogFileDetail == 2 or only when indicated
    size_t start[2] = {static_cast<size_t>(m_pSimulation->m_nTimeId), 0};
    if (m_pSimulation->m_nLogFileDetail == 2) {
        start[0] = static_cast<size_t>(m_pSimulation->m_nTimeLogId);
        m_pSimulation->m_nTimeLogId++;
    }

    const size_t count[2] = {1, static_cast<size_t>(m_pSimulation->m_nCrossSectionsNumber)};

    //! Write variable data
    for (const auto & strOutputVariable : m_pSimulation->m_vOutputVariables) {
        vector<double> data = m_pSimulation->vGetVariable(strOutputVariable);
        status = nc_put_vara_double(m_ncId, m_mVariableIds.at(strOutputVariable),  start, count, data.data());
    }

    if (status != 0) {
        cerr << "Error while creating output variables" << endl;
    }
}



//======================================================================================================================
//! Close the NetCDF file
//======================================================================================================================
void CDataWriter::nCloseNetCDFFile(CSimulation *m_pSimulation) {
    int status = nc_close(m_ncId);
    if (status != NC_NOERR) {
        std::cerr << "Error closing the NetCDF file: " << nc_strerror(status) << std::endl;
    }
}
