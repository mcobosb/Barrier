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
#include <cstdlib> // for strtol() and strtod()
#include <data_reader.h>
#include <fstream>
using std::ifstream;

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

#include <iomanip>
using std::put_time;

#include "simulation.h"
#include "data_writer.h"
#include "main.h"
#include "error_handling.h"
#include "utils.h"

// Error handling
#define NC_CHECK(status) if (status != NC_NOERR) { std::cerr << "Error: " << nc_strerror(status) << std::endl; return 1; }

//======================================================================================================================
//! Class constructor
//======================================================================================================================
CDataWriter::CDataWriter() {

    // Creating the attributes of the main variables
    m_mVariableDefinitions["A"]["description"] = "Cross-sectional flooding area";
    m_mVariableDefinitions["A"]["longname"] = "Area";
    m_mVariableDefinitions["A"]["units"] = "m2";

    m_mVariableDefinitions["Q"]["description"] = "Cross-sectional averaged water discharge";
    m_mVariableDefinitions["Q"]["longname"] = "water discharge";
    m_mVariableDefinitions["Q"]["units"] = "m3/s";

    m_mVariableDefinitions["q"]["description"] = "Fluvial contributions along the estuary";
    m_mVariableDefinitions["q"]["longname"] = "water flux contributions";
    m_mVariableDefinitions["q"]["units"] = "m3/s";

    m_mVariableDefinitions["Rh"]["description"] = "Hydraulic radius of the flooding area";
    m_mVariableDefinitions["Rh"]["longname"] = "hydraulic radius";
    m_mVariableDefinitions["Rh"]["units"] = "m";

    m_mVariableDefinitions["I1"]["description"] = "Term I1 of the balance"; //!TODO 021: include a best description
    m_mVariableDefinitions["I1"]["longname"] = "I1";
    m_mVariableDefinitions["I1"]["units"] = "";

    m_mVariableDefinitions["eta"]["description"] = "Free surface elevation";
    m_mVariableDefinitions["eta"]["longname"] = "flooding water depth";
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

    m_mVariableDefinitions["Qb"]["description"] = "Bedload sediment transport";
    m_mVariableDefinitions["Qb"]["longname"] = "bedload sediment transport";
    m_mVariableDefinitions["Qb"]["units"] = "m3/s";

    m_mVariableDefinitions["Qs"]["description"] = "Suspended sediment transport";
    m_mVariableDefinitions["Qs"]["longname"] = "suspended sediment transport";
    m_mVariableDefinitions["Qs"]["units"] = "m3/s";

    m_mVariableDefinitions["Qt"]["description"] = "Total sediment transport";
    m_mVariableDefinitions["Qt"]["longname"] = "total sediment transport";
    m_mVariableDefinitions["Qt"]["units"] = "m3/s";

    m_mVariableDefinitions["rho"]["description"] = "Fluid density (salinity, bedload and suspended sediment)";
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
int CDataWriter::nDefineNetCDFFile(const CSimulation* m_pSimulation) {
    //! Create a new NetCDF file
    int status = nc_create(m_pSimulation->m_strOutFile.c_str(), NC_NETCDF4 | NC_CLOBBER, &m_ncId);
    NC_CHECK(status);

    //! Define the dimensions: x and t
    status = nc_def_dim(m_ncId, "x", static_cast<size_t>(m_pSimulation->m_nCrossSectionsNumber), &n_XId);
    NC_CHECK(status);
    size_t time_len = NC_UNLIMITED;  // Unlimited time
    status = nc_def_dim(m_ncId, "time", time_len, &n_TId);
    NC_CHECK(status);

    // Define the variable with two dimensions (time, x)
    n_DimensionsIds[0] = n_TId;
    n_DimensionsIds[1] = n_XId;

    for (const auto & m_vOutputVariable : m_pSimulation->m_vOutputVariables) {
        const char *outputVariableName = m_vOutputVariable.c_str();
        status = nc_def_var(m_ncId, outputVariableName, NC_FLOAT, 2, n_DimensionsIds, &n_VariableId);

        // Long name
        nc_put_att_text(m_ncId, n_VariableId, "long_name", 30, strGetVariableMetadata(outputVariableName, "longname").c_str());

        // Units
        nc_put_att_text(m_ncId, n_VariableId, "units", 15, strGetVariableMetadata(outputVariableName, "units").c_str());

        // Description
        nc_put_att_text(m_ncId, n_VariableId, "description", 50, strGetVariableMetadata(outputVariableName, "description").c_str());

        //! Include the variable Id to the dictionary
        m_mVariableIds[outputVariableName] = n_VariableId;
        NC_CHECK(status);
    }

    //! Include global attributes
    nc_put_att_text(m_ncId, NC_GLOBAL, "Program", 40, NAME);
    nc_put_att_text(m_ncId, NC_GLOBAL, "Author", 20, AUTHOR);
    nc_put_att_text(m_ncId, NC_GLOBAL, "Version", 17, VERSION);

    // time_t m_tSysTime = time(nullptr);
    // nc_put_att_text(m_ncId, NC_GLOBAL, "Running on", 20, put_time(localtime(&m_tSysTime), "%H:%M on %A %d %B %Y"));

    // End the definition of NetCDF file
    status = nc_enddef(m_ncId);
    NC_CHECK(status);
    return status;
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
int CDataWriter::nSetOutputData(CSimulation *m_pSimulation) const {

    size_t start[2] = {static_cast<size_t>(m_pSimulation->m_nTimeId), 0};
    size_t count[2] = {1, static_cast<size_t>(m_pSimulation->m_nCrossSectionsNumber)};

    //! Write variable data
    for (const auto & strOutputVariable : m_pSimulation->m_vOutputVariables) {
        vector<double> data = m_pSimulation->vGetVariable(strOutputVariable);
        int status = nc_put_vara_double(m_ncId, m_mVariableIds.at(strOutputVariable),  start, count, data.data());
        NC_CHECK(status);
    }
    return false;
    // std::cout << "Datos del tiempo " << t << " escritos correctamente.\n";
}
