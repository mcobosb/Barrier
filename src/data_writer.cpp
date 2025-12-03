/*!
*
 * \file data_writer.cpp
 * \brief Write output NetCDF file with simulation results
 * \details Handles creation and writing of NetCDF output files containing
 *          temporal evolution of hydraulic variables (area, flow, velocity, etc.)
 *          along the estuarine channel.
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

    m_mVariableDefinitions["I1"]["description"] = "Integral I1 = int(B dh) for source term balance";
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

//======================================================================================================================
//! Class destructor
//======================================================================================================================
CDataWriter::~CDataWriter() {}


//======================================================================================================================
//! The definition of NetCDF file
//======================================================================================================================
void CDataWriter::nDefineNetCDFFile(const CSimulation* m_pSimulation) {
    //! Create a new NetCDF file
    int status = nc_create(m_pSimulation->m_strOutFile.c_str(), NC_NETCDF4 | NC_CLOBBER, &m_ncId);
    if (status != NC_NOERR) {
        std::cerr << "Error creating NetCDF file: " << nc_strerror(status) << std::endl;
        return;
    }

    //! Define the dimensions: x and t
    //==================================================================================================================
    nc_def_dim(m_ncId, "x", static_cast<size_t>(m_pSimulation->m_nCrossSectionsNumber), &n_XId);
    nc_def_var(m_ncId, "x", NC_DOUBLE, 1, &n_XId, &n_VariableId);

    // Usar strlen() para obtener la longitud real de las cadenas
    const char* x_long_name = "along estuary coordinate";
    const char* x_units = "m";
    nc_put_att_text(m_ncId, n_VariableId, "long_name", strlen(x_long_name), x_long_name);
    nc_put_att_text(m_ncId, n_VariableId, "units", strlen(x_units), x_units);

    if (m_pSimulation->m_nLogFileDetail == 2 || m_pSimulation->bGetSaveAllTimesteps()) {
        size_t time_len = NC_UNLIMITED;  // Unlimited time
        nc_def_dim(m_ncId, "time", time_len, &n_TId);
        
        // También definir la variable tiempo con NC_UNLIMITED
        nc_def_var(m_ncId, "time", NC_DOUBLE, 1, &n_TId, &n_VariableId);
        const char* time_long_name = "time";
        const char* time_units = "s";
        nc_put_att_text(m_ncId, n_VariableId, "long_name", strlen(time_long_name), time_long_name);
        nc_put_att_text(m_ncId, n_VariableId, "units", strlen(time_units), time_units);
    }
    else {
        nc_def_dim(m_ncId, "time", static_cast<size_t>(m_pSimulation->m_vOutputTimes.size()), &n_TId);
        nc_def_var(m_ncId, "time", NC_DOUBLE, 1, &n_TId, &n_VariableId);
        
        // Usar strlen() para la variable tiempo también
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
        
        status = nc_def_var(m_ncId, outputVariableName, NC_DOUBLE, 2, n_DimensionsIds, &varId);
        if (status != NC_NOERR) {
            std::cerr << "Error defining variable '" << outputVariableName 
                      << "': " << nc_strerror(status) << std::endl;
            continue;
        }

        // Obtener metadatos y usar su longitud real
        std::string longname = strGetVariableMetadata(outputVariableName, "longname");
        std::string units = strGetVariableMetadata(outputVariableName, "units");
        std::string description = strGetVariableMetadata(outputVariableName, "description");

        // Usar la longitud real de cada cadena
        if (!longname.empty()) {
            nc_put_att_text(m_ncId, varId, "long_name", longname.length(), longname.c_str());
        }
        
        if (!units.empty()) {
            nc_put_att_text(m_ncId, varId, "units", units.length(), units.c_str());
        }
        
        if (!description.empty()) {
            nc_put_att_text(m_ncId, varId, "description", description.length(), description.c_str());
        }

        //! Guardar el ID específico de cada variable
        m_mVariableIds[outputVariableName] = varId;
    }

    //! Include global attributes - usar longitud real
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

    // Obtener fecha de inicio de CDataReader (debe pasarse a CSimulation)
    ostringstream time_offset_oss;
    time_offset_oss << setfill('0') << setw(4) << m_pSimulation->nGetSimStartYear() << "-"
                    << setw(2) << m_pSimulation->nGetSimStartMonth() << "-"
                    << setw(2) << m_pSimulation->nGetSimStartDay() << " "
                    << setw(2) << m_pSimulation->nGetSimStartHour() << ":"
                    << setw(2) << m_pSimulation->nGetSimStartMin() << ":"
                    << setw(2) << m_pSimulation->nGetSimStartSec();

    const string time_offset = time_offset_oss.str();
    nc_put_att_text(m_ncId, NC_GLOBAL, "time_offset", time_offset.length(), time_offset.c_str());

    //! ✅ AÑADIR: Descripción del time_offset
    const char* time_offset_description = "Simulation start time used as time offset for all time variables";
    nc_put_att_text(m_ncId, NC_GLOBAL, "time_offset_description", strlen(time_offset_description), time_offset_description);

    //! ✅ AÑADIR: Unidades de tiempo relativo
    const char* time_units = "seconds since simulation start";
    nc_put_att_text(m_ncId, NC_GLOBAL, "time_units", strlen(time_units), time_units);

    // End the definition of NetCDF file
    status = nc_enddef(m_ncId);
    if (status != NC_NOERR) {
        std::cerr << "Error ending NetCDF definition: " << nc_strerror(status) << std::endl;
        return;
    }

    //! Write coordinates and times
    nc_put_var_double(m_ncId, n_XId, m_pSimulation->m_vCrossSectionX.data());
    // Solo escribir todos los tiempos si NO se usa modo unlimited (detail != 2 y no saveAllTimesteps)
    if (m_pSimulation->m_nLogFileDetail != 2 && !m_pSimulation->bGetSaveAllTimesteps()) {
        nc_put_var_double(m_ncId, n_TId, m_pSimulation->m_vOutputTimes.data());
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
    if (m_pSimulation->m_nLogFileDetail == 2 || m_pSimulation->bGetSaveAllTimesteps()) {
        start[0] = static_cast<size_t>(m_pSimulation->m_nTimeLogId);
        
        // Escribir el valor del tiempo actual en la dimensión unlimited
        size_t time_start = static_cast<size_t>(m_pSimulation->m_nTimeLogId);
        size_t time_count = 1;
        double current_time = m_pSimulation->m_dCurrentTime;
        status = nc_put_vara_double(m_ncId, n_TId, &time_start, &time_count, &current_time);
        if (status != NC_NOERR) {
            std::cerr << "Error writing time variable: " << nc_strerror(status) << std::endl;
        }
        
        m_pSimulation->m_nTimeLogId++;
    }

    const size_t count[2] = {1, static_cast<size_t>(m_pSimulation->m_nCrossSectionsNumber)};

    //! Write variable data
    for (const auto& strOutputVariable : m_pSimulation->m_vOutputVariables) {
        if (m_mVariableIds.find(strOutputVariable) == m_mVariableIds.end()) {
            std::cerr << "Variable ID for '" << strOutputVariable << "' not found!" << std::endl;
            continue;
        }

        vector<double> data = m_pSimulation->vGetVariable(strOutputVariable);

        if (data.empty()) {
            std::cerr << "No data available for variable '" << strOutputVariable << "'." << std::endl;
            continue;
        }

        status = nc_put_vara_double(
            m_ncId,
            m_mVariableIds.at(strOutputVariable),
            start,
            count,
            data.data()
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



//======================================================================================================================
//! Close the NetCDF file
//======================================================================================================================
void CDataWriter::nCloseNetCDFFile()
{
    // Verificar que el archivo esté abierto antes de cerrarlo
    if (m_ncId < 0) {
        std::cout << "NetCDF file already closed or never opened" << std::endl;
        return;
    }
        
    int status = nc_close(m_ncId);
    if (status != NC_NOERR) {
        std::cerr << "Error closing NetCDF file: " << nc_strerror(status) << std::endl;
    }
    // Marcar como cerrado independientemente del resultado
    m_ncId = -1;
}
