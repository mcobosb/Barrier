
#include <sstream>
using std::ostringstream;

#include <string>
using std::to_string;

#include "error_handling.h"

/**
 * @brief Convert error code to human-readable message
 * 
 * Error codes:
 * - RTN_OK (0): Successful completion
 * - RTN_ERR_INI: Configuration file read error (missing file, syntax error)
 * - RTN_ERR_SV_DIR: Invalid output directory path
 * - RTN_ERR_BADLY_FORMAT_COLON: Malformed string with colon separator
 * - RTN_CLOSE_NETCDF: NetCDF file close error (disk full, permissions)
 * - RTN_EXTRA_LINES: Unexpected lines in config file
 * - RTN_ERR_TIMEUNITS: Invalid time units in forcing data
 * 
 * @param nErr Error code from error_handling.h
 * @return Error description string
 * 
 * @note Used by:
 * - CSimulation::bDoSimulationEnd() for final status
 * - Exception handlers for debugging
 * - Log files for post-processing analysis
 */
string strGetErrorText(int const nErr)
{
    string strErr;

    switch (nErr)
    {
        case RTN_OK:
            strErr = "Run successfully end";
            break;
        case RTN_ERR_INI:
            strErr = "Error reading configuration file";
            break;
        case RTN_ERR_SV_DIR:
            strErr = "Error in directory name";
            break;
        case RTN_ERR_BADLY_FORMAT_COLON:
            strErr = "Badly formatted string colon found";
            break;
        case RTN_CLOSE_NETCDF:
            strErr = "Error closing the NetCDF file";
            break;
        case RTN_EXTRA_LINES:
            strErr = "Remove extra lines in configuration file";
            break;
        case RTN_ERR_TIMEUNITS:
            strErr = "error in time units";
            break;
        default:
            // should never get here
            strErr = "unknown error";
            break;
    }

    return strErr;
}


