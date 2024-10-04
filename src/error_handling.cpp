
#include <sstream>
using std::ostringstream;

#include <string>
using std::to_string;

#include "error_handling.h"

//===============================================================================================================================
//! Returns an error message given an error code
//===============================================================================================================================
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


