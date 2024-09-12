
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
        // case RTN_USER_ABORT:
        //    strErr = "aborted by user";
        //    break;
        // case RTN_ERR_BADPARAM:
        //    strErr = "error in command-line parameter";
        //    break;
        case RTN_ERR_INI:
            strErr = "error reading initialization file";
            break;
        case RTN_ERR_SV_DIR:
            strErr = "error in directory name";
            break;
        case RTN_ERR_BADLY_FORMAT_COLON:
            strErr = "badly formatted string colon found";
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


