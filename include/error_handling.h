/*!
*
 * \brief
 * \details TODO 001 This is a more detailed description of the Error Handling
 * \author Manuel Cobos Budia

 * \date 2024
 * \copyright GNU General Public License
 *
 * \file error_handling.h
 * \brief Contains ErrorHandling definitions
 *
 */
#ifndef ERROR_HANDLING_H
#define ERROR_HANDLING_H

#include <string>
using std::string;


// Return codes
static constexpr int RTN_OK = 0;
static constexpr int RTN_ERR_INI = 1;
static constexpr int RTN_ERR_SV_DIR = 2;
static constexpr int RTN_ERR_BADLY_FORMAT_COLON = 3;
static constexpr int RTN_CLOSE_NETCDF = 4;
static constexpr int RTN_EXTRA_LINES = 5;
static constexpr int RTN_ERR_LOGFILE = 10;
static constexpr int RTN_ERR_OUTFILE = 11;
static constexpr int RTN_ERR_TIMEUNITS = 36;

static constexpr int ERR = 0;

string strGetErrorText(int);


#endif // ERROR_HANDLING_H
