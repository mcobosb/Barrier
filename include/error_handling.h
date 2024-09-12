/*!
*
 * \brief
 * \details TODO 001 This is a more detailed description of the Error Handling
 * \author Manuel Cobos Budia

 * \date 2024
 * \copyright GNU General Public License
 *
 * \file data_reader.h
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
static constexpr int RTN_ERR_LOGFILE = 10;
static constexpr int RTN_ERR_OUTFILE = 11;
static constexpr int RTN_ERR_TIMEUNITS = 36;

string const ERROR_NOTICE = " with error code ";
string const ERR = "*** ERROR ";

string strGetErrorText(int);


#endif // ERROR_HANDLING_H
