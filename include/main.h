/*!
*
 * Main function
 * \details
 * \author Manuel Cobos Budia

 * \date 2024
 * \copyright GNU General Public License
 *
 * \file main.h
 * \brief
 *
 */

//======================================================================================================================
// TODO 001: A more detail description of main function.
// TODO 002: Improve error handling functionalities.
// TODO 003: Add computing progress function.
// TODO 004: Save Predicted and Corrected variables for m_nLogDetail == 1 and 2.
// TODO 005: Include hydrographs input. Add the interpolated value at node i in term Gv[0].
// TODO 006: Include density, salinity and sediment load
// TODO 007: Options for another time units

//======================================================================================================================


#ifndef MAIN_H
#define MAIN_H

#include <string>
using std::string;


#include <sstream>
using std::ostream;
using std::ostringstream;



// Global constant
//! gravitational constant
inline double G = 9.81;
//! minimum cross-section area
inline double DRY_AREA = 1;
//! minimum cross-section water flow (to ensure that U = Q/A ~ 0)
inline double DRY_Q = 1e-2;

// Intel x86, byte order is little-endian, 32-bit
string const PLATFORM = "Intel x86/GNU C++";



inline auto NAME = "Saint-Venant Equations Solver";
inline auto AUTHOR =  "Manuel Cobos (GDFA, University of Granada)";
inline auto VERSION = "v0.3.0 - 20240912";
string const PROGRAM_NAME = "Saint Venant Equations Solver version 0.3.0 (12 Sep 2024)";
string const PROGRAM_NAME_SHORT = "SV";


string const NOTE = "      Note ";
string const COPYRIGHT = "(C) 2024 Manuel Cobos";
string const LINE = "-------------------------------------------------------------------------------";
string const DISCLAIMER1 = "This program is distributed in the hope that it will be useful. but WITHOUT ANY";
string const DISCLAIMER2 = "WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A";
string const DISCLAIMER3 = "PARTICULAR PURPOSE. See the GNU General Public License for more details. You";
string const DISCLAIMER4 = "should have received a copy of the GNU General Public License along with this";
string const DISCLAIMER5 = "program; if not. contact the Free Software Foundation. Inc.. 675 Mass Ave.";
string const DISCLAIMER6 = "Cambridge. MA 02139. USA.";

string const START_NOTICE = "    - Started on ";
string const INITIALIZING_NOTICE = "    - Initializing";
string const SIMULATING = "\r  - Simulating ";
string const RUN_END_NOTICE = "    - Run ended at ";

// Not likely that user will need to change these
static constexpr int BUF_SIZE = 2048;                                     // Max length (inc. terminating NULL) of any C-type string

// clock_t is a signed long: see <time.h>
constexpr long CLOCK_T_MIN = LONG_MIN;
constexpr double CLOCK_T_RANGE = static_cast<double>(LONG_MAX) - static_cast<double>(CLOCK_T_MIN);

inline float CDX_MURILLO = 0.6;

#endif // MAIN_H