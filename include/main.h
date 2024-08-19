/*   TODO 031 Get raster slice output working with multiple slices

*/

#ifndef MAIN_H
#define MAIN_H
/*===============================================================================================================================
*/
#include "simulation.h"

#include <climits>

#include <sstream>
using std::ostream;
using std::ostringstream;
//================================


//===================================================== hard-wired constants ====================================================
char const COLON = ':';
char const COMMA = ',';
char const DASH = '-';
char const PATH_SEPARATOR = '/'; // Works for Windows too!
char const QUOTE1 = ';';
char const QUOTE2 = '#';
char const SLASH = '/';
char const SPACE = ' ';
char const TILDE = '~';


string const VERSION = "0.0.1";
// Log file detail level
int const NO_LOG_FILE = 0;
int const LOG_FILE_HIGH_DETAIL = 1;

float const INITIAL_CONDITION = 0.00;
string const INITIAL_Q_FILENAME = "-";
float const QCI = 1.00;
float const ETACI =	0.02;
float const COURANTNO =	0.20;
string const OUTPUT_FILENAME = "test_001.nc";
float const SAVE_TIMESTEP = 1.00;
float const FINAL_TIME = 100;
bool const INCLUDE_DENSITY = false;
float const SALINITY_BETA = 0.00;
float const KH = 50.00;
string const SALINITY_FILENAME = "Initial_along_estuary_salinity.csv";
string const SEDIMENT_FILENAME = "sediments.csv";
bool const MCCOMARCK_LIMITER_FLUX = true;
bool const SURFACE_GRADIENT_METHOD = true;
bool const SOURCE_TERM_BALANCE = true;
bool const BETA = false;
bool const DRY_BED = true;
bool const MURILLO_CONDITION = true;
int const INITIAL_BOUNDARY_CONDITION = 1;
int const FINAL_BOUNDARY_CONDITION = 0;
float const FIX_Q = 0;
string const TIDE_FILENAME = "tides.csv";
float const DT_MIN = 10000000000;
int const FORMULA_LIMITER_FLUX = 4;
int const PSI_FORMULA = 2.00;
float const DELTA_VALUE = 0.20;

// clock_t is a signed long: see <time.h>
long const CLOCK_T_MIN = LONG_MIN;
double const CLOCK_T_RANGE = static_cast<double>(LONG_MAX) - static_cast<double>(CLOCK_T_MIN);

// Intel x86, byte order is little-endian, 32-bit
string const PLATFORM = "Intel x86/GNU C++";

string const PROGRAM_NAME = "Saint Venant Equations Solver version 0.0.1 (15 Aug 2024)";
string const PROGRAM_NAME_SHORT = "SV";
string const SV_INI = "configuration.ini";

string const NOTE = "      Note ";
string const COPYRIGHT = "(C) 2024 Manuel Cobos";
string const LINE = "-------------------------------------------------------------------------------";
string const DISCLAIMER1 = "This program is distributed in the hope that it will be useful. but WITHOUT ANY";
string const DISCLAIMER2 = "WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A";
string const DISCLAIMER3 = "PARTICULAR PURPOSE. See the GNU General Public License for more details. You";
string const DISCLAIMER4 = "should have received a copy of the GNU General Public License along with this";
string const DISCLAIMER5 = "program; if not. contact the Free Software Foundation. Inc.. 675 Mass Ave.";
string const DISCLAIMER6 = "Cambridge. MA 02139. USA.";

string const START_NOTICE = "- Started on ";
string const INITIALIZING_NOTICE = "- Initializing";
string const RUN_END_NOTICE = "- Run ended at ";

string const ERROR_NOTICE = " with error code ";
string const ERR = "*** ERROR ";

// Return codes
int const RTN_OK = 0;
int const RTN_ERR_INI = 1;
int const RTN_ERR_SVDIR = 2;
int const RTN_ERR_LOGFILE = 10;
int const RTN_ERR_OUTFILE = 11;


// Not likely that user will need to change these
int const BUF_SIZE = 2048;                                     // Max length (inc. terminating NULL) of any C-type string

string const OUTEXT = ".out";
string const LOGEXT = ".log";

bool bIsStringValidDouble(string &);
bool bIsStringValidInt(string &);

#endif // MAIN_H