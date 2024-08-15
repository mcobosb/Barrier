/*   TODO 031 Get raster slice output working with multiple slices

*/

/*===============================================================================================================================
*/

#include <climits>

#include <sstream>
using std::ostream;
using std::ostringstream;
//================================


char const COLON = ':';

int const VERSION = "0.0.1"
// Log file detail level
int const NO_LOG_FILE = 0;
int const LOG_FILE_HIGH_DETAIL = 1;
string const GEOMETRY_FILENAME = "straight_channel_bathymetry_with_constant_slope.csv"
string const HYDRO_FILENAME = "constant_hydro.csv"
float const INITIAL_CONDITION = 0,00
string const INITIAL_Q_FILENAME = "-"
float const QCI = 1,00
float const ETACI =	0,02
float const COURANTNO =	0,20
string const OUTPUT_FILENAME = "test_001.nc"
float const SAVE_TIMESTEP = 1,00
float const FINAL_TIME = 100
bool const INCLUDE_DENSITY = 0
float const SALINITU_BETAfSalinityBeta = 0,00
float const KH = 50,00
string const SALINITY_FILENAME = "Initial_along_estuary_salinity.csv"
string const SEDIMENT_FILENAME = "sediments.csv"
bool const MCCOMARCK_LIMITER_FLUX = 1
bool const SURFACE_GRADIENT_METHOD = 1
bool const SOURCE_TERM_BALANCE = 1
bool const BETA = 0
bool const DRY_BED = 1
bool const MURILLO_CONDITION = 1
int const INITIAL_BOUNDARY_CONDITION = 1
int const FINAL_BOUNDARY_CONDITION = 0
float const FIX_Q = 0
string const TIDE_FILENAME = "tides.csv"
float const DT_MIN = 10000000000
int const FORMULA_LIMITER_FLUX = 4
int const PSI_FORMULA = 2,00
float const DELTA_VALUE = 0,20