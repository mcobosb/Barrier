#ifndef MAIN_H
#define MAIN_H
/*===============================================================================================================================
*/

#include <string>
using std::string;

#include <climits>

#include <sstream>
using std::ostream;
using std::ostringstream;
//================================



// Global constant
//! gravitational constant
inline double G = 9.81;
//! minimum cross-section area
inline double DRY_AREA = 1e-8;
//! minimum cross-section water flow (to ensure that U = Q/A ~ 0)
inline double DRY_Q = 1e-5;


// float const INITIAL_CONDITION = 0.00;
// string const INITIAL_Q_FILENAME = "initial_along_estuary_water_flow.csv";
// float const QCI = 1.00;
// float const ETACI =	0.02;
// float const COURANT_NUMBER = 0.20;
inline float CDX_MURILLO = 0.6;
string const OUTPUT_FILENAME = "test_001.nc";
// float const SAVE_TIMESTEP = 1.00;
// float const FINAL_TIME = 100;
// bool const DO_DENSITY = false;
float const SALINITY_BETA = 0.00;
float const KH = 50.00;
string const SALINITY_FILENAME = "initial_along_estuary_salinity.csv";
string const SEDIMENT_FILENAME = "sediments_properties.csv";
// bool const DO_MCCOMARCK_LIMITER_FLUX = true;
// bool const DO_SURFACE_GRADIENT_METHOD = true;
// bool const DO_SOURCE_TERM_BALANCE = true;
// bool const DO_BETA = false;
// bool const DO_DRY_BED = true;
// bool const DO_MURILLO_CONDITION = true;
// int const INITIAL_BOUNDARY_CONDITION = 1;
// int const FINAL_BOUNDARY_CONDITION = 0;
// float const FIX_Q = 0;
// string const TIDE_FILENAME = "tides.csv";
// float const DT_MIN = 10000000000.0;
int const FORMULA_LIMITER_FLUX = 4;
int const PSI_FORMULA = 2.00;
float const DELTA_VALUE = 0.20;


#endif // MAIN_H