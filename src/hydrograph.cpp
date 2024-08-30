/*!
*
 * \file hydro.cpp
 * \brief
 * \details TODO 001 A more detailed description of these routines.
 * \author Manuel Cobos Budia

 * \date 28/08/2024
 * \copyright GNU General Public License
 *
 */

/*==============================================================================================================================

==============================================================================================================================*/
//#include <cstdlib> // for strtol() and strtod()
//#include <fstream>
//using std::ifstream;
//
//#include <sstream>
//using std::stringstream;
//
//#include <iostream>
//using std::cerr;
//using std::cout;
//using std::endl;
//using std::ios;
//
#include <string>
using std::to_string;
//
//#include <algorithm>
//using std::find;

#include <hydrograph.h>

//===============================================================================================================================
//! The CHydro constructor
//===============================================================================================================================
CHydro::CHydro(){
    m_dHydroXLocation =
    m_dHydroYLocation = 0.0;

    m_nCrossSectionNo = 0;
}

//===============================================================================================================================
//! The CHydro destructor
//===============================================================================================================================
CHydro::~CHydro() = default;


//===============================================================================================================================
//! Append a hydrogram object to the CSimulation object
//===============================================================================================================================
void CHydro::dAppend2Vector(string strItem, double dValue){
    if (strItem == "time")
        m_vHydroTime.push_back(dValue);
    else if (strItem == "water flow")
        m_vHydroWaterFlow.push_back(dValue);

}