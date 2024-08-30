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
//! The CHydrograph constructor
//===============================================================================================================================
CHydrograph::CHydrograph(){
    m_dHydrographXLocation =
    m_dHydrographYLocation = 0.0;

    m_nNearestCrossSectionNo = 0;
}

//===============================================================================================================================
//! The CHydrograph destructor
//===============================================================================================================================
CHydrograph::~CHydrograph() = default;


//===============================================================================================================================
//! Append a hydrograph object to the CSimulation object
//===============================================================================================================================
void CHydrograph::dAppend2Vector(const string& strItem, const double dValue){
    if (strItem == "time")
        m_vHydroTime.push_back(dValue);
    else if (strItem == "water flow")
        m_vHydroWaterFlow.push_back(dValue);
}


//===============================================================================================================================
//! Set the hydrograph X coordinate location
//===============================================================================================================================
void CHydrograph::dSetHydrographXLocation(double dValue) {
    m_dHydrographXLocation = dValue;
}

//===============================================================================================================================
//! Get the hydrograph X coordinate location
//===============================================================================================================================
double CHydrograph::dGetHydrographXLocation() {
    return m_dHydrographXLocation;
}

//===============================================================================================================================
//! Set the hydrograph Y coordinate location
//===============================================================================================================================
void CHydrograph::dSetHydrographYLocation(double dValue) {
    m_dHydrographYLocation = dValue;
}


//===============================================================================================================================
//! Get the hydrograph Y coordinate location
//===============================================================================================================================
double CHydrograph::dGetHydrographYLocation() {
    return m_dHydrographYLocation;
}