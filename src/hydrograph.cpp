/*!
*
 * \file hydrograph.cpp
 * \brief
 * \details TODO 001 A more detailed description of these routines.
 * \author Manuel Cobos Budia

 * \date 28/08/2024
 * \copyright GNU General Public License
 *
 */

/*==============================================================================================================================

==============================================================================================================================*/
#include <string>
using std::to_string;

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
void CHydrograph::dSetHydrographXLocation(const double dHydroXLocation) {
    m_dHydrographXLocation = dHydroXLocation;
}

//===============================================================================================================================
//! Get the hydrograph X coordinate location
//===============================================================================================================================
double CHydrograph::dGetHydrographXLocation() const {
    return m_dHydrographXLocation;
}

//===============================================================================================================================
//! Set the hydrograph Y coordinate location
//===============================================================================================================================
void CHydrograph::dSetHydrographYLocation(const double dHydroYLocation) {
    m_dHydrographYLocation = dHydroYLocation;
}


//===============================================================================================================================
//! Get the hydrograph Y coordinate location
//===============================================================================================================================
double CHydrograph::dGetHydrographYLocation() const {
    return m_dHydrographYLocation;
}

//! Getter for vector variables
vector<double> CHydrograph::vGetTime() {
    return m_vHydroTime;
}

vector<double> CHydrograph::vGetQ() {
    return m_vHydroWaterFlow;
}