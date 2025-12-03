/*!
*
 * \file hydrograph.cpp
 * \brief Lateral inflow (tributary/point source) handling
 * \details Manages time series of lateral water inflows entering the main channel.
 *          Interpolates flow rates and locates nearest cross-section for each source.
 * \author Manuel Cobos Budia

 * \date 28/08/2024
 * \copyright GNU General Public License
 *
 */

/*==============================================================================================================================

==============================================================================================================================*/
#include <string>
#include <iostream>  // Para std::cerr y std::endl
#include <vector>    // Para std::vector (buena práctica incluirlo explícitamente)

using std::to_string;
using std::string;
using std::vector;
using std::cerr;
using std::cout;
using std::endl;

#include "hydrograph.h"

//===============================================================================================================================
//! The CHydrograph constructor
//===============================================================================================================================
CHydrograph::CHydrograph() : m_dHydrographXLocation(0.0), 
                            m_dHydrographYLocation(0.0),
                            m_nNearestCrossSectionNo(-1) {
    // Reservar espacio inicial para evitar reasignaciones frecuentes
    m_vHydroTime.reserve(1000);         // Reservar para ~1000 timesteps
    m_vHydroWaterFlow.reserve(1000);    // Reservar para ~1000 timesteps
}

//===============================================================================================================================
//! The CHydrograph destructor
//===============================================================================================================================
CHydrograph::~CHydrograph() = default;

//===============================================================================================================================
//! Append a hydrograph object to the CSimulation object
//===============================================================================================================================
void CHydrograph::dAppend2Vector(const string& strItem, double dValue) {
    if (strItem == "time" || strItem == "TIME") {
        m_vHydroTime.push_back(dValue);
    }
    else if (strItem == "flow" || strItem == "FLOW" || strItem == "Q" || strItem == "water flow") {
        m_vHydroWaterFlow.push_back(dValue);
    }
    else {
        cerr << "Warning: Unknown item '" << strItem << "' in hydrograph data" << endl;
    }
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
const vector<double>& CHydrograph::vGetTime() const {
    return m_vHydroTime;
}

const vector<double>& CHydrograph::vGetQ() const {
    return m_vHydroWaterFlow;
}

size_t CHydrograph::size() const {
    return m_vHydroTime.size();
}

bool CHydrograph::isValid() const {
    return m_vHydroTime.size() == m_vHydroWaterFlow.size() && !m_vHydroTime.empty();
}

void CHydrograph::clear() {
    m_vHydroTime.clear();
    m_vHydroWaterFlow.clear();
}

void CHydrograph::addTimeFlowPair(double time, double flow) {
    m_vHydroTime.push_back(time);
    m_vHydroWaterFlow.push_back(flow);
}

double CHydrograph::getFlowAtTime(double time) const {
    if (!isValid()) {
        cerr << "Error: Invalid hydrograph data" << endl;
        return 0.0;
    }
    
    // Interpolación lineal simple
    for (size_t i = 0; i < m_vHydroTime.size() - 1; ++i) {
        if (time >= m_vHydroTime[i] && time <= m_vHydroTime[i + 1]) {
            double ratio = (time - m_vHydroTime[i]) / (m_vHydroTime[i + 1] - m_vHydroTime[i]);
            return m_vHydroWaterFlow[i] + ratio * (m_vHydroWaterFlow[i + 1] - m_vHydroWaterFlow[i]);
        }
    }
    
    // Si está fuera del rango, devolver el valor más cercano
    if (time < m_vHydroTime.front()) {
        return m_vHydroWaterFlow.front();
    } else {
        return m_vHydroWaterFlow.back();
    }
}