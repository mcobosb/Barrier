/*!
*
 * \file hydrograph.cpp
 * \brief Lateral inflow (tributary/point source) handling
 * \details Manages time series of lateral water inflows entering the main channel.
 *          Interpolates flow rates and locates nearest cross-section for each source.
 * \author Manuel Cobos Budia

 * \date 2026
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

/**
 * @brief Construct a new CHydrograph object for lateral inflow time series
 * 
 * Initializes:
 * - Location coordinates to (0.0, 0.0)
 * - Nearest cross-section to -1 (unassigned)
 * - Reserves 1000 elements for time/flow vectors (optimization)
 * 
 * @note Reserving space prevents frequent reallocation during CSV parsing
 */
CHydrograph::CHydrograph() : m_dHydrographXLocation(0.0), 
                            m_dHydrographYLocation(0.0),
                            m_nNearestCrossSectionNo(-1) {
    // Reserve space to avoid frequent reallocations during data loading
    m_vHydroTime.reserve(1000);         // ~1000 timesteps typical for monthly data
    m_vHydroWaterFlow.reserve(1000);    // ~1000 flow values
}

/**
 * @brief Destructor (default implementation)
 */
CHydrograph::~CHydrograph() = default;

/**
 * @brief Append data item to hydrograph time series
 * 
 * Case-insensitive column recognition:
 * - "time" or "TIME" → m_vHydroTime
 * - "flow", "FLOW", "Q", "water flow" → m_vHydroWaterFlow
 * 
 * @param strItem Column name from CSV header
 * @param dValue Numerical value to append
 * 
 * @warning Prints warning for unrecognized column names
 * @note Called by data parser during CSV reading
 */
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

/**
 * @brief Set X coordinate of lateral inflow location
 * @param dHydroXLocation X position (m) along main channel axis
 */
void CHydrograph::dSetHydrographXLocation(const double dHydroXLocation) {
    m_dHydrographXLocation = dHydroXLocation;
}

/**
 * @brief Get X coordinate of lateral inflow location
 * @return X position (m)
 */
double CHydrograph::dGetHydrographXLocation() const {
    return m_dHydrographXLocation;
}

/**
 * @brief Set Y coordinate of lateral inflow location
 * @param dHydroYLocation Y position (m) perpendicular to main axis
 */
void CHydrograph::dSetHydrographYLocation(const double dHydroYLocation) {
    m_dHydrographYLocation = dHydroYLocation;
}

/**
 * @brief Get Y coordinate of lateral inflow location
 * @return Y position (m)
 */
double CHydrograph::dGetHydrographYLocation() const {
    return m_dHydrographYLocation;
}

/**
 * @brief Get time vector (const reference for efficiency)
 * @return Vector of time values (seconds since simulation start)
 */
const vector<double>& CHydrograph::vGetTime() const {
    return m_vHydroTime;
}

/**
 * @brief Get flow rate vector (const reference for efficiency)
 * @return Vector of discharge values (m³/s)
 */
const vector<double>& CHydrograph::vGetQ() const {
    return m_vHydroWaterFlow;
}

/**
 * @brief Get number of time-flow pairs
 * @return Size of time series
 */
size_t CHydrograph::size() const {
    return m_vHydroTime.size();
}

/**
 * @brief Validate hydrograph data consistency
 * @return true if time and flow vectors have same size and are non-empty
 */
bool CHydrograph::isValid() const {
    return m_vHydroTime.size() == m_vHydroWaterFlow.size() && !m_vHydroTime.empty();
}

/**
 * @brief Clear all hydrograph data
 * @note Useful for reloading or resetting simulation
 */
void CHydrograph::clear() {
    m_vHydroTime.clear();
    m_vHydroWaterFlow.clear();
}

/**
 * @brief Add a time-flow data pair
 * @param time Time value (seconds)
 * @param flow Discharge (m³/s)
 */
void CHydrograph::addTimeFlowPair(double time, double flow) {
    m_vHydroTime.push_back(time);
    m_vHydroWaterFlow.push_back(flow);
}

/**
 * @brief Interpolate flow rate at arbitrary time using linear interpolation
 * 
 * Interpolation rules:
 * - Within data range: Linear interpolation between adjacent points
 * - Before first point: Return first value (constant extrapolation)
 * - After last point: Return last value (constant extrapolation)
 * 
 * @param time Query time (seconds since simulation start)
 * @return Interpolated flow rate (m³/s), or 0.0 if invalid data
 * 
 * @note O(n) search - for frequent calls, consider binary search or caching
 * @warning Returns 0.0 and prints error if isValid() == false
 */
double CHydrograph::getFlowAtTime(double time) const {
    if (!isValid()) {
        cerr << "Error: Invalid hydrograph data" << endl;
        return 0.0;
    }
    
    // Simple linear interpolation
    for (size_t i = 0; i < m_vHydroTime.size() - 1; ++i) {
        if (time >= m_vHydroTime[i] && time <= m_vHydroTime[i + 1]) {
            double ratio = (time - m_vHydroTime[i]) / (m_vHydroTime[i + 1] - m_vHydroTime[i]);
            return m_vHydroWaterFlow[i] + ratio * (m_vHydroWaterFlow[i + 1] - m_vHydroWaterFlow[i]);
        }
    }
    
    // Outside range: return nearest value (constant extrapolation)
    if (time < m_vHydroTime.front()) {
        return m_vHydroWaterFlow.front();
    } else {
        return m_vHydroWaterFlow.back();
    }
}