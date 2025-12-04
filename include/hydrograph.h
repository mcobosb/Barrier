/*!
*
 * \class CHydrograph
 * \details Represents lateral inflow sources (tributaries, point sources) along the channel.
 *          Handles interpolation of flow rates over time.
 * \author Manuel Cobos Budia

 * \date 2026
 * \copyright GNU General Public License
 *
 * \file hydrograph.h
 * \brief Contains DataReader definitions
 *
 */
#ifndef HYDROGRAPH_H
#define HYDROGRAPH_H

#include <string>
using std::string;

#include <vector>
using std::vector;


class CHydrograph {


private:
    //! X coordinate of input hydro
    double m_dHydrographXLocation;

    //! Y coordinate of input hydro
    double m_dHydrographYLocation;


    //! Time from starting simulation in hours and water discharge
    vector<double> m_vHydroTime;

    //! Water discharge at time
    vector<double> m_vHydroWaterFlow;

public:
    CHydrograph();
    ~CHydrograph();

    //! Nearest cross-section estuary to flow
    int m_nNearestCrossSectionNo;

    [[nodiscard]] double dGetHydrographXLocation() const;
    [[nodiscard]] double dGetHydrographYLocation() const;

    void dSetHydrographXLocation(double dHydroXLocation);
    void dSetHydrographYLocation(double dHydroYLocation);
    void dAppend2Vector(const string& strItem, double dValue);

    //! Vector getter - devolver referencias constantes para evitar copias
    [[nodiscard]] const vector<double>& vGetTime() const;
    [[nodiscard]] const vector<double>& vGetQ() const;

    //! Get the size of the time series
    [[nodiscard]] size_t size() const;

    //! Check if the vectors are synchronized (same size)
    [[nodiscard]] bool isValid() const;

    //! Clear all data
    void clear();

    //! Add a time-flow pair
    void addTimeFlowPair(double time, double flow);

    //! Get interpolated flow at given time
    [[nodiscard]] double getFlowAtTime(double time) const;
};

#endif // HYDROGRAPH_H