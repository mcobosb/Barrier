/*!
*
 * \class CHydrograph
 * \details TODO 001 This is a more detailed description of the CDataReader class
 * \author Manuel Cobos Budia

 * \date 2024
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

    //! Vector getter
    vector<double> vGetTime();
    vector<double> vGetQ() ;

};

#endif // HYDROGRAPH_H