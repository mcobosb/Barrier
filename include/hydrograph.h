/*!
*
 * \class CHydrograph
 * \brief
 * \details TODO 001 This is a more detailed description of the CDataReader class
 * \author Manuel Cobos Budia

 * \date 2024
 * \copyright GNU General Public License
 *
 * \file hydro.h
 * \brief Contains DataReader definitions
 *
 */
#ifndef HYDROGRAPH_H
#define HYDROGRAPH_H

#include <simulation.h>
#include <string>
using std::string;

#include <vector>

//class CSimulation;
//class CEstuary;

class CHydrograph {

//    //! The CSimulation class is a friend of the CDataReader class
//    friend class CSimulation;
//
//    //! The CEstuary class is a friend of the CDataReader class
//    friend class CEstuary;

private:
    //! X coordinate of input hydro
    double m_dHydrographXLocation;

    //! Y coordinate of input hydro
    double m_dHydrographYLocation;

    //! Nearest cross-section estuary to flow
    int m_nCrossSectionNo;

    //! Time from starting simulation in hours and water discharge
    vector<double> m_vHydroTime;

    //! Water discharge at time
    vector<double> m_vHydroWaterFlow;

public:
    CHydrograph();
    ~CHydrograph()();

    double GetHydrographXLocation();
    double GetHydrographYLocation();
    int GetCrossSectionNo();
    vector<double> GetHydrographTime();
    vector<double> GetHydrographWaterFlow();

    void SetHydrographXLocation(double dHydroXLocation);
    void SetHydrographYLocation(double dHydroYLocation);
    void SetCrossSectionNo(int nCrossSectionNo);
    void dAppend2Vector(string strItem, double dValue);

};

#endif // HYDROGRAPH_H