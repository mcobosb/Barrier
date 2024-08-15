/*!
*
 * \class CSimulation
 * \brief This class runs CoastalME simulations
 * \details TODO 001 This is a more detailed description of the CSimulation class
 * \author Manuel Cobos Budia

 * \date 2024
 * \copyright GNU General Public License
 *
 * \file simulation.h
 * \brief Contains CSimulation definitions
 *
 */

#ifndef SIMULATION_H
#define SIMULATION_H
/*===============================================================================================================================

===============================================================================================================================*/
#include <ctime>
using std::localtime;
using std::time;
using std::time_t;

#include <fstream>
using std::ofstream;

#include <string>
using std::string;


class CSimulation
{
private:
    //! The CME folder
    string m_strCMEDir;

    //! Folder for the CME .ini file
    string m_strCMEIni;

    //! An email addresx to which to send end-of-simulation messages
    string m_strMailAddress;

    // Utility routines
    static void AnnounceStart(void);
    void AnnounceLicence(void);
    void AnnounceReadBasementDEM(void) const;
    static void AnnounceAddLayers(void);
    static void AnnounceReadRasterFiles(void);
    static void AnnounceReadVectorFiles(void);
    void AnnounceReadLGIS(void) const;
    void AnnounceReadICGIS(void) const;
    void AnnounceReadIHGIS(void) const;
    static void AnnounceInitializing(void);

public:
    ofstream LogStream;

    CSimulation(void);
    ~CSimulation(void);

    //! Runs the simulation
    int nDoSimulation(int, char const* []);

    //! Carries out end-of-simulation tidying (error messages etc.)
    void DoSimulationEnd(int const);
};
#endif // SIMULATION_H
}