/*!
*
 * \class CScreenPresenter
 * \brief
 * \details TODO 001 This is a more detailed description of the ScreenPresenter class
 * \author Manuel Cobos Budia

 * \date 2024
 * \copyright GNU General Public License
 *
 * \file screen_presenter.h
 * \brief Contains CScreenPresenter definitions
 *
 */

#include <string>
using std::string;

#include <vector>
using std::vector;

#include <sstream>
using std::ostream;
using std::ostringstream;

#include <fstream>
using std::ofstream;

#include "simulation.h"


class CScreenPresenter
    {
private:

    //! The SV folder
    string m_strSVDir;

    //! The main output file stream
    ofstream OutStream;

    //! System finish-simulation time
    time_t m_tSysEndTime{};

    //! Last value returned by clock()
    double m_dClkLast{};
    //! Total elapsed CPU time
    double m_dCPUClock{};

    // Utility routines
    void StartClock(CSimulation *m_pSimulation);

    bool bFindExeDir(char const*);
    // bool bTimeToQuit(void);
    // static int nDoTimeUnits(string const*);
    // int nDoSimulationTimeMultiplier(string const*);
    // static double dGetTimeMultiplier(string const*);
    // static bool bParseDate(string const*, int&, int&, int&);
    // static bool bParseTime(string const*, int&, int&, int&);
    // void DoTimestepTotals(void);
    static string strGetBuild();
    static string strGetComputerName();

    CSimulation *m_pSimulation;


public:
    ofstream LogStream;

    CScreenPresenter();
    ~CScreenPresenter();

    //! Carries out init-of-simulation tidying (error messages etc.)
    void StartingRun(int, char const* [], CSimulation* m_pSimulation);
    void EndingRun();

    void AnnounceStart();
    static void AnnounceLicence(const CSimulation *m_pSimulation) ;
    // void AnnounceReadBasementDEM() const;
    // static void AnnounceAddLayers();
    // static void AnnounceReadRasterFiles();
    // static void AnnounceReadVectorFiles();
    // void AnnounceReadLGIS() const;
    // void AnnounceReadICGIS() const;
    // void AnnounceReadIHGIS() const;
    static void AnnounceInitializing();
    // void AnnounceReadInitialSuspSedGIS() const;
    // void AnnounceReadInitialFineUnconsSedGIS(int const) const;
    // void AnnounceReadInitialSandUnconsSedGIS(int const) const;
    // void AnnounceReadInitialCoarseUnconsSedGIS(int const) const;
    // void AnnounceReadInitialFineConsSedGIS(int const) const;
    // void AnnounceReadInitialSandConsSedGIS(int const) const;
    // void AnnounceReadInitialCoarseConsSedGIS(int const) const;
    // void AnnounceReadDeepWaterWaveValuesGIS() const;
    // void AnnounceReadSedimentEventInputValuesGIS() const;
    // void AnnounceReadFloodLocationGIS(void) const;
    // void AnnounceReadTideData(void) const;
    // static void AnnounceReadSCAPEShapeFunctionFile(void);
    // static void AnnounceAllocateMemory(void);
    // static void AnnounceIsRunning(void);
    // static void AnnounceSimEnd(void);
};
