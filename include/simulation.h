/*!
*
 * \class CSimulation
 * \brief This class runs Saint Venant simulations
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

#include <vector>
using std::pmr::vector;

class CSimulation
{
private:

    bool bReadConfigurationFile(void);
    bool bOpenLogFile(void);
    static bool bParseDate(string const*, int&, int&, int&);
    static bool bParseTime(string const*, int&, int&, int&);
    static vector<string> *VstrSplit(string const*, char const, vector<string>*);
    static vector<string> VstrSplit(string const*, char const);


    //! The SV folder
    string m_strSVDir;

    //! Folder for the SV .ini file
    string m_strSVIni;

    //! An email address to which to send end-of-simulation messages
    string m_strMailAddress;

    //! Folder in which the CME data file is found
    string m_strDataPathName;

    //! The name of this simulation
    string m_strRunName;

    //! Name of output log file
    string m_strLogFile;

    //! Path for all output files
    string m_strOutPath;

    //! Name of main output file
    string m_strOutFile;

    //! Name of along channel geometry file
    string m_strAlongChannelGeometryFilename;
    //! Name of cross sections file
    string m_strCrossSectionsFilename;
    //! Name of hydro file
    string m_strHydroFilename;

    //! Start time of the simulation (seconds)
    int m_nSimStartSec;

    //! Start time of the simulation (minutes)
    int m_nSimStartMin;

    //! Start time of the simulation (hours)
    int m_nSimStartHour;

    //! Start date of the simulation (day)
    int m_nSimStartDay;

    //! Start date of the simulation (month)
    int m_nSimStartMonth;

    //! Start date of the simulation (year)
    int m_nSimStartYear;

    //! The level of detail in the log file output. Can be LOG_FILE_LOW_DETAIL, LOG_FILE_MIDDLE_DETAIL, or LOG_FILE_HIGH_DETAIL
    int m_nLogFileDetail;

    //! The number of the current iteration (time step)
    unsigned long m_ulIter;

    //! The main output file stream
    ofstream OutStream;

    // Utility routines
    // Utility routines
    static void AnnounceStart(void);
    void AnnounceLicence(void);
    // void AnnounceReadBasementDEM(void) const;
    // static void AnnounceAddLayers(void);
    // static void AnnounceReadRasterFiles(void);
    // static void AnnounceReadVectorFiles(void);
    // void AnnounceReadLGIS(void) const;
    // void AnnounceReadICGIS(void) const;
    // void AnnounceReadIHGIS(void) const;
    static void AnnounceInitializing(void);
    // void AnnounceReadInitialSuspSedGIS(void) const;
    // void AnnounceReadInitialFineUnconsSedGIS(int const) const;
    // void AnnounceReadInitialSandUnconsSedGIS(int const) const;
    // void AnnounceReadInitialCoarseUnconsSedGIS(int const) const;
    // void AnnounceReadInitialFineConsSedGIS(int const) const;
    // void AnnounceReadInitialSandConsSedGIS(int const) const;
    // void AnnounceReadInitialCoarseConsSedGIS(int const) const;
    // void AnnounceReadDeepWaterWaveValuesGIS(void) const;
    // void AnnounceReadSedimentEventInputValuesGIS(void) const;
    // void AnnounceReadFloodLocationGIS(void) const;
    // void AnnounceReadTideData(void) const;
    // static void AnnounceReadSCAPEShapeFunctionFile(void);
    // static void AnnounceAllocateMemory(void);
    // static void AnnounceIsRunning(void);
    // static void AnnounceSimEnd(void);
    void StartClock(void);
    bool bFindExeDir(char const*);
    // bool bTimeToQuit(void);
    // static int nDoTimeUnits(string const*);
    // int nDoSimulationTimeMultiplier(string const*);
    // static double dGetTimeMultiplier(string const*);
    // static bool bParseDate(string const*, int&, int&, int&);
    // static bool bParseTime(string const*, int&, int&, int&);
    // void DoTimestepTotals(void);
    static string strGetBuild(void);
    // static string strGetComputerName(void);
    static string strGetErrorText(int const);
    static string strTrim(string const*);
    static string strTrimLeft(string const*);
    static string strTrimRight(string const*);


    //! System start-simulation time
    time_t m_tSysStartTime;

    //! System finish-simulation time
    time_t m_tSysEndTime;

    //! Last value returned by clock()
    double m_dClkLast;
    //! Total elapsed CPU time
    double m_dCPUClock;

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