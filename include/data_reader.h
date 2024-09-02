/*!
*
 * \class CDataReader
 * \brief
 * \details TODO 001 This is a more detailed description of the CDataReader class
 * \author Manuel Cobos Budia

 * \date 2024
 * \copyright GNU General Public License
 *
 * \file data_reader.h
 * \brief Contains DataReader definitions
 *
 */
#ifndef DATA_READER_H
#define DATA_READER_H

#include <simulation.h>
#include <string>
using std::string;

class CSimulation;
class CEstuary;

class CDataReader {

    //! The CSimulation class is a friend of the CDataReader class
    friend class CSimulation;

    //! The CEstuary class is a friend of the CDataReader class
    friend class CEstuary;

private:

    //! Folder for the SV .ini file
    string m_strSVIni;

    //! An email address to which to send end-of-simulation messages
    string m_strMailAddress;

    //! The name of this simulation
    string m_strRunName;

    //! Name of output log file
    string m_strLogFile;

    //! Path for all output files
    string m_strOutPath;

    //! Name of main output file
    string m_strOutFile;

    //! Name of along channel geometry file
    string m_strAlongChannelDataFilename;

    //! Name of cross-sections file
    string m_strCrossSectionsFilename;

    //! Name of the sediment properties
    string m_strSedimentPropertiesFilename;

    //! Name of downward boundary tidal elevation file
    string m_strTidalFilename;

    //! Name of downward boundary water flow file
    string m_strWaterFlowFilename;


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

    //! Multiplier for duration units, to convert to hours
    double m_dDurationUnitsMultiplier;

    //! The duration units for this simulation
    string m_strDurationUnits;

    //! The level of detail in the log file output. Can be LOG_FILE_LOW_DETAIL, LOG_FILE_MIDDLE_DETAIL, or LOG_FILE_HIGH_DETAIL
    int m_nLogFileDetail;

    //! The number of the current iteration (time step)
    unsigned long m_ulIter;

    static string strTrim(string const*);
    static string strTrimLeft(string const*);
    static string strTrimRight(string const*);
    static bool bParseDate(string const*, int&, int&, int&);
    static bool bParseTime(string const*, int&, int&, int&);
    int nSimulationTimeMultiplier(string const*);
    double dGetTimeMultiplier(string const*);
    int nDoTimeUnits(string const*);

    // A pointer to the CSimulation object
    // CSimulation* m_pSimulation;

    // A pointer to the CEstuary object
    // CEstuary* m_pEstuary;

public:

     CDataReader();
    ~CDataReader();

    // CSimulation* pGetSimulation();

    // Log file detail level
    static constexpr int NO_LOG_FILE = 0;
    static constexpr int LOG_FILE_HIGH_DETAIL = 1;


    string const SV_INI = "configuration.ini";

    string const OUT_EXT = ".out";
    string const LOG_EXT = ".log";

    bool bReadConfigurationFile(CSimulation* m_pSimulation);
    bool bReadAlongChannelDataFile(CSimulation* m_pSimulation);
    bool bReadCrossSectionGeometryFile(CSimulation* m_pSimulation);
    bool bReadDownwardBoundaryConditionFile(CSimulation* m_pSimulation);
    bool bReadHydrographsFile(CSimulation* m_pSimulation);
    bool bOpenLogFile(CSimulation* m_pSimulation);

};
#endif // DATA_READER_H
