/*!
*
 * \class CDataReader
 * \brief
 * \details Description of CDataReader class which contains the methods for reading input files.
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

class CDataReader {

    //! The CSimulation class is a friend of the CDataReader class
    friend class CSimulation;

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

    //! The number of the current iteration (time step)
    unsigned long m_ulIter;

    //! Trailing whitespace
    static string strTrim(string const*);
    static string strTrimLeft(string const*);
    static string strTrimRight(string const*);

    //! Formatting Date and Time
    static bool bParseDate(string const*, int&, int&, int&);
    static bool bParseTime(string const*, int&, int&, int&);

    //! Constants to determine simulation time-steps
    int nSimulationTimeMultiplier(string const*);
    static double dGetTimeMultiplier(string const*);
    static int nDoTimeUnits(string const*);


public:

     CDataReader();
    ~CDataReader();

    // Log file detail level
    static constexpr int NO_LOG_FILE = 0;
    static constexpr int LOG_FILE_LOW_DETAIL = 1;
    static constexpr int LOG_FILE_HIGH_DETAIL = 2;


    string const SV_INI = "configuration.ini";

    string const OUT_EXT = ".nc";
    string const LOG_EXT = ".log";

    //! Read configuration file with global data
    bool bReadConfigurationFile(CSimulation* m_pSimulation);

    //! Read along channel data and the initial estuarine condition
    bool bReadAlongChannelDataFile(CSimulation* m_pSimulation) const;

    //! Read every cross-section geometry
    bool bReadCrossSectionGeometryFile(CSimulation* m_pSimulation) const;

    //! Read time-series with upward and downward boundary conditions
    static bool bReadUpwardBoundaryConditionFile(CSimulation* m_pSimulation);
    static bool bReadDownwardBoundaryConditionFile(CSimulation* m_pSimulation);

    //! Read input hydrographs
    bool bReadHydrographsFile(CSimulation* m_pSimulation) const;

    //! Open log file for writing
    static bool bOpenLogFile(CSimulation* m_pSimulation);

};
#endif // DATA_READER_H
