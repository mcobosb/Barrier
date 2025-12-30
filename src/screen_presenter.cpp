/*!
 *
 * \file screen_presenter.cpp
 * \brief Console output and progress display
 * \details Manages terminal output for simulation progress, timing,
 *          system information, and user notifications.
 * \author Manuel Cobos Budia

 * \date 2026
 * \copyright GNU General Public License
 *
 */

/*==============================================================================================================================

==============================================================================================================================*/
// #include <assert.h>

#include <cmath>
using std::floor;

#include <ctime>
using std::clock;
using std::clock_t;
using std::difftime;
using std::localtime;
using std::time;

#include <cstdlib>
#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#include <sys/utsname.h>
#endif

#include <iostream>
using std::cerr;
using std::cout;
using std::endl;
using std::ios;

#include <iomanip>
using std::put_time;
using std::resetiosflags;
using std::setiosflags;
using std::setprecision;
using std::setw;

#include <string>
using std::to_string;

#include <sstream>
using std::stringstream;

#include <algorithm>
using std::all_of;
using std::transform;

#include "screen_presenter.h"
#include "utils.h"
#include "simulation.h"
#include "main.h"



/**
 * @brief Construct a new CScreenPresenter object for console I/O
 * 
 * Initializes system end time to 0 (will be set when simulation completes).
 */
CScreenPresenter::CScreenPresenter() {
   m_tSysEndTime = 0;
}

/**
 * @brief Destructor (default implementation)
 */
CScreenPresenter::~CScreenPresenter() = default;


/**
 * @brief Initialize console output and announce simulation start
 * 
 * Performs:
 * 1. Print banner with program name, version, build timestamp
 * 2. Start CPU and wall-clock timers
 * 3. Find executable directory
 * 4. Display license and disclaimer
 * 5. Print system info (hostname, start time)
 * 
 * @param nArg Number of command-line arguments (unused)
 * @param pcArgv Array of argument strings (used to find executable path)
 * @param m_pSimulation Pointer to main simulation object
 * 
 * @note Called once at beginning of main()
 */
void CScreenPresenter::StartingRun([[maybe_unused]] int nArg, char const* pcArgv[], CSimulation* m_pSimulation)
{
   // Print program banner
   AnnounceStart(m_pSimulation);

   // Start the clock ticking
   StartClock(m_pSimulation);

   //    return (RTN_ERR_SVDIR);
   bFindExeDir(pcArgv[0]);

   // OK, we are off, tell the user about the licence and the start time
   AnnounceLicence(m_pSimulation);
}



/**
 * @brief Display simulation end message with timestamp
 * 
 * Prints completion message to stdout and log file.
 * Records system end time if not already set (e.g., from exception).
 * 
 * @note Commented code shows historical features:
 * - Progress bar final update (100%)
 * - Elapsed time display
 * - Error message formatting
 * - Email notification (removed)
 */
void CScreenPresenter::EndingRun()
{
   // If we don't know the time that the run ended (e.g. because it did not finish correctly), then get it now
   // static double sdElapsed = 0;
   m_tSysEndTime = time(nullptr);

   // Calculate time elapsed and remaining
   // sdElapsed = difftime(m_tSysEndTime, pSimulation->m_tSysStartTime);

   // if (nRtn  == RTN_OK)
   //    // normal ending
   // {
   // cout << "\r    - Remaining Time: " << std::fixed << setprecision(3) << setw(6) << 0.000 << " s -  Progress: " << std::fixed << setprecision(3) << setw(6) << 100 << '%' << std::flush;
   // cout << "\r    - Elapsed Time: " << std::fixed << setprecision(3) << setw(6) << sdElapsed  << endl;
   cout << RUN_END_NOTICE << put_time(localtime(&m_tSysEndTime), "%H:%M on %A %d %B %Y") << endl;
   // }
   //
   // else
      // Aborting because of some error
   // {
   //    cerr << RUN_END_NOTICE << "iteration " << m_ulIter << ERROR_NOTICE << nRtn << " (" << strGetErrorText(nRtn) << ") on " << std::put_time(localtime(&m_tSysEndTime), "%T %A %d %B %Y") << endl;
   // }

   // Write the error message to the logfile and to stdout
   // if (LogStream && LogStream.is_open())
   // {
   //    LogStream << ERR << strGetErrorText(nRtn) << " (error code " << nRtn << ") on " << std::put_time(localtime(&m_tSysEndTime), "%T %A %d %B %Y") << endl;
   //    LogStream.flush();
   // }
   //
   // if (OutStream && OutStream.is_open())
   // {
   //    OutStream << ERR << strGetErrorText(nRtn) << " (error code " << nRtn << ") on " << std::put_time(localtime(&m_tSysEndTime), "%T %A %d %B %Y") << endl;
   //    OutStream.flush();
   // }
}

/**
 * @brief Get compilation timestamp string
 * 
 * Returns string with format: "(HH:MM:SS MMM DD YYYY build)"
 * or "(HH:MM:SS MMM DD YYYY DEBUG build)" for debug builds.
 * 
 * Uses compiler macros __TIME__ and __DATE__.
 * 
 * @return Build timestamp string
 * @note Useful for identifying executable version in bug reports
 */
string CScreenPresenter::strGetBuild()
{
   string strBuild("(");
   strBuild.append(__TIME__);
   strBuild.append(" ");
   strBuild.append(__DATE__);
#ifdef _DEBUG
   strBuild.append(" DEBUG)");
#endif
   strBuild.append(" build)");

   return strBuild;
}


/**
 * @brief Print program banner to console and log file
 * 
 * Displays: PROGRAM_NAME for PLATFORM (build_timestamp)
 * 
 * @param m_pSimulation Pointer to simulation object (for log file access)
 */
void CScreenPresenter::AnnounceStart(CSimulation* m_pSimulation)
{
   cout << endl
        << PROGRAM_NAME << endl
        << "      for " << PLATFORM << " (" << strGetBuild() << " build)" << endl
        << "      " << COPYRIGHT << endl;

   m_pSimulation->LogStream << endl
        << PROGRAM_NAME << endl
        << "      for " << PLATFORM << " (" << strGetBuild() << " build)" << endl
        << "      " << COPYRIGHT << endl;
}

/**
 * @brief Initialize CPU and wall-clock timers
 * 
 * Sets up two timing mechanisms:
 * 1. CPU time: clock() for processor time (excludes I/O waits)
 * 2. Wall time: time() for actual elapsed time
 * 
 * @param m_pSimulation Pointer to simulation (stores start time)
 * 
 * @note CPU time can rollover on very long runs (handled via m_dClkLast)
 * @warning If clock() fails, CPU timing is disabled (m_dCPUClock = -1)
 */
void CScreenPresenter::StartClock(CSimulation* m_pSimulation)
{
   // Initialize CPU clock
   if (static_cast<clock_t>(-1) == clock())
   {
      // There's a problem with the clock, but continue anyway
      m_pSimulation->LogStream << NOTE << "CPU time not available" << endl;
      m_dCPUClock = -1;
   }
   else
   {
      // All OK, so get the time in m_dClkLast (this is needed to check for clock rollover on long runs)
      m_dClkLast = static_cast<double>(clock());
      m_dClkLast -= CLOCK_T_MIN; // necessary if clock_t is signed to make m_dClkLast unsigned
   }

   // And now get the actual time we started
   m_pSimulation->m_tSysStartTime = time(nullptr);
}


/**
 * @brief Find directory containing Barrier executable
 * 
 * Extracts directory path from argv[0] by removing executable name.
 * 
 * @param pcArg argv[0] from main() (full path to executable)
 * @return true if directory found, false if path empty
 * 
 * @note Result stored in m_strSVDir (terminated with PATH_SEPARATOR)
 * @warning Does not validate if directory actually exists
 */
bool CScreenPresenter::bFindExeDir(char const* pcArg)
{
   string strTmp;

    // Use argv[0] as path
    strTmp = pcArg;


   // Neither approach has worked, so give up
   if (strTmp.empty())
      return false;

   // It's OK, so trim off the executable's name
   int nPos = static_cast<int>(strTmp.find_last_of(PATH_SEPARATOR));
   m_strSVDir = strTmp.substr(0, nPos + 1); // Note that this must be terminated with a backslash

   return true;
}

/**
 * @brief Get hostname of computer running simulation
 * 
 * Platform-specific implementations:
 * - Windows: GetComputerNameA() from windows.h
 * - Linux/Unix: uname() from sys/utsname.h
 * 
 * @return Hostname string, or empty string if detection fails
 * 
 * @note Used in simulation start banner for identifying compute resources
 */
string CScreenPresenter::strGetComputerName()
{  
    string strComputerName;

#ifdef _WIN32
    char buffer[MAX_COMPUTERNAME_LENGTH + 1];
    DWORD size = sizeof(buffer);
    if (GetComputerNameA(buffer, &size))
    {
        strComputerName = buffer;
    }
#else
    struct utsname uts;
    if (uname(&uts) == 0)
    {
        strComputerName = uts.nodename;
    }
#endif

    return strComputerName;
}

/**
 * @brief Display copyright, disclaimer, and simulation start info
 * 
 * Prints:
 * - Copyright notice
 * - 6-line disclaimer (no warranty, use at own risk)
 * - Hostname and start timestamp
 * - "Initializing..." message
 * 
 * @param m_pSimulation Pointer to simulation (for start time)
 * 
 * @note Required by GPL license to inform users of warranty status
 */
void CScreenPresenter::AnnounceLicence(const CSimulation* m_pSimulation) {
   cout << COPYRIGHT << endl
        << endl;
   cout << LINE << endl;
   cout << DISCLAIMER1 << endl;
   cout << DISCLAIMER2 << endl;
   cout << DISCLAIMER3 << endl;
   cout << DISCLAIMER4 << endl;
   cout << DISCLAIMER5 << endl;
   cout << DISCLAIMER6 << endl;
   cout << LINE << endl
        << endl;

   cout << START_NOTICE << strGetComputerName() << " at " << put_time(localtime(&m_pSimulation->m_tSysStartTime), "%H:%M on %A %d %B %Y") << endl;
   cout << INITIALIZING_NOTICE << endl;
}

/**
 * @brief Print final status message
 * 
 * @param strText Message to display (e.g., error description, success message)
 */
void CScreenPresenter::AnnounceEnding(const string& strText)
{
   cout << endl
        << strText << endl;
}