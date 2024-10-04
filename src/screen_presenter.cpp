/*!
 *
 * \file utils.cpp
 * \brief Utility routines
 * \details TODO 001 A more detailed description of this routine.
 * \author Manuel Cobos Budia

 * \date 2024
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

#include <iostream>
using std::cerr;
using std::cout;
using std::endl;
using std::ios;

// #include <psapi.h>
#include <io.h>            // For isatty()

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

#include <windows.h>       // Needed for GetModuleFileName()

#include "screen_presenter.h"
#include "utils.h"
#include "simulation.h"
#include "main.h"



//===============================================================================================================================
//! The CScreenPresenter constructor
//===============================================================================================================================
CScreenPresenter::CScreenPresenter() {
   m_tSysEndTime = 0;
}

//===============================================================================================================================
//! The CScreenPresenter destructor
//===============================================================================================================================
CScreenPresenter::~CScreenPresenter() = default;


//===============================================================================================================================
//! The nDoSimulation member function of CSimulation sets up and runs the simulation
//===============================================================================================================================
void CScreenPresenter::StartingRun(int nArg, char const* pcArgv[], CSimulation* m_pSimulation)
{
   // ================================================== initialization section ================================================
   // Hello, World!
   AnnounceStart(m_pSimulation);

   // Start the clock ticking
   StartClock(m_pSimulation);

   //    return (RTN_ERR_SVDIR);
   bFindExeDir(pcArgv[0]);

   // OK, we are off, tell the user about the licence and the start time
   AnnounceLicence(m_pSimulation);
}



//===============================================================================================================================
//! Notifies the user that the simulation has ended, asks for keypress if necessary, and if compiled under GNU can send an email
//===============================================================================================================================
void CScreenPresenter::EndingRun(CSimulation* pSimulation)
{
   // If we don't know the time that the run ended (e.g. because it did not finish correctly), then get it now
   static double sdElapsed = 0;
   m_tSysEndTime = time(nullptr);

   // Calculate time elapsed and remaining
   sdElapsed = difftime(m_tSysEndTime, pSimulation->m_tSysStartTime);

   // if (nRtn  == RTN_OK)
   //    // normal ending
   // {
   cout << "\r    - Elapsed Time: " << std::fixed << setprecision(3) << setw(6) << sdElapsed  << endl;
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

//===============================================================================================================================
//! Returns the date and time on which the program was compiled
//===============================================================================================================================
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


//===============================================================================================================================
//! Tells the user that we have started the simulation
//===============================================================================================================================
void CScreenPresenter::AnnounceStart(CSimulation* m_pSimulation)
{
   cout << endl
        << PROGRAM_NAME << " for " << PLATFORM << " " << strGetBuild() << endl;

   m_pSimulation->LogStream << endl
        << PROGRAM_NAME << " for " << PLATFORM << " " << strGetBuild() << endl;
}

//===============================================================================================================================
//! Starts the clock ticking
//===============================================================================================================================
void CScreenPresenter::StartClock(CSimulation* m_pSimulation)
{
   // First start the 'CPU time' clock ticking
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


//===============================================================================================================================
//! Finds the folder (directory) in which the CoastalME executable is located
//===============================================================================================================================
bool CScreenPresenter::bFindExeDir(char const* pcArg)
{
   string strTmp;
   char szBuf[BUF_SIZE] = "";

   if (0 != GetModuleFileName(nullptr, szBuf, BUF_SIZE))
      strTmp = szBuf;
   else
      // It failed, so try another approach
      strTmp = pcArg;


   // Neither approach has worked, so give up
   if (strTmp.empty())
      return false;

   // It's OK, so trim off the executable's name
   int nPos = static_cast<int>(strTmp.find_last_of(PATH_SEPARATOR));
   m_strSVDir = strTmp.substr(0, nPos + 1); // Note that this must be terminated with a backslash

   return true;
}

//===============================================================================================================================
//! Returns a string, hopefully giving the name of the computer on which the simulation is running
//===============================================================================================================================
string CScreenPresenter::strGetComputerName()
{
   char* strComputerName;
   // Being compiled to run under Windows, either by MS VC++, Borland C++, or Cygwin
   strComputerName = getenv("COMPUTERNAME");

   return strComputerName;
}

//===============================================================================================================================
//! Tells the user about the licence
//===============================================================================================================================
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


void CScreenPresenter::AnnounceEnding(const string& strText)
{
   cout << endl
        << strText << endl;
}