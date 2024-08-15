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
#include <assert.h>

#ifdef _WIN32
#include <windows.h>       // Needed for CalcProcessStats()
#include <psapi.h>
#include <io.h>            // For isatty()
#elif defined __GNUG__
#include <sys/resource.h>  // Needed for CalcProcessStats()
#include <unistd.h>        // For isatty()
#endif

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
using std::noshowpos;
using std::showpos;

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

#include <numeric>
using std::accumulate;
using std::inner_product;


#include "main.h"
#include "simulation.h"

//===============================================================================================================================
//! Tells the user that we have started the simulation
//===============================================================================================================================
void CSimulation::AnnounceStart(void)
{
   cout << endl
        << PROGRAM_NAME << " for " << PLATFORM << " " << strGetBuild() << endl;
}

//===============================================================================================================================
//! Starts the clock ticking
//===============================================================================================================================
void CSimulation::StartClock(void)
{
   // First start the 'CPU time' clock ticking
   if (static_cast<clock_t>(-1) == clock())
   {
      // There's a problem with the clock, but continue anyway
      LogStream << NOTE << "CPU time not available" << endl;
      m_dCPUClock = -1;
   }
   else
   {
      // All OK, so get the time in m_dClkLast (this is needed to check for clock rollover on long runs)
      m_dClkLast = static_cast<double>(clock());
      m_dClkLast -= CLOCK_T_MIN; // necessary if clock_t is signed to make m_dClkLast unsigned
   }

   // And now get the actual time we started
   m_tSysStartTime = time(nullptr);
}

//===============================================================================================================================
//! Finds the folder (directory) in which the CoastalME executable is located
//===============================================================================================================================
bool CSimulation::bFindExeDir(char const* pcArg)
{
   string strTmp;
   char szBuf[BUF_SIZE] = "";

#ifdef _WIN32
   if (0 != GetModuleFileName(NULL, szBuf, BUF_SIZE))
      strTmp = szBuf;
   else
      // It failed, so try another approach
      strTmp = pcArg;
#else
   //   char* pResult = getcwd(szBuf, BUF_SIZE);          // Used to use this, but what if cwd is not the same as the CoastalME dir?

   if (-1 != readlink("/proc/self/exe", szBuf, BUF_SIZE))
      strTmp = szBuf;
   else
      // It failed, so try another approach
      strTmp = pcArg;
#endif

   // Neither approach has worked, so give up
   if (strTmp.empty())
      return false;

   // It's OK, so trim off the executable's name
   int nPos = static_cast<int>(strTmp.find_last_of(PATH_SEPARATOR));
   m_strCMEDir = strTmp.substr(0, nPos + 1); // Note that this must be terminated with a backslash

   return true;
}

//===============================================================================================================================
//! Tells the user about the licence
//===============================================================================================================================
void CSimulation::AnnounceLicence(void)
{
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

   cout << START_NOTICE << strGetComputerName() << " at " << put_time(localtime(&m_tSysStartTime), "%T on %A %d %B %Y") << endl;
   cout << INITIALIZING_NOTICE << endl;
}
