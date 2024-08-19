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

#include "main.h"
#include "simulation.h"


//===============================================================================================================================
//! Returns the date and time on which the program was compiled
//===============================================================================================================================
string CSimulation::strGetBuild(void)
{
   string strBuild("(");
   strBuild.append(__TIME__);
   strBuild.append(" ");
   strBuild.append(__DATE__);
#ifdef _DEBUG
   strBuild.append(" DEBUG");
#endif
   strBuild.append(" build)");

   return strBuild;
}


//===============================================================================================================================
//! From http://stackoverflow.com/questions/236129/split-a-string-in-c They implement (approximately) Python's split() function. This first version puts the results into a pre-constructed string vector. It ignores empty items
//===============================================================================================================================
vector<string> *CSimulation::VstrSplit(string const *s, char const delim, vector<string> *elems)
{
   stringstream ss(*s);
   string item;
   while (getline(ss, item, delim))
   {
      if (!item.empty())
         elems->push_back(item);
   }
   return elems;
}

//===============================================================================================================================
//! From http://stackoverflow.com/questions/236129/split-a-string-in-c They implement (approximately) Python's split() function. This second version returns a new string vector (it calls the first version)
//===============================================================================================================================
vector<string> CSimulation::VstrSplit(string const *s, char const delim)
{
   vector<string> elems;
   VstrSplit(s, delim, &elems);
   return elems;
}

//===============================================================================================================================
//! Opens the log file
//===============================================================================================================================
bool CSimulation::bOpenLogFile(void)
{
   // Open in binary mode if just checking random numbers
#ifdef RANDCHECK
   if (m_nLogFileDetail == 0)
   {
      LogStream.open("/dev/null", ios::out | ios::binary | ios::trunc);
      cout << "Warning: log file is not writting" << endl;
   }
   else
      LogStream.open(m_strLogFile.c_str(), ios::out | ios::binary | ios::trunc);

#else
   if (m_nLogFileDetail == 0)
   {
      LogStream.open("/dev/null", ios::out | ios::trunc);
      cout << "Warning: log file is not writting" << endl;
   }
   else
      LogStream.open(m_strLogFile.c_str(), ios::out | ios::trunc);
#endif

   if (!LogStream)
   {
      // Error, cannot open log file
      cerr << ERR << "cannot open " << m_strLogFile << " for output" << endl;
      return false;
   }

   return true;
}
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

   if (0 != GetModuleFileName(NULL, szBuf, BUF_SIZE))
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

   // cout << START_NOTICE << strGetComputerName() << " at " << put_time(localtime(&m_tSysStartTime), "%T on %A %d %B %Y") << endl;
   cout << INITIALIZING_NOTICE << endl;
}


//===============================================================================================================================
//! Returns an error message given an error code
//===============================================================================================================================
string CSimulation::strGetErrorText(int const nErr)
{
   string strErr;

   switch (nErr)
   {
   // case RTN_USER_ABORT:
   //    strErr = "aborted by user";
   //    break;
   // case RTN_ERR_BADPARAM:
   //    strErr = "error in command-line parameter";
   //    break;
   case RTN_ERR_INI:
      strErr = "error reading initialization file";
      break;
   case RTN_ERR_SVDIR:
      strErr = "error in directory name";
      break;
   default:
      // should never get here
         strErr = "unknown error";
   }

   return strErr;
}


//===============================================================================================================================
//! Trims whitespace from the left side of a string, does not change the original string
//===============================================================================================================================
string CSimulation::strTrimLeft(string const *strIn)
{
   // Trim leading spaces
   size_t nStartpos = strIn->find_first_not_of(" \t");
   if (nStartpos == string::npos)
      return *strIn;
   else
      return strIn->substr(nStartpos);
}

//===============================================================================================================================
//! Trims whitespace from the right side of a string, does not change the original string
//===============================================================================================================================
string CSimulation::strTrimRight(string const *strIn)
{
   string strTmp(*strIn);

   // Remove any stray carriage returns (can happen if file was edited in Windows)
   strTmp.erase(std::remove(strTmp.begin(), strTmp.end(), '\r'), strTmp.end());

   // Trim trailing spaces
   size_t nEndpos = strTmp.find_last_not_of(" \t");
   if (nEndpos == string::npos)
      return strTmp;
   else
      return strTmp.substr(0, nEndpos + 1);
}

//===============================================================================================================================
//! Trims whitespace from both sides of a string, does not change the original string
//===============================================================================================================================
string CSimulation::strTrim(string const *strIn)
{
   string strTmp = *strIn;

   // Remove any stray carriage returns (can happen if file was edited in Windows)
   strTmp.erase(std::remove(strTmp.begin(), strTmp.end(), '\r'), strTmp.end());

   // Trim trailing spaces
   size_t nPos = strTmp.find_last_not_of(" \t");

   if (nPos != string::npos)
      strTmp.resize(nPos+1);

   // Trim leading spaces
   nPos = strTmp.find_first_not_of(" \t");

   if (nPos != string::npos)
      strTmp = strTmp.substr(nPos);

   return strTmp;
}


//===============================================================================================================================
//! Parses a date string into days, months, and years, and checks each of them
//===============================================================================================================================
bool CSimulation::bParseDate(string const *strDate, int &nDay, int &nMonth, int &nYear)
{
#ifdef _WIN32
   // For Windows, make sure has backslashes, not Unix-style slashes
   vector<string> VstrTmp = VstrSplit(strDate, SLASH);
#else
   vector<string> VstrTmp = VstrSplit(strDate, SLASH);
#endif

   if (VstrTmp.size() < 3)
   {
      cerr << "date string must include day, month, and year '" << strDate << "'" << endl;
      return false;
   }

   // Sort out day
   if (! bIsStringValidInt(VstrTmp[0]))
   {
      cerr << "invalid integer for day in date '" << strDate << "'" << endl;
      return false;
   }

   nDay = stoi(VstrTmp[0]);

   if ((nDay < 1) || (nDay > 31))
   {
      cerr << "day must be between 1 and 31 in date '" << strDate << "'" << endl;
      return false;
   }

   // Sort out month
   if (! bIsStringValidInt(VstrTmp[1]))
   {
      cerr << "invalid integer for month in date '" << strDate << "'" << endl;
      return false;
   }

   nMonth = stoi(VstrTmp[1]);

   if ((nMonth < 1) || (nMonth > 12))
   {
      cerr << "month must be between 1 and 12 in date '" << strDate << "'" << endl;
      return false;
   }

   // Sort out year
   if (! bIsStringValidInt(VstrTmp[2]))
   {
      cerr << "invalid integer for year in date '" << strDate << "'" << endl;
      return false;
   }

   nYear = stoi(VstrTmp[2]);

   if (nYear < 0)
   {
      cerr << "year must be > 0 in date '" << strDate << "'" << endl;
      return false;
   }

   return true;
}

//===============================================================================================================================
//! Parses a time string into hours, minutes, and seconds, and checks each of them
//===============================================================================================================================
bool CSimulation::bParseTime(string const *strTime, int &nHour, int &nMin, int &nSec)
{
   vector<string> VstrTmp = VstrSplit(strTime, DASH);

   if (VstrTmp.size() < 3)
   {
      cerr << "time string must include hours, minutes, and seconds '" << strTime << "'" << endl;
      return false;
   }

   // Sort out hour
   if (! bIsStringValidInt(VstrTmp[0]))
   {
      cerr << "invalid integer for hours in time '" << strTime << "'" << endl;
      return false;
   }

   nHour = stoi(VstrTmp[0]);

   if ((nHour < 0) || (nHour > 23))
   {
      cerr << "hour must be between 0 and 23 in time '" << strTime << "'" << endl;
      return false;
   }

   // Sort out minutes
   if (! bIsStringValidInt(VstrTmp[1]))
   {
      cerr << "invalid integer for minutes in time '" << strTime << "'" << endl;
      return false;
   }

   nMin = stoi(VstrTmp[1]);

   if ((nMin < 0) || (nMin > 59))
   {
      cerr << "minutes must be betwen 0 and 59 in time '" << strTime << "'" << endl;
      return false;
   }

   // Sort out seconds
   if (! bIsStringValidInt(VstrTmp[2]))
   {
      cerr << "invalid integer for seconds in time '" << strTime << "'" << endl;
      return false;
   }

   nSec = stoi(VstrTmp[2]);

   if ((nSec < 0) || (nSec > 59))
   {
      cerr << "seconds must be between 0 and 59 in time '" << strTime << "'" << endl;
      return false;
   }

   return true;
}


//===============================================================================================================================
//! Checks to see if a string can be read as a valid double number. Does not find trailing (i.e.post-number) rubbish, but then neither does strtod(). From https://stackoverflow.com/questions/392981/how-can-i-convert-string-to-double-in-c
//===============================================================================================================================
bool bIsStringValidDouble(string& str)
{
   std::istringstream iStr(str);
   double dDummy;

   if (!(iStr >> dDummy))
      return false;

   return true;
}

//===============================================================================================================================
//! Checks to see if a string can be read as a valid integer, from https://stackoverflow.com/questions/2844817/how-do-i-check-if-a-c-string-is-an-int
//===============================================================================================================================
bool bIsStringValidInt(string& str)
{
   // Trim leading whitespace
   size_t nPos = str.find_first_not_of(" \t");
   if (nPos != string::npos)
      str = str.substr(nPos);

   // If the first character is the sign, remove it
   if ((str[0] == '-') || (str[0] == '+'))
      str.erase(0, 1);

   // Now check that the string contains only numbers
   return (str.find_first_not_of("0123456789") == string::npos);
}