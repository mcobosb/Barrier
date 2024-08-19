
#include <cassert>

#include <iostream>
using std::cerr;
using std::cout;
using std::endl;
using std::ios;

#include <iomanip>
using std::setprecision;

#include <string>
using std::to_string;

#include "main.h"
#include "simulation.h"

//===============================================================================================================================
//! The CSimulation constructor
//===============================================================================================================================
CSimulation::CSimulation (void)
{

   m_ulIter = 0;
   //! The CME folder
   string m_strSVDir;

   //! Folder for the CME .ini file
   string m_strSVIni;

   //! An email addresx to which to send end-of-simulation messages
   string m_strMailAddress;

   //! Name of along channel geometry file
   string m_strAlongChannelGeometryFilename;
   //! Name of cross sections file
   string m_strCrossSectionsFilename;
   //! Name of hydro file
   string m_strHydroFilename;
}

//===============================================================================================================================
//! The CSimulation destructor
//===============================================================================================================================
CSimulation::~CSimulation (void)
{
   // Close output files if open
   // if (LogStream && LogStream.is_open())
   // {
   //    LogStream.flush();
   //    LogStream.close();
   // }

   // if (OutStream && OutStream.is_open())
   // {
   //    OutStream.flush();
   //    OutStream.close();
   // }
}


//===============================================================================================================================
//! The nDoSimulation member function of CSimulation sets up and runs the simulation
//===============================================================================================================================
int CSimulation::nDoSimulation(int nArg, char const* pcArgv[])
{
    // ================================================== initialization section ================================================
    // Hello, World!
    AnnounceStart();

    // Start the clock ticking
    StartClock();

    // Find out the folder in which the SV executable sits, in order to open the .ini file (they are assumed to be in the same folder)
    if (! bFindExeDir(pcArgv[0]))
       return (RTN_ERR_SVDIR);

    // Deal with command-line parameters
    // int nRet = nHandleCommandLineParams (nArg, pcArgv);
    // if (nRet != RTN_OK)
    //    return (nRet);

    // OK, we are off, tell the user about the licence and the start time
    AnnounceLicence();

    // Read the .ini file and get the name of the run-data file, and path for output etc.
    if (! bReadConfigurationFile())
       return (RTN_ERR_INI);

   // Open log file
   if (! bOpenLogFile())
      return (RTN_ERR_LOGFILE);
    return RTN_OK;
}


//===============================================================================================================================
//! Notifies the user that the simulation has ended, asks for keypress if necessary, and if compiled under GNU can send an email
//===============================================================================================================================
void CSimulation::DoSimulationEnd(int const nRtn)
{
   // If we don't know the time that the run ended (e.g. because it did not finish correctly), then get it now
   if (m_tSysEndTime == 0)
      m_tSysEndTime = time(nullptr);

   switch (nRtn)
   {
   case (RTN_OK):
      // normal ending
      cout << RUN_END_NOTICE << std::put_time(localtime(&m_tSysEndTime), "%T %A %d %B %Y") << endl;
      break;

   default:
      // Aborting because of some error
      cerr << RUN_END_NOTICE << "iteration " << m_ulIter << ERROR_NOTICE << nRtn << " (" << strGetErrorText(nRtn) << ") on " << std::put_time(localtime(&m_tSysEndTime), "%T %A %d %B %Y") << endl;

      // Write the error message to the logfile and to stdout
      if (LogStream && LogStream.is_open())
      {
         LogStream << ERR << strGetErrorText(nRtn) << " (error code " << nRtn << ") on " << std::put_time(localtime(&m_tSysEndTime), "%T %A %d %B %Y") << endl;
         LogStream.flush();
      }

      if (OutStream && OutStream.is_open())
      {
         OutStream << ERR << strGetErrorText(nRtn) << " (error code " << nRtn << ") on " << std::put_time(localtime(&m_tSysEndTime), "%T %A %d %B %Y") << endl;
         OutStream.flush();
      }
   }
}