/*!
 *
 * \file read_input.cpp
 * \brief Reads non-GIS input files
 * \details TODO 001 A more detailed description of these routines.
 * \author Manuel Cobos Budia

 * \date 15/08/2024
 * \copyright GNU General Public License
 *
 */

/*==============================================================================================================================

==============================================================================================================================*/
#include <cstdlib> // for strtod()
#include <fstream>
using std::ifstream;

#include <iostream>
using std::cerr;
using std::cout;
using std::endl;
using std::ios;

#include <string>
using std::to_string;

#include <algorithm>
using std::find;

#include "main.h"
#include "simulation.h"

bool CSimulation::bReadConfigurationFile(void)
{
    // Create an ifstream object
    ifstream InStream;

    // Try to open run details file for input
    InStream.open(SV_INI.c_str(), ios::in);

    // Did it open OK?
    if (!InStream.is_open())
    {
        // Error: cannot open run details file for input
        cerr << ERR << "cannot open " << SV_INI << " for input" << endl;
        return false;
    }

    int nLine = 0;
    int i = 0;
    size_t nPos;
    string strRec, strErr;

    while (getline(InStream, strRec))
    {
        nLine++;

        // Trim off leading and trailing whitespace
        strRec = strTrim(&strRec);

        // If it is a blank line or a comment then ignore it
        if ((! strRec.empty()) && (strRec[0] != QUOTE1) && (strRec[0] != QUOTE2))
        {
            // It isn't so increment counter
            i++;

            // Find the colon: note that lines MUST have a colon separating data from leading description portion
            nPos = strRec.find(COLON);
            if (nPos == string::npos)
            {
                // Error: badly formatted (no colon)
                cerr << ERR << "on line " << to_string(nLine) << "badly formatted (no ':') in " << m_strDataPathName << endl
                     << strRec << endl;
                return false;
            }

            // Strip off leading portion (the bit up to and including the colon)
            string strRH = strRec.erase(0, nPos+1);

            // Remove leading whitespace after the colon
            strRH = strTrimLeft(&strRH);

            // Look for trailing comments, if found then terminate string at that point and trim off any trailing whitespace
            bool bFound = true;
            while (bFound)
            {
                bFound = false;

                nPos = strRH.rfind(QUOTE1);
                if (nPos != string::npos)
                {
                    strRH.resize(nPos);
                    bFound = true;
                }

                nPos = strRH.rfind(QUOTE2);
                if (nPos != string::npos)
                {
                    strRH.resize(nPos);
                    bFound = true;
                }

                // Trim trailing spaces
                strRH = strTrimRight(&strRH);
            }

            bool bFirst = true;
            int nRet = 0,
             nHour = 0,
             nMin = 0,
             nSec = 0,
             nDay = 0,
             nMonth = 0,
             nYear = 0;
            double dMult = 0;
            string strTmp;
            vector<string> VstrTmp;

            switch (i)
            {
                // ---------------------------------------------- Run Information -----------------------------------------------------
                case 1:
                // Text output file names, don't change case
                if (strRH.empty())
                    strErr = "line " + to_string(nLine) + ": output file names";
                else
                {
                    m_strRunName = strRH;

                    m_strOutFile = m_strOutPath;
                    m_strOutFile.append(strRH);
                    m_strOutFile.append(OUTEXT);

                    m_strLogFile = m_strOutPath;
                    m_strLogFile.append(strRH);
                    m_strLogFile.append(LOGEXT);
                }
                break;

            case 2:
                // Content of log file, 0 = no log file, 1 = least detail, 3 = most detail
            	if (! bIsStringValidInt(strRH))
            	{
               		strErr = "line " + to_string(nLine) + ": invalid integer for log file detail level '" + strRH + "' in " + m_strDataPathName;
               		break;
            	}

            	m_nLogFileDetail = stoi(strRH);

            	if ((m_nLogFileDetail < NO_LOG_FILE) || (m_nLogFileDetail > LOG_FILE_HIGH_DETAIL))
               		strErr = "line " + to_string(nLine) + ": log file detail level";
            	break;

         	case 3:
            	// Get the start date/time of the simulation, format is [hh-mm-ss dd/mm/yyyy]
            	VstrTmp = VstrSplit(&strRH, SPACE);

            	// Both date and time here?
            	if (VstrTmp.size() < 2)
            	{
               		strErr = "line " + to_string(nLine) + ": must have both date and time for simulation start in '" + m_strDataPathName + "'";
               		break;
            	}

            	// OK, first sort out the time
            	if (! bParseTime(&VstrTmp[0], nHour, nMin, nSec))
            	{
               		strErr = "line " + to_string(nLine) + ": could not understand simulation start time in '" + m_strDataPathName + "'";
               		break;
            	}

            	// Next sort out the date
            	if (! bParseDate(&VstrTmp[1], nDay, nMonth, nYear))
            	{
               		strErr = "line " + to_string(nLine) + ": could not understand simulation start date in '" + m_strDataPathName + "'";
               		break;
            	}

            	// Store simulation start time and date
            	m_nSimStartSec = nSec;
            	m_nSimStartMin = nMin;
            	m_nSimStartHour = nHour;
            	m_nSimStartDay = nDay;
            	m_nSimStartMonth = nMonth;
            	m_nSimStartYear = nYear;

         	case 6:
            	// Get the along channel geometry file name
				if (strRH.empty())
               		strErr = "line " + to_string(nLine) + ": along channel geometry file name";
            	else
            	{
               		m_strAlongChannelGeometryFilename = strRH;
				}
			case 7:
            	// Get the cross section channel geometry file name
				if (strRH.empty())
               		strErr = "line " + to_string(nLine) + ": cross sections file name";
            	else
            	{
               		m_strCrossSectionsFilename = strRH;
				}

			case 8:
            	// Get the hydro file name
				if (strRH.empty())
               		strErr = "line " + to_string(nLine) + ": hydro file name";
            	else
            	{
               		m_strHydroFilename = strRH;
				}
			}
    	}
	}
	// Close file
    InStream.close();

    // Finally, need to check that we have at least one raster file, so that we know the grid size and units (and preferably also the projection)
    // bool bNoRasterFiles = true;
	return true;
}
