/*!
 *
 * \file data_reader.cpp
 * \brief Reads non-GIS input files
 * \details TODO 001 A more detailed description of these routines.
 * \author Manuel Cobos Budia

 * \date 15/08/2024
 * \copyright GNU General Public License
 *
 */

/*==============================================================================================================================

==============================================================================================================================*/
#include <cstdlib> // for strtol() and strtod()
#include <data_reader.h>
#include <fstream>
using std::ifstream;

#include <sstream>
using std::stringstream;

#include <iostream>
using std::cerr;
using std::cout;
using std::endl;
using std::ios;

#include <string>
using std::to_string;
using std::stoi;
using std::string;

#include <vector>
using std::vector;

#include <algorithm>
using std::find;

#include "simulation.h"
#include "data_reader.h"
#include "error_handling.h"
#include "utils.h"

//===============================================================================================================================
//! The CDataReader constructor
//===============================================================================================================================
CDataReader::CDataReader() {
	// CEstuary m_pEstuary;
	// m_pSimulation = new CSimulation;
	m_ulIter = 0;

	//! Date properties
	m_nSimStartSec =
	m_nSimStartMin =
	m_nSimStartHour =
	m_nSimStartDay =
	m_nSimStartMonth =
	m_nSimStartYear =
	m_nLogFileDetail = 0;

	m_dDurationUnitsMultiplier = 0.0;



}

//======================================================================================================================
//! The CDataReader destructor
//======================================================================================================================
CDataReader::~CDataReader() = default;


//===============================================================================================================================
//! Opens the log file
//===============================================================================================================================
bool CDataReader::bOpenLogFile(CSimulation* m_pSimulation)
{
	if (m_nLogFileDetail != 1)
	{
		m_pSimulation->LogStream.open("/dev/null", ios::out | ios::trunc);

	}
	else
		m_pSimulation->LogStream.open(m_pSimulation->m_strLogFile.c_str(), ios::out | ios::trunc);


	if (!m_pSimulation->LogStream)
	{
		cout << "Warning: log file is not writting" << endl;
	}

	return true;
}

//======================================================================================================================
//! Read configuration.ini file
//======================================================================================================================
bool CDataReader::bReadConfigurationFile(CSimulation* m_pSimulation)
{
    // Create an ifstream object
    ifstream InStream;

    // Try to open run details file for input
    InStream.open(SV_INI.c_str(), ios::in);

	// Create the raster grid object
	// m_pSimulation = new CSimulation;

    // Did it open OK?
    if (!InStream.is_open())
    {
        // Error: cannot open run details file for input
        cerr << ERR << "cannot open " << SV_INI << " for input" << endl;
        return true;
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
                cerr << ERR << "on line " << to_string(nLine) << "badly formatted (no ':') in " << SV_INI << endl
                << strRec << endl;
            	return true;
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

	                    m_pSimulation->m_strOutFile = m_strOutPath;
	                    m_pSimulation->m_strOutFile.append(strRH);
	                    m_pSimulation->m_strOutFile.append(OUT_EXT);

	                    m_pSimulation->m_strLogFile = m_strOutPath;
	                    m_pSimulation->m_strLogFile.append(strRH);
	                    m_pSimulation->m_strLogFile.append(LOG_EXT);
	                }
	                break;

				case 2:
	                // Content of log file, 0 = no log file, 1 = least detail, 3 = most detail
            		if (! bIsStringValidInt(strRH))
            		{
               			strErr = "line " + to_string(nLine) + ": invalid integer for log file detail level '" + strRH + "' in " + SV_INI;
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
               			strErr = "line " + to_string(nLine) + ": must have both date and time for simulation start in '" + SV_INI + "'";
               			break;
            		}

            		// OK, first sort out the time
            		if (! bParseTime(&VstrTmp[0], nHour, nMin, nSec))
            		{
               			strErr = "line " + to_string(nLine) + ": could not understand simulation start time in '" + SV_INI + "'";
               			break;
            		}

            		// Next sort out the date
            		if (! bParseDate(&VstrTmp[1], nDay, nMonth, nYear))
            		{
               			strErr = "line " + to_string(nLine) + ": could not understand simulation start date in '" + SV_INI + "'";
               			break;
            		}

            		// Store simulation start time and date
            		m_nSimStartSec = nSec;
            		m_nSimStartMin = nMin;
            		m_nSimStartHour = nHour;
            		m_nSimStartDay = nDay;
            		m_nSimStartMonth = nMonth;
            		m_nSimStartYear = nYear;
            		break;

	            case 4: {
		            // Duration of simulation (in hours, days, months, or years): sort out multiplier and user units, as used in the per-timestep output
            		strRH = strToLower(&strRH);

            		nRet = nSimulationTimeMultiplier(&strRH);

	            	//! TODO 020: Options for another time unit
	            	m_pSimulation->m_dTimeFactor = 1.0;

            		if (nRet != RTN_OK)
            		{
            			strErr = "line " + to_string(nLine) + ": units for duration of simulation";
            			break;
            		}

            		// And now calculate the duration of the simulation in hours: first find whitespace between the number and the unit
            		nPos = strRH.rfind(SPACE);
            		if (nPos == string::npos)
            		{
            			strErr = "line " + to_string(nLine) + ": format of duration simulation line";
            			break;
            		}

            		// Cut off rh bit of string
            		strRH.resize(nPos);

            		// Remove trailing spaces
            		strRH = strTrimRight(&strRH);

            		// Calculate the duration of the simulation in hours
            		double dSimDuration = strtod(strRH.c_str(), nullptr) * m_dDurationUnitsMultiplier*m_pSimulation->m_dTimeFactor;

            		if (dSimDuration <= 0)
            			strErr = "line " + to_string(nLine) + ": duration of simulation must be > 0";
            		else
            			m_pSimulation->dSetSimulationDuration(dSimDuration);

            		break;
	            }

				case 5: {
		            // Timestep of simulation (in hours or days)
            		strRH = strToLower(&strRH);

            		dMult = dGetTimeMultiplier(&strRH);
	            	//! TODO 020: Options for another time unit
	            	m_pSimulation->m_dTimeFactor = 1.0;

            		if (static_cast<int>(dMult) == TIME_UNKNOWN)
            		{
            			strErr = "line " + to_string(nLine) + ": units for simulation timestep";
            			break;
            		}

            		// we have the multiplier, now calculate the timestep in hours: look for the whitespace between the number and unit
            		nPos = strRH.rfind(SPACE);
            		if (nPos == string::npos)
            		{
            			strErr = "line " + to_string(nLine) + ": format of simulation timestep";
            			break;
            		}

            		// cut off rh bit of string
            		strRH.resize(nPos);

            		// remove trailing spaces
            		strRH = strTrimRight(&strRH);

            		// Check that this is a valid double
            		if (! bIsStringValidDouble(strRH))
            		{
            			strErr = "line " + to_string(nLine) + ": invalid floating point number for timestep '" + strRH + "' in " + SV_INI;
            			break;
            		}

            		double dTimeStep = strtod(strRH.c_str(), nullptr) * dMult*m_pSimulation->m_dTimeFactor; // in hours

            		if (dTimeStep <= 0)
            			strErr = "line " + to_string(nLine) + ": timestep of simulation must be > 0";

            		if (dTimeStep >= m_pSimulation->dGetSimulationDuration())
            			strErr = "line " + to_string(nLine) + ": timestep of simulation must be < the duration of the simulation";

            		m_pSimulation->dSetSimulationTimestep(dTimeStep);
            		break;
	            }

            	case 6: {
		            // Get the output variables
	            	if (strRH.empty())
	            		strErr = "line " + to_string(nLine) + ": along channel geometry file name";
	            	else {
	            		if (strRH == "full") {
	            			m_pSimulation->m_vOutputVariables = {"A", "Ap", "Ac", "Q", "Qp", "Qc", "Rh", "B", "eta", "beta", "I1", "I2", "rho", "U", "c", "S", "Qb", "Qs", "Qt", "xl", "xr"};
	            		}
	            		else {
	            			vector<string> vOutputVariables;
	            			// Obtain the new line
	            			stringstream sline(strRH);
	            			string token;

	            			// Using getline for spliting the string by commas
	            			while (getline(sline, token, ',')) {
	            				m_pSimulation->strAddOutputVariable(token);
	            			}
	            		}
	            	}
	            }

         		case 7: {
         			// Get the along channel geometry file name
         			if (strRH.empty())
         				strErr = "line " + to_string(nLine) + ": along channel geometry file name";
         			else
         			{
         				m_strAlongChannelDataFilename = strRH;
         				m_strAlongChannelDataFilename.append(".csv");
         			}
         			break;
         		}

				case 8: {
		            // Get the cross-sections channel geometry file name
            		if (strRH.empty())
            			strErr = "line " + to_string(nLine) + ": cross sections file name";
            		else
            		{
            			m_strCrossSectionsFilename = strRH;
            			m_strCrossSectionsFilename.append(".csv");
            		}
            		break;
	            }

	            case 9: {
		            // Get the initial along-channel estuarine condition
	            	if (strRH.empty())
	            		strErr = "line " + to_string(nLine) + ":  initial along-channel estuarine condition";
	            	else
	            	{
	            		// Convert string to int
	            		m_pSimulation->nSetInitialEstuarineCondition(strtol(strRH.c_str(), nullptr, 10));
	            	}
	            	break;
	            }

            	case 10: {
	            	// Get the upward estuarine boundary condition
	            	if (strRH.empty())
	            		strErr = "line " + to_string(nLine) + ": upward estuarine boundary condition";
	            	else
	            	{
	            		// Convert string to int
	            		m_pSimulation->nSetUpwardEstuarineCondition(strtol(strRH.c_str(), nullptr, 10));
	            	}
	            	break;
            	}

            	case 11: {
	            	// Get the UEBC filename

	            	// if (strRH.empty())
	            	// 	strErr = "line " + to_string(nLine) + ": tidal filename";

	            	if (m_pSimulation->nGetUpwardEstuarineCondition() == 1 || m_pSimulation->nGetUpwardEstuarineCondition() == 2)
	            	{
	            		m_pSimulation->m_strUpwardBoundaryConditionFilename = strRH;
	            		m_pSimulation->m_strUpwardBoundaryConditionFilename.append(".csv");
	            	}
	            	else {
	            		m_pSimulation->m_strUpwardBoundaryConditionFilename = "";
	            	}

	            	break;
            	}

            	case 12: {
	            	// Get the downward estuarine boundary condition
	            	if (strRH.empty())
	            		strErr = "line " + to_string(nLine) + ": downward estuarine boundary condition";
	            	else
	            	{
	            		int nEstuaryCondition = strtol(strRH.c_str(), nullptr, 10);
	            		if (nEstuaryCondition != 0 && nEstuaryCondition != 1 && nEstuaryCondition != 2) {
	            			//! TODO 007: return an error code
	            		}
	            		else {
	            			m_pSimulation->nSetDownwardEstuarineCondition(nEstuaryCondition);
	            		}

	            	}
	            	break;
            	}

            	case 13: {
		            // Get the tidal filename [if DEBC = E]
	            	// if (strRH.empty())
	            	// 	strErr = "line " + to_string(nLine) + ": - DEBC file name ";

	            	if (m_pSimulation->nGetDownwardEstuarineCondition() == 1 || m_pSimulation->nGetDownwardEstuarineCondition() == 2) {
	            		m_pSimulation->m_strDownwardBoundaryConditionFilename = strRH;
	            		m_pSimulation->m_strDownwardBoundaryConditionFilename.append(".csv");
		            }
	            	else {
	            		m_pSimulation->m_strDownwardBoundaryConditionFilename = "";
	            	}

	            	break;
            	}


				case 14:
            		// Get the hydro file name
					if (strRH.empty())
               			strErr = "line " + to_string(nLine) + ": hydro file name";
					else
            		{
               			m_strHydroFilename = strRH;
						m_strHydroFilename.append(".csv");
					}

            	case 15: {
            		// Get the courant number
            		if (strRH.empty())
            			strErr = "line " + to_string(nLine) + ": courant number";
            		else
            			m_pSimulation->dSetCourantNumber(strtod(strRH.c_str(), nullptr));
            		break;
            	}

            	case 16: {
					// Use McComarck limiter flux?
					strRH = strToLower(&strRH);

					if (strRH.empty())
						strErr = "line " + to_string(nLine) + ": Compute water density?";

					if  (strRH.find("y") != string::npos)
						m_pSimulation->bSetDoMcComarckLimiterFlux(true);
					else
						m_pSimulation->bSetDoMcComarckLimiterFlux(false);
					break;
            	}

            	case 17: {
					// Equation for the limiter flux
					if (strRH.empty())
						strErr = "line " + to_string(nLine) + ": equation for the limiter flux ";

					if (m_pSimulation->bGetDoMcComarckLimiterFlux())
					{
						m_pSimulation->nSetEquationLimiterFlux(strtol(strRH.c_str(), nullptr, 10));
					}
					break;
            	}

            	case 18: {
					// Psi formula
					if (strRH.empty())
						strErr = "line " + to_string(nLine) + ": psi formula ";

					if (m_pSimulation->bGetDoMcComarckLimiterFlux())
					{
						m_pSimulation->nSetPsiFormula(strtol(strRH.c_str(), nullptr, 10));
					}
					break;
            	}

            	case 19: {
					//! Delta value
					if (strRH.empty())
						strErr = "line " + to_string(nLine) + ": delta value ";

					if (m_pSimulation->bGetDoMcComarckLimiterFlux())
					{
						m_pSimulation->dSetDeltaValue(strtod(strRH.c_str(), nullptr));
					}
					break;
            	}

            	case 20: {
						// Use Surface Gradient method?
						strRH = strToLower(&strRH);

						if (strRH.empty())
							strErr = "line " + to_string(nLine) + ": surface gradient method";

						if  (strRH.find("y") != string::npos)
							m_pSimulation->bSetDoSurfaceGradientMethod(true);
						else
							m_pSimulation->bSetDoSurfaceGradientMethod(false);
						break;
            	}

            	case 21: {
						// Use Source Term balance?
						strRH = strToLower(&strRH);

						if (strRH.empty())
							strErr = "line " + to_string(nLine) + ": source Term balance";

						if  (strRH.find("y") != string::npos)
							m_pSimulation->bSetDoSurfaceTermBalance(true);
						else
							m_pSimulation->bSetDoSurfaceTermBalance(false);
						break;
            	}

            	case 22: {
						// Use beta coefficient?
						strRH = strToLower(&strRH);

						if (strRH.empty())
							strErr = "line " + to_string(nLine) + ": beta coefficient";

						if  (strRH.find("y") != string::npos)
							m_pSimulation->bSetDoBetaCoefficient(true);
						else
							m_pSimulation->bSetDoBetaCoefficient(false);
						break;
            	}

            	case 23: {
						// Use Dry bed?
						strRH = strToLower(&strRH);

						if (strRH.empty())
							strErr = "line " + to_string(nLine) + ": dry bed";

						if  (strRH.find("y") != string::npos)
							m_pSimulation->bSetDoDryBed(true);
						else
							m_pSimulation->bSetDoDryBed(false);
						break;
            	}

            	case 24: {
						// Use Murillo condition?
						strRH = strToLower(&strRH);

						if (strRH.empty())
							strErr = "line " + to_string(nLine) + ": Murillo condition";

						if  (strRH.find("y") != string::npos)
							m_pSimulation->bSetDoMurilloCondition(true);
						else
							m_pSimulation->bSetDoMurilloCondition(false);
						break;
            	}

            	case 25: {
						// Compute water density?
						strRH = strToLower(&strRH);

						if (strRH.empty())
							strErr = "line " + to_string(nLine) + ": Compute water density?";

						if  (strRH.find("y") != string::npos)
							m_pSimulation->bSetDoWaterDensity(true);
						else
							m_pSimulation->bSetDoWaterDensity(false);
						break;
            	}

            	case 26: {
						// Get the beta constant for salinity if compute water density
						if (strRH.empty())
							strErr = "line " + to_string(nLine) + ": beta constant for salinity";

						if (m_pSimulation->bGetDoWaterDensity())
						{
							m_pSimulation->dSetBetaSalinityConstant(strtod(strRH.c_str(), nullptr));
						}
						break;
            	}

            	case 27: {
						// Get the longitudinal dispersion constant, KH if compute water density
						if (strRH.empty())
							strErr = "line " + to_string(nLine) + ": longitudinal dispersion constant, KH";

						if (m_pSimulation->bGetDoWaterDensity())
						{
							m_pSimulation->dSetLongitudinalDispersionConstant(strtod(strRH.c_str(), nullptr));
						}
						break;
            	}

            	case 28: {
						// Get the sediment properties file name
						if (strRH.empty())
							strErr = "line " + to_string(nLine) + ": sediment properties file name";

						if (m_pSimulation->bGetDoWaterDensity())
						{
							m_strSedimentPropertiesFilename = strRH;
							m_strSedimentPropertiesFilename.append(".csv");
						}
						break;
            	}

            	default: {
						// More lines in the configuration file
						cerr << endl
							 << ERR << " remove extra lines" << endl;
						InStream.close();
						return true;
					}

			}
        	// Did an error occur?
        	if (! strErr.empty())
        	{
        		// Error in input to run details file
        		cerr << endl
					 << ERR << strErr << ".\nPlease edit " << SV_INI << " and change the following text:" << endl
					 << "    - '" << strRec << "'" << endl
					 << endl;
        		InStream.close();
        		return true;
        	}
    	}
	}
	// Close file
    InStream.close();

    // Finally, need to check that we have at least one raster file, so that we know the grid size and units (and preferably also the projection)
    // bool bNoRasterFiles = true;
	return false;
}


//======================================================================================================================
//!	Read Along Channel geometry file
//======================================================================================================================
bool CDataReader::bReadAlongChannelDataFile(CSimulation* m_pSimulation) {
	// Create an ifstream object
	ifstream InStream;

	// Try to open run details file for input
	InStream.open(m_strAlongChannelDataFilename.c_str(), ios::in);

	// Did it open OK?
	if (!InStream.is_open())
	{
		// Error: cannot open run details file for input
		cerr << ERR << "cannot open " << m_strAlongChannelDataFilename << " for input" << endl;
		return true;
	}

	int nCrossSectionNumber = 0;
	// int i = 0;
	// size_t nPos;
	string strRec, strErr;

	while (getline(InStream, strRec)) {

		// Trim off leading and trailing whitespace
		strRec = strTrim(&strRec);

		// If it is a blank line or a comment then ignore it
		if ((! strRec.empty()) && (strRec[0] != QUOTE1) && (strRec[0] != QUOTE2)) {

			// Create a new cross-section object and append to estuary
			m_pSimulation->AddCrossSection();

			// Update section number
			m_pSimulation->estuary[nCrossSectionNumber].nSetSectionNumber(nCrossSectionNumber);

			// Obtain the new line
			stringstream sline(strRec);

			string token;
			int j = 0;

			// Using getline for spliting the string by commas
			while (getline(sline, token, ',')) {
				double dValue = strtod(token.c_str(), nullptr);

				if (j == 0) {
					m_pSimulation->estuary[nCrossSectionNumber].dSetX(dValue);
				}

				if (j == 1) {
					m_pSimulation->estuary[nCrossSectionNumber].dSetZ(dValue);
				}

				if (j == 2) {
					m_pSimulation->estuary[nCrossSectionNumber].dSetManningNumber(dValue);
				}

				if (j == 3) {
					m_pSimulation->estuary[nCrossSectionNumber].dSetXUTM(dValue);
				}

				if (j == 4) {
					m_pSimulation->estuary[nCrossSectionNumber].dSetYUTM(dValue);
				}

				if (j == 5) {
					m_pSimulation->estuary[nCrossSectionNumber].dSetRightRBAngle(dValue);
				}

				if (j == 6) {
					m_pSimulation->estuary[nCrossSectionNumber].dSetLeftRBAngle(dValue);
				}

				if (j == 7)
				{
					if  (m_pSimulation->nGetInitialEstuarineCondition() == 0)
					{
						continue;
					}
					else if (m_pSimulation->nGetInitialEstuarineCondition() == 1) {
						m_pSimulation->m_vCrossSectionQ.push_back(dValue);
					}
					else
					{
						m_pSimulation->m_vCrossSectionElevation.push_back(dValue);
					}
				}

				// Increment counter
				j++;

			}

			// Increment counter
			nCrossSectionNumber++;
		}

	}
	m_pSimulation->m_nCrossSectionsNumber = nCrossSectionNumber - 1;
	return false;
}

//======================================================================================================================
//!	Read Cross-Section geometry file
//======================================================================================================================
bool CDataReader::bReadCrossSectionGeometryFile(CSimulation* m_pSimulation) {
	// Create an ifstream object
	ifstream InStream;

	// Try to open run details file for input
	InStream.open(m_strCrossSectionsFilename.c_str(), ios::in);

	// Did it open OK?
	if (!InStream.is_open())
	{
		// Error: cannot open run details file for input
		cerr << ERR << "cannot open " << m_strCrossSectionsFilename << " for input" << endl;
		return true;
	}

	int nLine = 0;

	int nLastElevationLine = 1;
	// int i = 0;
	// size_t nPos;
	string strRec, strErr;
	int nCrossSectionNumber = 0;
	// CEstuary* CEstuary;

	while (getline(InStream, strRec)) {

		bool bSetElevationSectionsNumber = true;

		// Trim off leading and trailing whitespace
		strRec = strTrim(&strRec);
		// If it is a blank line or a comment then ignore it
		if ((! strRec.empty()) && (strRec[0] != QUOTE1) && (strRec[0] != QUOTE2)) {
			// It isn't so increment counter
			nLine++;

			stringstream sline(strRec);

			string token;
			int j = 0;

			// Using getline for spliting the string by commas
			while (getline(sline, token, ',')) {

				double dValue = strtod(token.c_str(), nullptr);

				if (j == 0 && m_pSimulation->estuary[nCrossSectionNumber].dGetX() != dValue) {
					if (bSetElevationSectionsNumber) {
						m_pSimulation->estuary[nCrossSectionNumber].nSetElevationSectionsNumber(nLine - nLastElevationLine);
						bSetElevationSectionsNumber = false;
						nLastElevationLine = nLine;
					}
					nCrossSectionNumber++;
				}

				if (j == 2) {
					string strItem = "elevation";
					m_pSimulation->estuary[nCrossSectionNumber].dAppend2Vector(strItem, dValue);
				}

				if (j == 3) {
					string strItem = "width";
					m_pSimulation->estuary[nCrossSectionNumber].dAppend2Vector(strItem, dValue);
				}

				if (j == 4) {
					string strItem = "area";
					m_pSimulation->estuary[nCrossSectionNumber].dAppend2Vector(strItem, dValue);
				}

				if (j == 5) {
					string strItem = "perimeter";
					m_pSimulation->estuary[nCrossSectionNumber].dAppend2Vector(strItem, dValue);
				}

				if (j == 6) {
					string strItem = "hydraulic radius";
					m_pSimulation->estuary[nCrossSectionNumber].dAppend2Vector(strItem, dValue);
				}

				if (j == 7) {
					string strItem = "sigma";
					m_pSimulation->estuary[nCrossSectionNumber].dAppend2Vector(strItem, dValue);
				}

				if (j == 8) {
					string strItem = "left river bank location";
					m_pSimulation->estuary[nCrossSectionNumber].dAppend2Vector(strItem, dValue);
				}

				if (j == 9) {
					string strItem = "right river bank location";
					m_pSimulation->estuary[nCrossSectionNumber].dAppend2Vector(strItem, dValue);
				}

				if (j == 10) {
					string strItem = "beta";
					m_pSimulation->estuary[nCrossSectionNumber].dAppend2Vector(strItem, dValue);
				}

				// Increase counter
				j++;

			}
		}

	}
	// Number of elevation sections for the last Cross-Section
	m_pSimulation->estuary[nCrossSectionNumber].nSetElevationSectionsNumber(nLine - nLastElevationLine + 1);

	return false;
}


//======================================================================================================================
//!	Read Upward Boundary Condition file
//======================================================================================================================
bool CDataReader::bReadUpwardBoundaryConditionFile(CSimulation* m_pSimulation) {

	int nCrossSectionsNumber = m_pSimulation->m_nCrossSectionsNumber - 1;
	// Create an ifstream object
	ifstream InStream;

	// Try to open run details file for input
	InStream.open(m_pSimulation->m_strUpwardBoundaryConditionFilename.c_str(), ios::in);

	// Did it open OK?
	if (!InStream.is_open())
	{
		// Error: cannot open run details file for input
		cerr << ERR << "cannot open " << m_pSimulation->m_strUpwardBoundaryConditionFilename << " for input" << endl;
		return true;
	}

	int nLine = 0;
	int i = 0;
	size_t nPos;
	string strRec, strErr;

	while (getline(InStream, strRec)) {
		nLine++;

		// Trim off leading and trailing whitespace
		strRec = strTrim(&strRec);

		// If it is a blank line or a comment then ignore it
		if ((! strRec.empty()) && (strRec[0] != QUOTE1) && (strRec[0] != QUOTE2)) {
			stringstream string_line(strRec);

			string token;
			int j = 0;

			// Using getline for spliting the string line by commas
			while (getline(string_line, token, ',')) {
				double dValue = strtod(token.c_str(), nullptr);
				if (j == 0) {
					//! TODO 010: Setter and getter
					m_pSimulation->m_vUpwardBoundaryConditionTime.push_back(dValue);
				}

				if (j == 1) {
					m_pSimulation->m_vUpwardBoundaryConditionValue.push_back(dValue);
				}

				// Increment counter
				j++;
			}

			// Increment counter
			i++;
		}
	}
	return false;
}


//======================================================================================================================
//!	Read Downward Boundary Condition file
//======================================================================================================================
bool CDataReader::bReadDownwardBoundaryConditionFile(CSimulation* m_pSimulation) {

	if (m_pSimulation->nGetDownwardEstuarineCondition() != 0)
	{
		int nCrossSectionsNumber = m_pSimulation->m_nCrossSectionsNumber - 1;
		// Create an ifstream object
		ifstream InStream;

		// Try to open run details file for input
		InStream.open(m_pSimulation->m_strDownwardBoundaryConditionFilename.c_str(), ios::in);

		// Did it open OK?
		if (!InStream.is_open())
		{
			// Error: cannot open run details file for input
			cerr << ERR << "cannot open " << m_pSimulation->m_strDownwardBoundaryConditionFilename << " for input" << endl;
			return true;
		}

		int nLine = 0;
		int i = 0;
		size_t nPos;
		string strRec, strErr;

		while (getline(InStream, strRec)) {
			nLine++;

			// Trim off leading and trailing whitespace
			strRec = strTrim(&strRec);

			// If it is a blank line or a comment then ignore it
			if ((! strRec.empty()) && (strRec[0] != QUOTE1) && (strRec[0] != QUOTE2)) {
				stringstream string_line(strRec);

				string token;
				int j = 0;

				// Using getline for spliting the string line by commas
				while (getline(string_line, token, ',')) {
					double dValue = strtod(token.c_str(), nullptr);
					if (j == 0) {
						//! TODO 010: Setter and getter
						m_pSimulation->m_vDownwardBoundaryConditionTime.push_back(dValue);
					}

					if (j == 1) {
						m_pSimulation->m_vDownwardBoundaryConditionValue.push_back(dValue);
					}

					// Increment counter
					j++;
				}

				// Increment counter
				i++;
			}
		}
	}
	return false;
}


//======================================================================================================================
//! Read Hydro input file
//======================================================================================================================
bool CDataReader::bReadHydrographsFile(CSimulation* m_pSimulation) {
	// Create an ifstream object
	ifstream InStream;

	// Try to open run details file for input
	InStream.open(m_strHydroFilename.c_str(), ios::in);

	// Did it open OK?
	if (!InStream.is_open())
	{
		// Error: cannot open run details file for input
		cerr << ERR << "cannot open " << m_strHydroFilename << " for input" << endl;
		return true;
	}

	int nLine = 0;
	int i = 0;
	size_t nPos;
	string strRec, strErr;

	int nHydrographNo = -1;
	bool bReadHydrographsNo = false;
	bool bHydrographLocation = true;

	while (getline(InStream, strRec)) {
		nLine++;

		// Trim off leading and trailing whitespace
		strRec = strTrim(&strRec);

		// If it is a blank line or a comment then ignore it
		if ((! strRec.empty()) && (strRec[0] != QUOTE1) && (strRec[0] != QUOTE2)) {

			if (!bReadHydrographsNo) {
				int nValue = strtol(strRec.c_str(), nullptr, 10);
				m_pSimulation->nSetHydrographsNumber(nValue);

				for (int ii = 0; ii < nValue; ii++) {
					// Create a new hydrograph object
					m_pSimulation->AddHydrograph();
				}

				bReadHydrographsNo = true;
				continue;
			}

			stringstream string_line(strRec);

			string token;
			int j = 0;

			// Using getline for spliting the string line by commas
			while (getline(string_line, token, ',')) {
				double dValue = strtod(token.c_str(), nullptr);
				if (bHydrographLocation) {
					if (j == 0) {
						m_pSimulation->hydrographs[nHydrographNo].dSetHydrographXLocation(dValue);
					}

					if (j == 1) {
						m_pSimulation->hydrographs[nHydrographNo].dSetHydrographYLocation(dValue);
						bHydrographLocation = false;
					}
				}
				else {
					if (j == 0) {
						string strItem = "time";
						m_pSimulation->hydrographs[nHydrographNo].dAppend2Vector(strItem, dValue);
					}

					if (j == 1) {
						string strItem = "water flow";
						m_pSimulation->hydrographs[nHydrographNo].dAppend2Vector(strItem, dValue);
					}
				}

				// Increment counter
				j++;

			}
		}
		else {
			if (bReadHydrographsNo) {
				bHydrographLocation = true;
				nHydrographNo++;
			}
		}

		// Increment counter
		i++;
	}
	return false;
}


//===============================================================================================================================
//! Given a string containing time units, this returns the appropriate multiplier
//===============================================================================================================================
double CDataReader::dGetTimeMultiplier(string const *strIn)
{
	// First decide what the time units are
	int nTimeUnits = nDoTimeUnits(strIn);

	// Then return the correct multiplier, since m_dTimeStep is in hours
	switch (nTimeUnits)
	{
		default:
			return TIME_UNKNOWN;

		case TIME_SECONDS:
				return 1; // Multiplier for hours

		case TIME_HOURS:
			return 1; // Multiplier for hours

		case TIME_DAYS:
			return 24; // Multiplier for days -> hours

		case TIME_MONTHS:
			return 24 * 30.416667; // Multiplier for months -> hours (assume 30 + 5/12 day months, no leap years)

		case TIME_YEARS:
			return 24 * 365.25; // Multiplier for years -> hours
	}
}


//===============================================================================================================================
//! Given a string containing time units, this sets up the appropriate multiplier and display units for the simulation
//===============================================================================================================================
int CDataReader::nSimulationTimeMultiplier(string const *strIn)
{
    // First decide what the time units are
    int nTimeUnits = nDoTimeUnits(strIn);

    // Next set up the correct multiplier, since m_dTimeStep is in hours
    switch (nTimeUnits)
    {
    	//! TODO 023: multiplier for seconds
	    case TIME_SECONDS:
    		m_dDurationUnitsMultiplier = 1; // Multiplier for seconds
    		m_strDurationUnits = "seconds";
    		break;

        case TIME_HOURS:
            m_dDurationUnitsMultiplier = 1; // Multiplier for hours
	        m_strDurationUnits = "hours";
	        break;

        case TIME_DAYS:
            m_dDurationUnitsMultiplier = 24; // Multiplier for days -> hours
	        m_strDurationUnits = "days";
	        break;

        case TIME_MONTHS:
            m_dDurationUnitsMultiplier = 24 * 30.416667; // Multiplier for months -> hours (assume 30 + 5/12 day months, no leap years)
	        m_strDurationUnits = "months";
	        break;

        case TIME_YEARS:
            m_dDurationUnitsMultiplier = 24 * 365.25; // Multiplier for years -> hours
	        m_strDurationUnits = "years";
	        break;
    	default:
    		return RTN_ERR_TIMEUNITS;
    }

    return RTN_OK;
}

//===============================================================================================================================
//! This finds time units in a string
//===============================================================================================================================
int CDataReader::nDoTimeUnits(string const *strIn)
{
	if (strIn->find("second") != string::npos)
		return TIME_SECONDS;
	else if (strIn->find("hour") != string::npos)
        return TIME_HOURS;
	else if (strIn->find("day") != string::npos)
        return TIME_DAYS;
	else if (strIn->find("month") != string::npos)
        return TIME_MONTHS;
	else if (strIn->find("year") != string::npos)
        return TIME_YEARS;
	else
        return TIME_UNKNOWN;
}


//===============================================================================================================================
//! Trims whitespace from the left side of a string, does not change the original string
//===============================================================================================================================
string CDataReader::strTrimLeft(string const *strIn)
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
string CDataReader::strTrimRight(string const *strIn)
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
string CDataReader::strTrim(string const *strIn)
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
bool CDataReader::bParseDate(string const *strDate, int &nDay, int &nMonth, int &nYear)
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
bool CDataReader::bParseTime(string const *strTime, int &nHour, int &nMin, int &nSec)
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
