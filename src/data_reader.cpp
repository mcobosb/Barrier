/*!
 *
 * \file data_reader.cpp
 * \brief Reads configuration and input data files
 * \details Parses configuration files (.conf), cross-section geometry,
 *          along-channel data, boundary conditions, and hydrograph data.
 *          Validates input and initializes simulation parameters.
 * \author Manuel Cobos Budia

 * \date 2026
 * \copyright GNU General Public License
 *
 */

/*==============================================================================================================================

==============================================================================================================================*/
#include <cstdlib>
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

#include <limits>        // Para std::numeric_limits
#include <stdexcept>     // Para std::runtime_error
#include <cmath>         // Para std::abs

#include <algorithm>
using std::find;

#include "simulation.h"
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
	m_nSimStartYear = 0;

	m_dDurationUnitsMultiplier = 0.0;



}

//======================================================================================================================
//! The CDataReader destructor
//======================================================================================================================
CDataReader::~CDataReader() = default;


//===============================================================================================================================
//! Opens the log file
//===============================================================================================================================
void CDataReader::bOpenLogFile(CSimulation* m_pSimulation)
{
	if (m_pSimulation->m_nLogFileDetail == 0)
	{
		m_pSimulation->LogStream.open("/dev/null", ios::out | ios::trunc);

	}
	else
		m_pSimulation->LogStream.open(m_pSimulation->m_strLogFile.c_str(), ios::out | ios::trunc);
}


//======================================================================================================================
//!	Read Along Channel geometry file
//======================================================================================================================
void CDataReader::bReadAlongChannelDataFile(CSimulation* m_pSimulation) const {
	// Create an object
	ifstream InStream;

	// Try to open run details file for input
	InStream.open(m_strAlongChannelDataFilename.c_str(), ios::in);

	// Did it open OK?
	if (!InStream.is_open())
	{
		// Error: cannot open run details file for input
		cerr << ERR << "cannot open " << m_strAlongChannelDataFilename << " for input" << endl;
		return;
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
			stringstream strLine(strRec);
			string token;
			int j = 0;
			// Using get line for splitting the string line by commas
			while (getline(strLine, token, ',')) {
				double dValue = strtod(token.c_str(), nullptr);
				if (j == 0) {
					m_pSimulation->estuary[nCrossSectionNumber].dSetX(dValue);
					m_pSimulation->m_vCrossSectionX.push_back(dValue);
				}
				if (j == 1) {
					m_pSimulation->estuary[nCrossSectionNumber].dSetZ(dValue);
				}
				if (j == 2) {
					m_pSimulation->estuary[nCrossSectionNumber].dSetManningNumber(dValue);
				}
				if (j == 3) {
					m_pSimulation->estuary[nCrossSectionNumber].dSetX_UTM(dValue);
				}
				if (j == 4) {
					m_pSimulation->estuary[nCrossSectionNumber].dSetY_UTM(dValue);
				}
				if (j == 5) {
					m_pSimulation->estuary[nCrossSectionNumber].dSetRightRBAngle(dValue);
				}
				if (j == 6) {
					m_pSimulation->estuary[nCrossSectionNumber].dSetLeftRBAngle(dValue);
				}
				if (j == 7) {
					m_pSimulation->estuary[nCrossSectionNumber].dSetBeta(dValue);
				}
				// Increment counter
				j++;
			}
			// Increment counter
			nCrossSectionNumber++;
		}
	}
	m_pSimulation->m_nCrossSectionsNumber = nCrossSectionNumber;
	if  (m_pSimulation->nGetInitialEstuarineCondition() == 0)
	{
		for (int i = 0; i < nCrossSectionNumber; i++)
		{
			m_pSimulation->m_vCrossSectionQ.push_back(0.0);
			m_pSimulation->m_vCrossSectionArea.push_back(0.0);
			m_pSimulation->m_vCrossSectionWaterElevation.push_back(0.0);
			m_pSimulation->m_vCrossSectionSalinity.push_back(0.0);
			// Inicializar temperatura si corresponde
			if (m_pSimulation->m_bDoWaterTemperature) {
				m_pSimulation->m_vCrossSectionTemperature.push_back(0.0);
			}
		}
	}
}

//======================================================================================================================
//! Leer archivo de condición de frontera aguas arriba para temperatura
//======================================================================================================================
void CDataReader::bReadUpwardTemperatureBoundaryConditionFile(CSimulation* m_pSimulation) {

	ifstream InStream;
	
	InStream.open(m_pSimulation->m_strUpwardTemperatureBoundaryConditionFilename.c_str(), ios::in);
	if (!InStream.is_open()) {
		cerr << ERR << "cannot open " << m_pSimulation->m_strUpwardTemperatureBoundaryConditionFilename << " for input" << endl;
		return;
	}
	string strRec;
	while (getline(InStream, strRec)) {
		strRec = strTrim(&strRec);
		if ((!strRec.empty()) && (strRec[0] != QUOTE1) && (strRec[0] != QUOTE2)) {
			stringstream string_line(strRec);
			string token;
			int j = 0;
			while (getline(string_line, token, ',')) {
				double dValue = strtod(token.c_str(), nullptr);
				if (j == 0) m_pSimulation->m_vUpwardTemperatureBoundaryConditionTime.push_back(dValue);
				if (j == 1) m_pSimulation->m_vUpwardTemperatureBoundaryConditionValue.push_back(dValue);
				j++;
			}
		}
	}
}

//======================================================================================================================
//! Leer archivo de condición de frontera aguas abajo para temperatura
//======================================================================================================================
void CDataReader::bReadDownwardTemperatureBoundaryConditionFile(CSimulation* m_pSimulation) {
	// Create an object
	ifstream InStream;

	InStream.open(m_pSimulation->m_strDownwardTemperatureBoundaryConditionFilename.c_str(), ios::in);
	if (!InStream.is_open()) {
		cerr << ERR << "cannot open " << m_pSimulation->m_strDownwardTemperatureBoundaryConditionFilename << " for input" << endl;
		return;
	}
	string strRec;
	while (getline(InStream, strRec)) {
		strRec = strTrim(&strRec);
		if ((!strRec.empty()) && (strRec[0] != QUOTE1) && (strRec[0] != QUOTE2)) {
			stringstream string_line(strRec);
			string token;
			int j = 0;
			while (getline(string_line, token, ',')) {
				double dValue = strtod(token.c_str(), nullptr);
				if (j == 0) m_pSimulation->m_vDownwardTemperatureBoundaryConditionTime.push_back(dValue);
				if (j == 1) m_pSimulation->m_vDownwardTemperatureBoundaryConditionValue.push_back(dValue);
				j++;
			}
		}
	}
}

//======================================================================================================================
//!	Read Cross-Section geometry file
//======================================================================================================================
void CDataReader::bReadCrossSectionGeometryFile(CSimulation* m_pSimulation) const {
	// Create an object
	ifstream InStream;

	// Try to open run details file for input
	InStream.open(m_strCrossSectionsFilename.c_str(), ios::in);

	// Did it open OK?
	if (!InStream.is_open())
	{
		// Error: cannot open run details file for input
		cerr << ERR << "cannot open " << m_strCrossSectionsFilename << " for input" << endl;
		return;
	}

	int nLine = 0;

	int nLastElevationLine = 1;
	// int i = 0;
	// size_t nPos;
	string strRec, strErr;
	int nCrossSectionNumber = 0;

	while (getline(InStream, strRec)) {
		// Trim off leading and trailing whitespace
		strRec = strTrim(&strRec);
		// If it is a blank line or a comment then ignore it
		if ((! strRec.empty()) && (strRec[0] != QUOTE1) && (strRec[0] != QUOTE2)) {
			bool bSetElevationSectionsNumber = true;
			// It isn't so increment counter
			nLine++;

			stringstream strLine(strRec);

			string token;
			int j = 0;

			// Using get line for splitting the string by commas
			while (getline(strLine, token, ',')) {

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
				// Increase counter
				j++;

			}
		}

	}
	// Number of elevation sections for the last Cross-Section
	m_pSimulation->estuary[nCrossSectionNumber].nSetElevationSectionsNumber(nLine - nLastElevationLine + 1);

	// Calculate I1 pressure integral for all cross-sections after reading geometry
	for (int i = 0; i <= nCrossSectionNumber; i++) {
		m_pSimulation->estuary[i].calculateI1();
	}

}


//======================================================================================================================
//!	Read Upward Boundary Condition file
//======================================================================================================================
void CDataReader::bReadUpwardBoundaryConditionFile(CSimulation* m_pSimulation) {

	// Create an object
	ifstream InStream;

	// Try to open run details file for input
	InStream.open(m_pSimulation->m_strUpwardBoundaryConditionFilename.c_str(), ios::in);

	// Did it open OK?
	if (!InStream.is_open())
	{
		// Error: cannot open run details file for input
		cerr << ERR << "cannot open " << m_pSimulation->m_strUpwardBoundaryConditionFilename << " for input" << endl;
		return;
	}

	int nLine = 0;
	int i = 0;
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

			// Using get line for splitting the string line by commas
			while (getline(string_line, token, ',')) {
				double dValue = strtod(token.c_str(), nullptr);
				if (j == 0) {
					//! Store boundary condition time value
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
}


//======================================================================================================================
//!	Read Downward Boundary Condition file
//======================================================================================================================
void CDataReader::bReadDownwardBoundaryConditionFile(CSimulation* m_pSimulation) {

	if (m_pSimulation->nGetDownwardEstuarineCondition() != 0)
	{
		// Create an object
		ifstream InStream;

		// Try to open run details file for input
		InStream.open(m_pSimulation->m_strDownwardBoundaryConditionFilename.c_str(), ios::in);

		// Did it open OK?
		if (!InStream.is_open())
		{
			// Error: cannot open run details file for input
			cerr << ERR << "cannot open " << m_pSimulation->m_strDownwardBoundaryConditionFilename << " for input" << endl;
			return;
		}

		int nLine = 0;
		int i = 0;
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

				// Using get line for splitting the string line by commas
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
}


//======================================================================================================================
//!	Read Along Channel Sediment properties file
//======================================================================================================================
void CDataReader::bReadAlongChannelSedimentsFile(CSimulation* m_pSimulation) const {
	// Create an object
	ifstream InStream;

	// Try to open run details file for input
	InStream.open(m_strSedimentPropertiesFilename.c_str(), ios::in);

	// Did it open OK?
	if (!InStream.is_open())
	{
		// Error: cannot open run details file for input
		cerr << ERR << "cannot open " << m_strSedimentPropertiesFilename << " for input" << endl;
		return;
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

			// Obtain the new line
			stringstream strLine(strRec);

			string token;
			int j = 0;

			// Using get line for splitting the string by commas
			while (getline(strLine, token, ',')) {
				double dValue = strtod(token.c_str(), nullptr);

				if (j == 1) {
					m_pSimulation->m_vCrossSectionDaveraged.push_back(dValue);
				}

				if (j == 2) {
					m_pSimulation->m_vCrossSectionD90.push_back(dValue);
				}

				if (j == 3) {
					m_pSimulation->m_vCrossSectionD50.push_back(dValue);
				}

				if (j == 4) {
					m_pSimulation->m_vCrossSectionSedimentSigma.push_back(dValue);
				}

				if (j == 5) {
					m_pSimulation->m_vCrossSectionRhos.push_back(dValue);
				}

				if (j == 6) {
					m_pSimulation->m_vCrossSectionThickness.push_back(dValue);
				}

				// Increment counter
				j++;

			}

			// Increment counter
			nCrossSectionNumber++;
		}

	}
}

//======================================================================================================================
//! Read Hydro input file
//======================================================================================================================
void CDataReader::bReadHydrographsFile(CSimulation* m_pSimulation) const {
	if (m_pSimulation->m_bHydroFile) {
		// Create an object
		ifstream InStream;

		// Try to open run details file for input
		InStream.open(m_strHydroFilename.c_str(), ios::in);

		// Did it open OK?
		if (!InStream.is_open())
		{
			// Error: cannot open run details file for input
			cerr << ERR << "cannot open " << m_strHydroFilename << " for input" << endl;
			return;
		}

		int nLine = 0;
		int i = 0;
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

				// Using get line for splitting the string line by commas
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
		int hydrographs_no = m_pSimulation->nGetHydrographsNumber();

		//! Find the nearest cross-section of every hydrograph
		double distance_to_node = 0.0;
    	double update_distance = 1e9;  // Valor muy grande en lugar de numeric_limits::max()
    	int cs_node = -1;
		double xh, yh, xc, yc;

		for (int j = 0; j < hydrographs_no; j++) {
			xh = m_pSimulation->hydrographs[j].dGetHydrographXLocation();
			yh = m_pSimulation->hydrographs[j].dGetHydrographYLocation();
			for (int k = 0; k < m_pSimulation->m_nCrossSectionsNumber; k++) {
				xc = m_pSimulation->estuary[k].dGetX_UTM();
				yc = m_pSimulation->estuary[k].dGetY_UTM();
				distance_to_node =	(xh - xc)*(xh - xc) + (yh - yc)*(yh - yc);
				if (distance_to_node < update_distance) {
					update_distance = distance_to_node;
					cs_node = k;
				}
				else if (distance_to_node == update_distance) {
					cs_node = k;
					break; // If the distance is the same, then we can break
				}
			}
			m_pSimulation->hydrographs[j].m_nNearestCrossSectionNo = cs_node;
		}
	}
}

//======================================================================================================================
//! Leer archivo único de forzamiento de heat flux (Tair, humedad relativa, viento)
//======================================================================================================================
void CDataReader::bReadHeatFluxFile(CSimulation* m_pSimulation) {
	if (m_pSimulation->m_strHeatFluxFile.empty()) return;
	std::ifstream InStream(m_pSimulation->m_strHeatFluxFile.c_str(), std::ios::in);
	if (!InStream.is_open()) {
		std::cerr << ERR << "cannot open " << m_pSimulation->m_strHeatFluxFile << " for input" << std::endl;
		return;
	}
	std::string strRec;
	while (getline(InStream, strRec)) {
		strRec = strTrim(&strRec);
		if ((!strRec.empty()) && (strRec[0] != QUOTE1) && (strRec[0] != QUOTE2)) {
			std::stringstream string_line(strRec);
			std::string token;
			int j = 0;
			double time = 0.0, tair = 0.0, rh = 0.0, wind = 0.0;
			while (getline(string_line, token, ',')) {
				double dValue = strtod(token.c_str(), nullptr);
				if (j == 0) time = dValue;
				if (j == 1) tair = dValue;
				if (j == 2) rh = dValue;
				if (j == 3) wind = dValue;
				j++;
			}
			if (j >= 4) {
				m_pSimulation->m_vHeatFluxTime.push_back(time);
				m_pSimulation->m_vHeatFluxAirTemp.push_back(tair);
				m_pSimulation->m_vHeatFluxRelHumidity.push_back(rh);
				m_pSimulation->m_vHeatFluxWind.push_back(wind);
			}
		}
	}
}


//===============================================================================================================================
//! Given a string containing time units, this returns the appropriate multiplier
//===============================================================================================================================
double CDataReader::dGetTimeMultiplier(string const *strIn)
{
	// Then return the correct multiplier, since m_dTimeStep is in hours
	switch (nDoTimeUnits(strIn))
	{
		default:
			return TIME_UNKNOWN;

		case TIME_SECONDS:
			return 1; // Multiplier for hours

		case TIME_HOURS:
			return 3600; // Multiplier for hours

		case TIME_DAYS:
			return 3600*24; // Multiplier for days -> hours

		case TIME_MONTHS:
			return 3600*24 * 30.416667; // Multiplier for months -> hours (assume 30 + 5/12 day months, no leap years)

		case TIME_YEARS:
			return 3600*24 * 365.25; // Multiplier for years -> hours
	}
}


//===============================================================================================================================
//! Given a string containing time units, this sets up the appropriate multiplier and display units for the simulation
//===============================================================================================================================
int CDataReader::nSimulationTimeMultiplier(string const *strIn)
{
    // Next set up the correct multiplier, since m_dTimeStep is in hours
    switch (nDoTimeUnits(strIn))
    {
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
   if (size_t nStartPosition = strIn->find_first_not_of(" \t"); nStartPosition == string::npos)
      return *strIn;
   else
      return strIn->substr(nStartPosition);
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
   if (size_t nEndPos = strTmp.find_last_not_of(" \t"); nEndPos == string::npos)
      return strTmp;
   else
      return strTmp.substr(0, nEndPos + 1);
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
   vector<string> vStrTmp = VstrSplit(strDate, SLASH);
#else
   vector<string> vStrTmp = VstrSplit(strDate, SLASH);
#endif

   if (vStrTmp.size() < 3)
   {
      cerr << "date string must include day, month, and year '" << strDate << "'" << endl;
      return false;
   }

   // Sort out day
   if (! bIsStringValidInt(vStrTmp[0]))
   {
      cerr << "invalid integer for day in date '" << strDate << "'" << endl;
      return false;
   }

   nDay = stoi(vStrTmp[0]);

   if ((nDay < 1) || (nDay > 31))
   {
      cerr << "day must be between 1 and 31 in date '" << strDate << "'" << endl;
      return false;
   }

   // Sort out month
   if (! bIsStringValidInt(vStrTmp[1]))
   {
      cerr << "invalid integer for month in date '" << strDate << "'" << endl;
      return false;
   }

   nMonth = stoi(vStrTmp[1]);

   if ((nMonth < 1) || (nMonth > 12))
   {
      cerr << "month must be between 1 and 12 in date '" << strDate << "'" << endl;
      return false;
   }

   // Sort out year
   if (! bIsStringValidInt(vStrTmp[2]))
   {
      cerr << "invalid integer for year in date '" << strDate << "'" << endl;
      return false;
   }

   nYear = stoi(vStrTmp[2]);

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
   vector<string> vStrTmp = VstrSplit(strTime, DASH);

   if (vStrTmp.size() < 3)
   {
      cerr << "time string must include hours, minutes, and seconds '" << strTime << "'" << endl;
      return false;
   }

   // Sort out hour
   if (! bIsStringValidInt(vStrTmp[0]))
   {
      cerr << "invalid integer for hours in time '" << strTime << "'" << endl;
      return false;
   }

   nHour = stoi(vStrTmp[0]);

   if ((nHour < 0) || (nHour > 23))
   {
      cerr << "hour must be between 0 and 23 in time '" << strTime << "'" << endl;
      return false;
   }

   // Sort out minutes
   if (! bIsStringValidInt(vStrTmp[1]))
   {
      cerr << "invalid integer for minutes in time '" << strTime << "'" << endl;
      return false;
   }

   nMin = stoi(vStrTmp[1]);

   if ((nMin < 0) || (nMin > 59))
   {
      cerr << "minutes must be between 0 and 59 in time '" << strTime << "'" << endl;
      return false;
   }

   // Sort out seconds
   if (! bIsStringValidInt(vStrTmp[2]))
   {
      cerr << "invalid integer for seconds in time '" << strTime << "'" << endl;
      return false;
   }

   nSec = stoi(vStrTmp[2]);

   if ((nSec < 0) || (nSec > 59))
   {
      cerr << "seconds must be between 0 and 59 in time '" << strTime << "'" << endl;
      return false;
   }

   return true;
}

//======================================================================================================================
//! Read .ini file to get input and output paths
//======================================================================================================================
bool CDataReader::bReadConfigurationPaths() {
    std::ifstream configFile(".ini");
    
    if (!configFile.is_open()) {
        std::cerr << "Error: Cannot open .ini file" << std::endl;
        return false;
    }
    
    std::string line;
    bool inputPathFound = false;
    bool outputPathFound = false;
    
    while (std::getline(configFile, line)) {
        // Remove whitespace
        line = strTrim(&line);
        
        // Skip empty lines and comments
        if (line.empty() || line[0] == '#' || line[0] == ';') {
            continue;
        }
        
        // Find the '=' character
        size_t equalPos = line.find('=');
        if (equalPos == std::string::npos) {
            continue;
        }
        
        std::string key = line.substr(0, equalPos);
        std::string value = line.substr(equalPos + 1);
        
        // Trim key and value
        key = strTrim(&key);
        value = strTrim(&value);
        
        // Remove quotes if present
        if (value.front() == '"' && value.back() == '"') {
            value = value.substr(1, value.length() - 2);
        }
        
        // Check for input path
        if (key == "input_path" || key == "INPUT_PATH") {
            m_strInputPath = value;
            inputPathFound = true;
        }
        // Check for output path
        else if (key == "output_path" || key == "OUTPUT_PATH") {
            m_strOutputBasePath = value;
            outputPathFound = true;
        }
    }
    
    configFile.close();
    
    if (!inputPathFound) {
        std::cerr << "Error: input_path not found in .ini" << std::endl;
        return false;
    }
    
    if (!outputPathFound) {
        std::cerr << "Error: output_path not found in .ini" << std::endl;
        return false;
    }
    
    std::cout << "      - Input path: " << m_strInputPath << std::endl;
    std::cout << "      - Output path: " << m_strOutputBasePath << std::endl;
    
    return true;
}

// Lee el último estado de un NetCDF y lo carga en la simulación
#include <netcdf.h>
#include <cstring>

void CDataReader::bRestoreStateFromNetCDF(CSimulation* m_pSimulation, const std::string& netcdfPath) const {
    int ncid;
    int retval = nc_open(netcdfPath.c_str(), NC_NOWRITE, &ncid);
    if (retval != NC_NOERR) {
        std::cerr << "Error abriendo NetCDF para restaurar estado: " << nc_strerror(retval) << std::endl;
        return;
    }

	// Leer dimensiones
	size_t len_x = 0, len_time = 0;
	int dimid_x, dimid_time;
	nc_inq_dimid(ncid, "x", &dimid_x);
	nc_inq_dimlen(ncid, dimid_x, &len_x);
	nc_inq_dimid(ncid, "time", &dimid_time);
	nc_inq_dimlen(ncid, dimid_time, &len_time);

	// Leer el último índice temporal
	size_t last_idx = len_time > 0 ? len_time - 1 : 0;

	// Mapear variables NetCDF a los vectores de la simulación
	auto restore_var = [&](const char* varname, std::vector<double>& target) {
		int varid;
		if (nc_inq_varid(ncid, varname, &varid) == NC_NOERR) {
			std::vector<double> buffer(len_x);
			size_t start[2] = {last_idx, 0};
			size_t count[2] = {1, len_x};
			if (nc_get_vara_double(ncid, varid, start, count, buffer.data()) == NC_NOERR) {
				target = buffer;
			}
		}
	};

	restore_var("Q", m_pSimulation->m_vCrossSectionQ);
	restore_var("A", m_pSimulation->m_vCrossSectionArea);
	restore_var("eta", m_pSimulation->m_vCrossSectionWaterElevation);
	restore_var("U", m_pSimulation->m_vCrossSectionU);
	restore_var("B", m_pSimulation->m_vCrossSectionWidth);
	restore_var("Rh", m_pSimulation->m_vCrossSectionHydraulicRadius);
	restore_var("rho", m_pSimulation->m_vCrossSectionDensity);
	restore_var("S", m_pSimulation->m_vCrossSectionSalinity);
	restore_var("c", m_pSimulation->m_vCrossSectionC);
	restore_var("Qs", m_pSimulation->m_vCrossSectionQs);
	restore_var("Qb", m_pSimulation->m_vCrossSectionQb);
	restore_var("Qt", m_pSimulation->m_vCrossSectionQt);
	restore_var("level", m_pSimulation->m_vCrossSectionWaterDepth);

	// Leer el tiempo actual
	// int time_varid;
	// if (nc_inq_varid(ncid, "time", &time_varid) == NC_NOERR) {
	// 	double t;
	// 	size_t start[1] = {last_idx};
	// 	size_t count[1] = {1};
	// 	if (nc_get_vara_double(ncid, time_varid, start, count, &t) == NC_NOERR) {
	// 		m_pSimulation->m_dCurrentTime = t;
	// 	}
	// }

	nc_close(ncid);
}
