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

/**
 * @brief Construct CDataReader object
 * 
 * Initializes:
 * - Iteration counter (m_ulIter = 0)
 * - Date/time members to 0
 * - Duration units multiplier to 0.0
 */
CDataReader::CDataReader() {
	m_ulIter = 0;

	// Date properties
	m_nSimStartSec =
	m_nSimStartMin =
	m_nSimStartHour =
	m_nSimStartDay =
	m_nSimStartMonth =
	m_nSimStartYear = 0;

	m_dDurationUnitsMultiplier = 0.0;
}

/**
 * @brief Destructor (default implementation)
 */
CDataReader::~CDataReader() = default;


/**
 * @brief Open log file for simulation output
 * 
 * Behavior:
 * - If log_level = 0: Redirect to /dev/null (no logging)
 * - Otherwise: Create log file with truncation (overwrites existing)
 * 
 * @param m_pSimulation Pointer to simulation (contains log filename and detail level)
 * 
 * @note Log file contains:
 * - Initialization messages
 * - Runtime warnings/errors
 * - Progress updates (if detail > 0)
 */
void CDataReader::bOpenLogFile(CSimulation* m_pSimulation)
{
	if (m_pSimulation->m_nLogFileDetail == 0)
	{
		m_pSimulation->LogStream.open("/dev/null", ios::out | ios::trunc);

	}
	else
		m_pSimulation->LogStream.open(m_pSimulation->m_strLogFile.c_str(), ios::out | ios::trunc);
}


/**
 * @brief Read along-channel geometry CSV file
 * 
 * CSV format (no header):
 *   X, Z, Manning_n, UTM_X, UTM_Y, RightAngle, LeftAngle, Beta
 * 
 * Columns:
 * 1. X: Along-channel coordinate (m, must be monotonic increasing)
 * 2. Z: Bed elevation (m, vertical datum)
 * 3. Manning_n: Roughness coefficient (s/m^{1/3})
 * 4. UTM_X: UTM easting (m)
 * 5. UTM_Y: UTM northing (m)
 * 6. RightAngle: Right bank azimuth (degrees from north)
 * 7. LeftAngle: Left bank azimuth (degrees from north)
 * 8. Beta: Momentum correction coefficient (typically 1.0-1.1)
 * 
 * Actions:
 * - Creates CCrossSection objects (one per row)
 * - Populates m_vCrossSectionX vector
 * - Sets m_nCrossSectionsNumber
 * - Initializes state vectors (Q, A, eta, S, T) to zero if IC type = 0 (calm)
 * 
 * @param m_pSimulation Pointer to simulation object
 * 
 * @note Comments start with # or '
 * @warning File must exist and have consistent columns per row
 * @see bReadCrossSectionGeometryFile() for cross-section hydraulics
 */
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
			// Initialize temperature if applicable
			if (m_pSimulation->m_bDoWaterTemperature) {
				m_pSimulation->m_vCrossSectionTemperature.push_back(0.0);
			}
		}
	}
}

/**
 * @brief Read upstream temperature boundary condition CSV file
 * 
 * CSV format (no header): time, temperature
 * - time: Seconds since simulation start
 * - temperature: Water temperature (°C)
 * 
 * @param m_pSimulation Pointer to simulation
 * 
 * @note Used when upstream BC type = 2 (time-varying)
 */
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

/**
 * @brief Read downstream temperature boundary condition CSV file
 * 
 * CSV format (no header): time, temperature
 * - time: Seconds since simulation start
 * - temperature: Ocean/downstream water temperature (°C)
 * 
 * @param m_pSimulation Pointer to simulation
 * 
 * @note Used when downstream BC type = 2 (time-varying)
 */
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

/**
 * @brief Read cross-section hydraulic properties CSV file
 * 
 * CSV format (no header, multiple rows per cross-section):
 *   SectionID, not_used, elevation, width, area, perimeter, Rh, sigma, left_y, right_y
 * 
 * Structure:
 * - Each cross-section has N rows (one per elevation level)
 * - SectionID matches X coordinate from along-channel file
 * - Rows grouped by SectionID (must be consecutive)
 * 
 * Columns:
 * 1. SectionID: X coordinate identifier
 * 2. (unused)
 * 3. elevation: Water depth from bed (m), must be monotonic
 * 4. width: Channel width at this elevation (m)
 * 5. area: Cross-sectional area (m²)
 * 6. perimeter: Wetted perimeter (m)
 * 7. Rh: Hydraulic radius = A/P (m)
 * 8. sigma: Width function (m)
 * 9. left_y: Left bank distance from thalweg (m)
 * 10. right_y: Right bank distance from thalweg (m)
 * 
 * Post-processing:
 * - Calls calculateI1() for each cross-section (pressure integral)
 * 
 * @param m_pSimulation Pointer to simulation
 * 
 * @note Typically generated from survey data using GIS tools or HEC-RAS
 * @warning Elevations must be strictly increasing within each cross-section
 * @see CCrossSection::calculateI1() for pressure integral computation
 */
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

	// // Calculate I1 pressure integral for all cross-sections after reading geometry
	// for (int i = 0; i <= nCrossSectionNumber; i++) {
	// 	m_pSimulation->estuary[i].calculateI1();
	// }

}


/**
 * @brief Read upstream hydraulic boundary condition CSV file
 * 
 * CSV format (no header): time, value
 * - time: Seconds since simulation start
 * - value: Discharge Q (m³/s) or elevation h (m), depending on BC type
 * 
 * @param m_pSimulation Pointer to simulation
 * 
 * @note BC type set in config:
 * - Type 1: value = constant discharge
 * - Type 2: value = water elevation time series
 */
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


/**
 * @brief Read downstream hydraulic boundary condition CSV file
 * 
 * CSV format (no header): time, value
 * - time: Seconds since simulation start
 * - value: Water elevation h (m) - typically tidal time series
 * 
 * @param m_pSimulation Pointer to simulation
 * 
 * @note Only reads if downstream BC type != 0 (not free)
 * @note Common use: Astronomical tide predictions from harmonic analysis
 */
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


/**
 * @brief Read sediment properties CSV file
 * 
 * CSV format (no header):
 *   SectionID, D_avg, D90, D50, sigma, rho_s, thickness
 * 
 * Columns:
 * 1. SectionID: Cross-section identifier
 * 2. D_avg: Average sediment diameter (m)
 * 3. D90: 90th percentile grain size (m)
 * 4. D50: Median grain size (m)
 * 5. sigma: Sediment gradation coefficient (D84/D50)
 * 6. rho_s: Sediment relative density (typically 2.65 for quartz)
 * 7. thickness: Active layer thickness (m)
 * 
 * @param m_pSimulation Pointer to simulation
 * 
 * @note Required if sediment transport enabled
 * @see calculate_sediment_transport() for usage in van Rijn formula
 */
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

/**
 * @brief Read lateral inflow hydrographs CSV file
 * 
 * File structure:
 * 1. First line: N (number of tributaries)
 * 2. For each tributary:
 *    - Blank line (separator)
 *    - Location line: UTM_X, UTM_Y
 *    - Data lines: time, discharge
 * 
 * Example:
 *   3
 *   
 *   450000, 4050000
 *   0, 10.5
 *   3600, 12.3
 *   ...
 *   
 *   451000, 4051000
 *   0, 5.2
 *   ...
 * 
 * Post-processing:
 * - Finds nearest cross-section for each tributary (Euclidean distance)
 * - Stores in m_nNearestCrossSectionNo
 * 
 * @param m_pSimulation Pointer to simulation
 * 
 * @note Discharge linearly interpolated during simulation
 * @warning UTM coordinates must match projection of main channel
 */
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

/**
 * @brief Read meteorological forcing CSV file for heat flux calculations
 * 
 * CSV format (no header):
 * - 5 columns: time, T_air, RH, wind, pressure (when RH data available)
 * - 4 columns: time, T_air, wind, pressure (when calculate_rh_from_temperature enabled)
 * 
 * Columns:
 * - time: Seconds since simulation start
 * - T_air: Air temperature (°C)
 * - RH: Relative humidity (%) [optional if calculating from temperature]
 * - wind: Wind speed (m/s)
 * - pressure: Atmospheric pressure (Pa)
 * 
 * @param m_pSimulation Pointer to simulation
 * 
 * @note Flexible format:
 * - If calculate_rh_from_temperature = true: expects 4 columns (no RH)
 * - If calculate_rh_from_temperature = false: expects 5 columns (with RH)
 * - Can also handle 5 columns when calculate_rh_from_temperature = true (RH data ignored)
 * 
 * @note Used in calculateRadiativeFluxes() for:
 * - Sensible heat flux: Q_S = ρ_a·c_p·C_H·U·(T_water - T_air)
 * - Latent heat flux: Q_L = ρ_a·L_v·C_E·U·(q_sat - q_air)
 * - Longwave radiation (if not measured directly)
 * 
 * @see calculateRadiativeFluxes() in simulation.cpp
 * @see calculateDailyMinTemperatures() for RH estimation preprocessing
 */
void CDataReader::bReadHeatFluxFile(CSimulation* m_pSimulation) {
	if (m_pSimulation->m_strHeatFluxFile.empty()) return;
	std::ifstream InStream(m_pSimulation->m_strHeatFluxFile.c_str(), std::ios::in);
	if (!InStream.is_open()) {
		std::cerr << ERR << "cannot open " << m_pSimulation->m_strHeatFluxFile << " for input" << std::endl;
		return;
	}
	
	bool calculate_rh = m_pSimulation->m_bCalculateRHFromTemperature;
	int expected_cols = calculate_rh ? 4 : 5;
	
	std::string strRec;
	while (getline(InStream, strRec)) {
		strRec = strTrim(&strRec);
		if ((!strRec.empty()) && (strRec[0] != QUOTE1) && (strRec[0] != QUOTE2)) {
			std::stringstream string_line(strRec);
			std::string token;
			int j = 0;
			double time = 0.0, tair = 0.0, rh = 0.0, wind = 0.0, pressure = 0.0;
			
			while (getline(string_line, token, ',')) {
				double dValue = strtod(token.c_str(), nullptr);
				
				if (calculate_rh) {
					// Format: time, T_air, wind, pressure (4 columns)
					if (j == 0) time = dValue;
					if (j == 1) tair = dValue;
					if (j == 2) wind = dValue;
					if (j == 3) pressure = dValue;
				} else {
					// Format: time, T_air, RH, wind, pressure (5 columns)
					if (j == 0) time = dValue;
					if (j == 1) tair = dValue;
					if (j == 2) rh = dValue;
					if (j == 3) wind = dValue;
					if (j == 4) pressure = dValue;
				}
				j++;
			}
			
			// Validate column count
			if ((calculate_rh && j >= 4) || (!calculate_rh && j >= 5)) {
				m_pSimulation->m_vHeatFluxTime.push_back(time);
				m_pSimulation->m_vHeatFluxAirTemp.push_back(tair);
				if (!calculate_rh) {
					m_pSimulation->m_vHeatFluxRelHumidity.push_back(rh);
				}
				m_pSimulation->m_vHeatFluxWind.push_back(wind);
				m_pSimulation->m_vHeatFluxAtmosphericPressure.push_back(pressure);
			} else {
				std::cerr << "Warning: Skipping line with " << j << " columns (expected " 
				          << expected_cols << ")" << std::endl;
			}
		}
	}
	
	// Log information about data reading
	if (calculate_rh && m_pSimulation->m_nLogFileDetail >= 1) {
		std::cout << "      - Heat flux file read: " << m_pSimulation->m_vHeatFluxTime.size() 
		          << " records (RH will be calculated from temperature)" << std::endl;
	} else if (m_pSimulation->m_nLogFileDetail >= 1) {
		std::cout << "      - Heat flux file read: " << m_pSimulation->m_vHeatFluxTime.size() 
		          << " records (RH from data)" << std::endl;
	}
}


/**
 * @brief Get time multiplier to convert specified units to seconds
 * 
 * Supported units (case-insensitive substring match):
 * - "second" → 1
 * - "hour" → 3600
 * - "day" → 86400
 * - "month" → 2628000 (assumes 30.416667 days)
 * - "year" → 31557600 (assumes 365.25 days)
 * 
 * @param strIn String containing time unit (e.g., "hours", "days")
 * @return Multiplier to convert to seconds, or TIME_UNKNOWN if not recognized
 * 
 * @note Used for converting boundary condition time series
 */
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


/**
 * @brief Set simulation duration units and multiplier
 * 
 * Sets:
 * - m_dDurationUnitsMultiplier: Conversion factor to hours
 * - m_strDurationUnits: Display string ("seconds", "hours", "days", etc.)
 * 
 * @param strIn String containing time unit
 * @return RTN_OK on success, RTN_ERR_TIMEUNITS if unit not recognized
 * 
 * @note Called during config file parsing
 */
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

/**
 * @brief Identify time unit from string
 * 
 * @param strIn String to search (case-insensitive)
 * @return TIME_SECONDS, TIME_HOURS, TIME_DAYS, TIME_MONTHS, TIME_YEARS, or TIME_UNKNOWN
 */
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


/**
 * @brief Trim leading whitespace from string (does not modify original)
 * @param strIn Input string
 * @return Trimmed copy
 */
string CDataReader::strTrimLeft(string const *strIn)
{
   // Trim leading spaces
   if (size_t nStartPosition = strIn->find_first_not_of(" \t"); nStartPosition == string::npos)
      return *strIn;
   else
      return strIn->substr(nStartPosition);
}

/**
 * @brief Trim trailing whitespace and carriage returns (does not modify original)
 * 
 * Removes:
 * - Trailing spaces and tabs
 * - Windows carriage returns (\r)
 * 
 * @param strIn Input string
 * @return Trimmed copy
 * 
 * @note Handles files edited in Windows and transferred to Linux
 */
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

/**
 * @brief Trim leading and trailing whitespace (does not modify original)
 * 
 * Combines strTrimLeft() and strTrimRight().
 * 
 * @param strIn Input string
 * @return Trimmed copy
 * 
 * @note Most commonly used trim function in CSV parsing
 */
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


/**
 * @brief Parse date string into day, month, year components
 * 
 * Expected format: DD/MM/YYYY (or DD-MM-YYYY, platform-dependent)
 * 
 * Validation:
 * - Day: 1-31
 * - Month: 1-12
 * - Year: > 0
 * - All components must be valid integers
 * 
 * @param strDate Input date string
 * @param nDay Output: day (1-31)
 * @param nMonth Output: month (1-12)
 * @param nYear Output: year (>0)
 * @return true if valid, false if parsing/validation fails
 * 
 * @note Does NOT validate day-month combinations (e.g., allows 31/02)
 * @warning Separator is platform-dependent (SLASH macro)
 */
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

/**
 * @brief Parse time string into hour, minute, second components
 * 
 * Expected format: HH-MM-SS (dash-separated)
 * 
 * Validation:
 * - Hour: 0-23
 * - Minute: 0-59
 * - Second: 0-59
 * - All components must be valid integers
 * 
 * @param strTime Input time string
 * @param nHour Output: hour (0-23)
 * @param nMin Output: minute (0-59)
 * @param nSec Output: second (0-59)
 * @return true if valid, false if parsing/validation fails
 * 
 * @note Uses DASH separator (defined in header)
 */
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

/**
 * @brief Read .ini file for input/output directory paths
 * 
 * .ini format:
 *   input_path = "/path/to/input/data"
 *   output_path = "/path/to/output/results"
 * 
 * Features:
 * - Ignores comments (lines starting with # or ;)
 * - Strips quotes from values
 * - Case-insensitive keys (INPUT_PATH also works)
 * 
 * @return true if both paths found, false if .ini missing or incomplete
 * 
 * @note Sets:
 * - m_strInputPath: Base directory for CSV files
 * - m_strOutputBasePath: Base directory for NetCDF output
 * 
 * @see YAML config for relative path resolution
 */
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

/**
 * @brief Restore simulation state from NetCDF restart file
 * 
 * Purpose: Enable warm-start simulations (continue from previous run)
 * 
 * Procedure:
 * 1. Open NetCDF file in read-only mode
 * 2. Read dimensions (x, time)
 * 3. Extract last timestep (index = time_len - 1)
 * 4. Restore all state variables:
 *    - Hydraulics: Q, A, eta, U, B, Rh, c
 *    - Transport: S (salinity), rho (density)
 *    - Sediment: Qb, Qs, Qt
 *    - Derived: level (water depth)
 * 5. Close file
 * 
 * @param m_pSimulation Pointer to simulation object
 * @param netcdfPath Path to NetCDF restart file
 * 
 * @note Variables not found are skipped (no error)
 * @warning Assumes spatial discretization matches (same number of cross-sections)
 * @warning Time not restored (commented out) - simulation starts at t=0 by default
 * 
 * Usage:
 *   config.yaml:
 *     continue_simulation: true
 *     continue_netcdf_path: "previous_run.nc"
 * 
 * @see nDefineNetCDFFile() in data_writer.cpp for output format
 */
void CDataReader::bRestoreStateFromNetCDF(CSimulation* m_pSimulation, const std::string& netcdfPath) const {
    int ncid;
    int retval = nc_open(netcdfPath.c_str(), NC_NOWRITE, &ncid);
    if (retval != NC_NOERR) {
        std::cerr << "Error opening NetCDF for state restoration: " << nc_strerror(retval) << std::endl;
        return;
    }

	// Read dimensions
	size_t len_x = 0, len_time = 0;
	int dimid_x, dimid_time;
	nc_inq_dimid(ncid, "x", &dimid_x);
	nc_inq_dimlen(ncid, dimid_x, &len_x);
	nc_inq_dimid(ncid, "time", &dimid_time);
	nc_inq_dimlen(ncid, dimid_time, &len_time);
	
	// VALIDATION: Check spatial discretization consistency
	if (len_x != static_cast<size_t>(m_pSimulation->m_nCrossSectionsNumber)) {
		std::cerr << "\n" << std::string(80, '=') << "\n";
		std::cerr << "ERROR: Spatial discretization mismatch!\n";
		std::cerr << "NetCDF file has " << len_x << " cross-sections\n";
		std::cerr << "Current geometry has " << m_pSimulation->m_nCrossSectionsNumber << " cross-sections\n";
		std::cerr << "\nPossible causes:\n";
		std::cerr << "  1. NetCDF file from different simulation setup\n";
		std::cerr << "  2. Geometry files (along_channel_data.csv) changed since NetCDF creation\n";
		std::cerr << "  3. Wrong NetCDF file specified in continue_netcdf_path\n";
		std::cerr << "\nSolution: Use NetCDF file with " << m_pSimulation->m_nCrossSectionsNumber 
		          << " cross-sections or update geometry files\n";
		std::cerr << std::string(80, '=') << "\n";
		nc_close(ncid);
		exit(EXIT_FAILURE);
	}

	// Read the last temporal index from saved state
	size_t last_idx = len_time > 0 ? len_time - 1 : 0;
	
	// VALIDATION: Check X coordinates consistency (optional but recommended)
	int x_varid;
	if (nc_inq_varid(ncid, "x", &x_varid) == NC_NOERR) {
		std::vector<double> netcdf_x(len_x);
		if (nc_get_var_double(ncid, x_varid, netcdf_x.data()) == NC_NOERR) {
			// Compare first, middle, and last coordinates
			const double tol = 1.0;  // 1 meter tolerance
			bool x_mismatch = false;
			int mismatch_idx = -1;
			
			const std::vector<size_t> check_indices = {0, len_x/2, len_x-1};
			for (size_t i : check_indices) {
				if (i < len_x && i < m_pSimulation->m_vCrossSectionX.size()) {
					double diff = fabs(netcdf_x[i] - m_pSimulation->m_vCrossSectionX[i]);
					if (diff > tol) {
						x_mismatch = true;
						mismatch_idx = i;
						break;
					}
				}
			}
			
			if (x_mismatch) {
				std::cerr << "\n" << std::string(80, '=') << "\n";
				std::cerr << "WARNING: X coordinate mismatch detected!\n";
				std::cerr << "Cross-section " << mismatch_idx << ":\n";
				std::cerr << "  NetCDF X = " << netcdf_x[mismatch_idx] << " m\n";
				std::cerr << "  Current X = " << m_pSimulation->m_vCrossSectionX[mismatch_idx] << " m\n";
				std::cerr << "  Difference = " << fabs(netcdf_x[mismatch_idx] - m_pSimulation->m_vCrossSectionX[mismatch_idx]) << " m\n";
				std::cerr << "\nThis suggests geometry files have changed since NetCDF was created.\n";
				std::cerr << "Continuing with current geometry, but results may be inconsistent.\n";
				std::cerr << std::string(80, '=') << "\n\n";
			}
		}
	}

	// Map NetCDF variables to simulation vectors
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

	// Read current time
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
