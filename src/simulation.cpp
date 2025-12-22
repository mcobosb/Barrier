/*!
 * \file simulation.cpp
 * \brief Implementation of CSimulation class for 1D Saint-Venant estuarine modeling
 * \author Manuel Cobos Budia
 * \date 2026
 * \copyright GNU General Public License
 */

// Standard library includes
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <string>

// Standard library using declarations
using std::cerr;
using std::cout;
using std::endl;
using std::fabs;
using std::ios;
using std::pow;
using std::setprecision;
using std::setw;
using std::sqrt;
using std::to_string;

// Project includes
#include "simulation.h"
#include "cross_section.h"
#include "data_reader.h"
#include "error_handling.h"
#include "hydrograph.h"
#include "main.h"
#include "screen_presenter.h"
#include "utils.h"
#include "yaml_reader.h"

// Forward declarations for helper functions
double calc_dynamic_albedo(double lat_deg, int day_of_year, double hour_of_day);



//======================================================================================================================
//! Constructor: Initializes all simulation parameters to default values
//======================================================================================================================
CSimulation::CSimulation() {
    // Time-related variables
    m_dSimDuration = 0.0;
    m_dSimTimestep = 0.0;
    m_dTimeFactor = 0.0;
    m_dTimestep = 0.0;
    m_dLambda = 0.0;
    m_dCurrentTime = 0.0;

    // Physical process flags
    m_bDoWaterTemperature = false;
    m_bManningDependsOnLevel = false;
    m_bDoWaterSalinity = false;
    m_bDoWaterDensity = false;
    m_bDoSedimentTransport = false;
    m_bDoDryBed = false;

    // Numerical method flags
    m_bDoMcCormackLimiterFlux = false;
    m_bDoSurfaceGradientMethod = false;
    m_bDoSourceTermBalance = false;
    m_bDoBetaCoefficient = false;
    m_bDoMurilloCondition = false;

    // Model state variables
    m_nPredictor = -1;
    m_nTimeLogId = 0;
    m_nStringError = 0;
    m_nCrossSectionsNumber = 0;
    m_nHydrographsNumber = 0;
    m_nTimeId = 0;

    // Boundary and initial conditions
    m_nInitialEstuarineCondition = 0;
    m_nUpwardEstuarineCondition = 0;
    m_nDownwardEstuarineCondition = 0;
    m_nUpwardSalinityCondition = 0;
    m_nDownwardSalinityCondition = 0;

    // I/O and control flags
    m_bSaveTime = false;
    m_bHydroFile = false;
    m_bReturnError = false;
    m_nLogFileDetail = 0;

    // Numerical parameters
    m_dCourantNumber = 0.0;
    m_nEquationMcCormackLimiterFlux = 0;
    m_nPsiFormula = 0;
    m_dDeltaValue = 0.0;

    // Physical constants and parameters
    m_dBetaSalinityConstant = 0.0;
    m_dLongitudinalDispersion = 0.0;
    m_dUpwardBoundaryValue = 0.0;
    m_dNextUpwardBoundaryValue = 0.0;
    m_dDownwardBoundaryValue = 0.0;
    m_dNextDownwardBoundaryValue = 0.0;
    m_nEquationSedimentTransport = 0;

    // Tide and geometry parameters
    m_dMaxAstronomicalTide = 0.0;
    m_vEtaWidthGradientThreshold.clear();

    // Simulation start date/time (default: 2024-01-01 00:00:00)
    m_nSimStartYear = 2024;
    m_nSimStartMonth = 1;
    m_nSimStartDay = 1;
    m_nSimStartHour = 0;
    m_nSimStartMin = 0;
    m_nSimStartSec = 0;

    // Collections
    m_vOutputVariables.clear();
    estuary.clear();
    hydrographs.clear();

    // System timing
    m_tSysStartTime = 0;
    m_tSysStartLoopTime = 0;
}

//======================================================================================================================
//! Destructor: Cleans up resources and closes open files
//======================================================================================================================
CSimulation::~CSimulation() {
    // Close log file if open
    if (LogStream && LogStream.is_open()) {
        LogStream.flush();
        LogStream.close();
    }
}


//===============================================================================================================================
//! Appends a CHydrograph objet to the CSimulation
//===============================================================================================================================
void CSimulation::AddHydrograph(){
    // Create the
    CHydrograph hydrograph;
    hydrographs.push_back(hydrograph);
}


//===============================================================================================================================
//! Appends a CCrossSection objet to the estuary
//===============================================================================================================================
void CSimulation::AddCrossSection(){
    // Create the
    CCrossSection crossSection;
    estuary.push_back(crossSection);
}

//===============================================================================================================================
//! The bDoSimulation member function of CSimulation sets up and runs the simulation
//===============================================================================================================================
void CSimulation::bDoSimulation(int nArg, char const* pcArgv[]){

    //! Show starting run message
    presenter.StartingRun(nArg, pcArgv, this);

    //! Detect configuration file type and read accordingly
    std::string configFile;
    
    if (nArg > 1) {
        configFile = pcArgv[1];
    } else {
        // No config file specified, try to find config.yaml in current directory
        if (std::filesystem::exists("config.yaml")) {
            configFile = "config.yaml";
        } else {
            std::cerr << "❌ Error: No configuration file detected. config.yaml not found" << std::endl;
            std::cerr << "Usage: " << pcArgv[0] << " <config_file>" << std::endl;
            m_bReturnError = true;
            return;
        }
    }
    std::filesystem::path configPath(configFile);
    std::string extension = configPath.extension().string();
    

    std::cout << "      - Detected configuration file" << std::endl;
    CYAMLReader yamlReader;
    if (!yamlReader.loadConfiguration(configFile, this)) {
        std::cerr << "❌ Error loading YAML configuration: " 
                    << yamlReader.getErrorMessage() << std::endl;
        m_bReturnError = true;
        return;
    }

    // Transfer file paths to reader and simulation
    reader.m_strAlongChannelDataFilename = yamlReader.m_strAlongChannelDataFilename;
    reader.m_strCrossSectionsFilename = yamlReader.m_strCrossSectionGeometryFilename;
    reader.m_strSedimentPropertiesFilename = yamlReader.m_strSedimentPropertiesFilename;
    reader.m_strHydroFilename = yamlReader.m_strHydrographsFilename;

    // Transferir el threshold de gradiente desde yamlReader
    m_nThresholddBdeta = yamlReader.m_nThresholddBdeta;
    
    // Read geometry and forcing files
    CDataReader::bOpenLogFile(this);
    reader.bReadAlongChannelDataFile(this);
    reader.bReadCrossSectionGeometryFile(this);
    if (m_nUpwardEstuarineCondition > 1) {
        CDataReader::bReadUpwardBoundaryConditionFile(this);
    }
    if (m_nDownwardEstuarineCondition > 1) {
        CDataReader::bReadDownwardBoundaryConditionFile(this);
    }
    // === Temperatura: Condiciones de frontera y forzamiento ===
    if (m_bDoWaterTemperature) {
        if (m_nUpwardTemperatureCondition > 1) {
            CDataReader::bReadUpwardTemperatureBoundaryConditionFile(this);
        }
        if (m_nDownwardTemperatureCondition > 1) {
            CDataReader::bReadDownwardTemperatureBoundaryConditionFile(this);
        }
        // Forzamiento: heat_flux_file (Tair, humedad relativa, viento)
        // REQUIRED for upstream type 3 (0D reservoir model)
        if (m_nUpwardTemperatureCondition == 3 && m_strHeatFluxFile.empty()) {
            std::cerr << ERR << "Upstream temperature condition type 3 (0D reservoir model) requires heat_flux_file" << std::endl;
            std::cerr << "      Please specify heat_flux_file with meteorological forcing (time, T_air, RH, wind, pressure)" << std::endl;
            exit(EXIT_FAILURE);
        }
        if (!m_strHeatFluxFile.empty()) {
            CDataReader::bReadHeatFluxFile(this);
        }
    }
    if (m_bDoSedimentTransport) {
        reader.bReadAlongChannelSedimentsFile(this);
    }
    if (m_bHydroFile) {
        reader.bReadHydrographsFile(this);
    }

    // Initialize simulation vectors
    initializeVectors();
    precomputeEstuaryData();

    // For reservoir temperature mode (type 3), initialize upstream temperature
    if (m_nUpwardTemperatureCondition == 3) {
        updateReservoirTemperature0D();
    }

    // Restore simulation state from NetCDF if continuing previous run
    if (m_bContinueSimulation && !m_strContinueNetcdfPath.empty()) {
        std::cout << "      - Loading initial NetCDF for continuing simulation: " << m_strContinueNetcdfPath << std::endl;
        reader.bRestoreStateFromNetCDF(this, m_strContinueNetcdfPath);
    }
    
    // Apply bathymetry smoothing before calculating slopes (if enabled)
    if (bGetDoSmoothBathymetry()) {
        smoothBathymetry();
    }

    calculateBedSlope();

    // Initialize hydraulic conditions (skipped if continuing from NetCDF)
    calculateAlongEstuaryInitialConditions();

    // Calculate initial hydraulic parameters if continuing simulation
    if (m_bContinueSimulation) {
        m_nPredictor = 1;
        calculateHydraulicParameters();
    }

    // Calculate UTM coordinates for river banks (always performed)
    calculateRiverBankUTMCoordinates();

    // Calculate eta threshold for width gradient in each cross-section
    for (int i = 0; i < m_nCrossSectionsNumber; ++i) {
        estuary[i].calculateEtaMaxWidthGradient(m_dMaxAstronomicalTide, m_nThresholddBdeta, m_vEtaWidthGradientThreshold[i]);
    }

    // Generate output filename and create NetCDF file
    std::string m_strOutFileName = generateOutputFileName();
    m_strOutFile += m_strOutFileName;
    
    std::cout << "      - Output file: " << m_strOutFile << std::endl;
    writer.nDefineNetCDFFile(this);

    // Apply dry bed conditions if enabled
    if (bGetDoDryBed())
        dryArea();

    // Save initial state
    m_nTimeId = 0;
    writer.nSetOutputData(this);
    m_nTimeId++;

    int m_nStep = 1;
    cout << "    - Running" << endl;
    m_tSysStartLoopTime = time(nullptr);

    //======================================================================================================================
    //! MAIN TIME-STEPPING LOOP: McCormack predictor-corrector scheme with TVD flux limiter
    //======================================================================================================================
    while (m_dCurrentTime <= m_dSimDuration) {
        // Update reservoir temperature for upstream boundary condition type 3 (0D tank model)
        if (m_nUpwardTemperatureCondition == 3) {
            // TODO: Integrate solar radiation calculation when input data is available
            // double current_hour = get_simulation_hour();
            // int current_day = get_day_of_year();
            // double dynamic_albedo = calc_dynamic_albedo(36.5, current_day, current_hour);
            // double H_sn = solar_radiation_input * (1.0 - dynamic_albedo);
            updateReservoirTemperature0D();
        }
        
        // Update boundary conditions and calculate adaptive timestep
        calculateBoundaryConditions();
        calculateTimestep();

        //==============================================================================================================
        //! STEP 1: PREDICTOR (forward in time, forward in space)
        //==============================================================================================================
        m_nPredictor = 1;
        calculateHydraulicParameters();

        // Calculate water density from salinity and temperature at time n
        if (bGetDoWaterSalinity() || bGetDoWaterTemperature()) {
            calculate_density();
        }

        // Calculate source terms (bed slope, friction, pressure gradients)
        calculate_GS_A_terms();
        calculateFlowTerms();
        calculateSourceTerms();

        // Advance hydrodynamics to predictor state (eta*, Q*)
        calculatePredictor();
        updatePredictorBoundaries();
        if (bGetDoDryBed()) dryArea();

        // Advance transport scalars to predictor state (S*, T*)
        if (m_bDoWaterTemperature) {
            calculateRadiativeFluxes();
            calculate_temperature_predictor();
        }
        if (bGetDoWaterSalinity()) {
            calculate_salinity_predictor();
        }

        //==============================================================================================================
        //! STEP 2: CORRECTOR (backward in time, backward in space)
        //==============================================================================================================
        m_nPredictor = 2;
        calculateHydraulicParameters();

        // Calculate water density from predictor state (S*, T*)
        if (bGetDoWaterSalinity() || bGetDoWaterTemperature()) {
            calculate_density();
        }

        // Recalculate source terms with predictor state
        calculate_GS_A_terms();
        calculateFlowTerms();
        calculateSourceTerms();

        // Advance hydrodynamics to corrector state (eta**, Q**)
        calculateCorrector();
        updateCorrectorBoundaries();
        if (bGetDoDryBed()) dryArea();

        // Advance transport scalars to corrector state (S**, T**)
        if (m_bDoWaterTemperature) {
            calculate_temperature_corrector();
        }
        if (bGetDoWaterSalinity()) {
            calculate_salinity_corrector();
        }

        //==============================================================================================================
        //! STEP 3: MERGE predictor and corrector with TVD flux limiter
        //==============================================================================================================
        mergePredictorCorrector();
        mergeTracerPredictorCorrector();

        m_nPredictor = 0;
        if (bGetDoDryBed()) dryArea();
        if (bGetDoSmoothSolution()) smoothSolution();

        // Display progress and save output if needed
        AnnounceProgress();
        if (m_bSaveTime || (m_nLogFileDetail == 2)) {
            writer.nSetOutputData(this);
        }

        // Advance time
        m_dCurrentTime += m_dTimestep;
        if (m_bSaveTime) m_nTimeId++;
        m_nStep++;
    }

    // Simulation completed
    cout << endl;
    writer.nCloseNetCDFFile();

    if (LogStream.is_open()) {
        LogStream.close();
    }
}


//======================================================================================================================
//! Initialize simulation state vectors to zero or default values
//======================================================================================================================
void CSimulation::initializeVectors() {
    const int nCrossSectionsNumber = m_nCrossSectionsNumber;
    const vector<double> vZeros(static_cast<size_t>(nCrossSectionsNumber), 0.0);

    // Hydrodynamic state variables (predictor-corrector scheme)
    m_vPredictedCrossSectionArea =
    m_vCorrectedCrossSectionArea =
    m_vPredictedCrossSectionQ =
    m_vCorrectedCrossSectionQ =
    m_vPredictedCrossSectionWaterDepth = 
    m_vCrossSectionHydraulicRadius =
    m_vCrossSectionDhDx =
    m_vCrossSectionWidth =
    m_vCrossSectionU =
    m_vCrossSectionC = vZeros;

    // Transport scalars (salinity and temperature)
    m_vCrossSectionSalinity =
    m_vPredictedCrossSectionS =
    m_vCorrectedCrossSectionS =
    m_vCrossSectionTemperature =
    m_vPredictedCrossSectionT =
    m_vCorrectedCrossSectionT = vZeros;

    // Sediment transport variables
    m_vCrossSectionQb =
    m_vCrossSectionQs =
    m_vCrossSectionQt =
    m_vCrossSectionDiamX = vZeros;

    // Geometric variables
    m_vCrossSectionLeftRBLocation =
    m_vCrossSectionRightRBLocation =
    m_vCrossSectionLeftRBLocation_UTM_X =
    m_vCrossSectionLeftRBLocation_UTM_Y =
    m_vCrossSectionRightRBLocation_UTM_X =
    m_vCrossSectionRightRBLocation_UTM_Y = vZeros;

    // Bed slope and friction terms
    m_vCrossSectionBedSlope =
    m_vCrossSectionBedSlopePredictor =
    m_vCrossSectionBedSlopeCorrector =
    m_vCrossSectionFrictionSlope =
    m_vCrossSectionDX =
    m_vCrossSectionManningNumber = vZeros;

    // Source terms for Saint-Venant equations
    m_vCrossSection_gAS0 =
    m_vCrossSection_gASf =
    m_vCrossSectionF0 =
    m_vCrossSectionF1 =
    m_vCrossSectionGv0 =
    m_vCrossSectionGv1 =
    m_vCrossSectionMurilloFactor =
    m_vLateralSourcesAtT = vZeros;

    // Temporal terms for transport equations
    m_vCrossSectionSalinityASt =
    m_vCrossSectionTemperatureASt = vZeros;

    // Density and baroclinic terms
    m_vCrossSectionDRhoDx =
    m_vPredictedCrossSectionDensity = vZeros;

    m_vEtaWidthGradientThreshold.resize(nCrossSectionsNumber, 0.0);

    // Initialize density with fresh water value
    const vector<double> vRhos(static_cast<size_t>(nCrossSectionsNumber), FRESH_WATER_DENSITY);
    m_vCrossSectionDensity = vRhos;

    // TVD flux limiter factors (size n+1 for interfaces)
    const vector<double> vOnes(static_cast<size_t>(nCrossSectionsNumber + 1), 1.0);
    m_vCrossSectionD1Factor =
    m_vCrossSectionD2Factor = vOnes;

    // Output time arrays
    const int nTimestepsNumber = static_cast<int>(m_dSimDuration / m_dSimTimestep) + 1;
    m_vOutputTimesIds.resize(nTimestepsNumber);
    m_vOutputTimes.resize(nTimestepsNumber);

    for (int i = 0; i < nTimestepsNumber; i++) {
        m_vOutputTimesIds[i] = i;
        m_vOutputTimes[i] = static_cast<double>(i) * m_dSimTimestep;
    }
}

//======================================================================================================================
//! Calculate bed slope (S0) and distance increments (dx) between cross-sections
//! Uses central differences for interior nodes and one-sided differences at boundaries
//! Also calculates upwind slopes for predictor-corrector source term balance
//======================================================================================================================
void CSimulation::calculateBedSlope() {
    double dX1, dX2;
    
    // Calculate bed slope S0 = dz/dx using finite differences
    for (int i = 0; i < m_nCrossSectionsNumber; i++) {
        if (i == 0) {
            // Upstream boundary: forward difference
            dX1 = estuary[1].dGetX() - estuary[0].dGetX();
            m_vCrossSectionBedSlope[i] = (estuary[0].dGetZ() - estuary[1].dGetZ()) / dX1;
            m_vCrossSectionDX[i] = dX1;
        }
        else if (i == m_nCrossSectionsNumber - 1) {
            // Downstream boundary: backward difference
            dX2 = estuary[i].dGetX() - estuary[i-1].dGetX();
            m_vCrossSectionBedSlope[i] = (estuary[i-1].dGetZ() - estuary[i].dGetZ()) / dX2;
            m_vCrossSectionDX[i] = dX2;
        }
        else {
            // Interior nodes: central difference
            dX1 = estuary[i+1].dGetX() - estuary[i].dGetX();
            dX2 = estuary[i].dGetX() - estuary[i-1].dGetX();
            m_vCrossSectionBedSlope[i] = (estuary[i-1].dGetZ() - estuary[i+1].dGetZ()) / (dX1 + dX2);
            m_vCrossSectionDX[i] = dX1;
        }
    }
    
    // Calculate upwind slopes for predictor-corrector source term balance
    // Predictor: forward difference (i -> i+1)
    // Corrector: backward difference (i-1 -> i)
    for (int i = 0; i < m_nCrossSectionsNumber; i++) {
        if (i != m_nCrossSectionsNumber - 1) {
            double dx_forward = estuary[i+1].dGetX() - estuary[i].dGetX();
            m_vCrossSectionBedSlopePredictor[i] = (estuary[i].dGetZ() - estuary[i+1].dGetZ()) / dx_forward;
        } else {
            double dx_backward = estuary[i].dGetX() - estuary[i-1].dGetX();
            m_vCrossSectionBedSlopePredictor[i] = (estuary[i-1].dGetZ() - estuary[i].dGetZ()) / dx_backward;
        }
        
        if (i != 0) {
            double dx_backward = estuary[i].dGetX() - estuary[i-1].dGetX();
            m_vCrossSectionBedSlopeCorrector[i] = (estuary[i-1].dGetZ() - estuary[i].dGetZ()) / dx_backward;
        } else {
            double dx_forward = estuary[i+1].dGetX() - estuary[i].dGetX();
            m_vCrossSectionBedSlopeCorrector[i] = (estuary[i].dGetZ() - estuary[i+1].dGetZ()) / dx_forward;
        }
    }
}



//======================================================================================================================
//! Calculate initial hydraulic conditions along the estuary and precompute numerical constants
//! Supports three initialization modes: 0=calm water, 1=given flow, 2=given elevation
//======================================================================================================================
void CSimulation::calculateAlongEstuaryInitialConditions() {
    // Allocate water depth and elevation vectors
    m_vCrossSectionWaterDepth.resize(m_nCrossSectionsNumber);
    m_vCrossSectionWaterElevation.resize(m_nCrossSectionsNumber);

    // Allocate TVD flux limiter working vectors (used in mergePredictorCorrector)
    m_vTVD_a1_med.resize(m_nCrossSectionsNumber);
    m_vTVD_a2_med.resize(m_nCrossSectionsNumber);
    m_vTVD_alfa1_med.resize(m_nCrossSectionsNumber);
    m_vTVD_alfa2_med.resize(m_nCrossSectionsNumber);
    m_vTVD_psi1_med.resize(m_nCrossSectionsNumber);
    m_vTVD_psi2_med.resize(m_nCrossSectionsNumber);
    m_vTVD_r1_med.resize(m_nCrossSectionsNumber);
    m_vTVD_r2_med.resize(m_nCrossSectionsNumber);
    m_vTVD_fi1_med.resize(m_nCrossSectionsNumber);
    m_vTVD_fi2_med.resize(m_nCrossSectionsNumber);
    m_vTVD_Factor1.resize(m_nCrossSectionsNumber);
    m_vTVD_Factor2.resize(m_nCrossSectionsNumber);

    // Allocate salinity transport working vectors
    m_vSalinity_KAS_forward.resize(m_nCrossSectionsNumber);
    m_vSalinity_KAS_backward.resize(m_nCrossSectionsNumber);
    m_vSalinity_AUS_diff.resize(m_nCrossSectionsNumber);

    // Precompute numerical constants that remain constant during simulation
    m_vManningNumberSquared.resize(m_nCrossSectionsNumber);
    m_vInvDX.resize(m_nCrossSectionsNumber);
    m_vDxSum.resize(m_nCrossSectionsNumber);
    m_vInvDxSum.resize(m_nCrossSectionsNumber);
    m_vGtimesDX.resize(m_nCrossSectionsNumber);
    
    for (int i = 0; i < m_nCrossSectionsNumber; i++) {
        // Precompute inverse spacing (1/Δx) for gradient calculations
        if (m_vCrossSectionDX[i] > 1e-10) {
            m_vInvDX[i] = 1.0 / m_vCrossSectionDX[i];
        } else {
            m_vInvDX[i] = 0.0;
        }
        
        // Precompute g*Δx for source terms
        m_vGtimesDX[i] = G * m_vCrossSectionDX[i];
        
        // Set minimum area threshold (10% of minimum geometric area)
        DRY_AREA = 0.1 * m_vEstuaryAreas[i][0];
    }
    
    // Precompute spacing sums for central differences
    for (int i = 1; i < m_nCrossSectionsNumber - 1; i++) {
        double dx1 = m_vPositionX[i+1] - m_vPositionX[i];
        double dx2 = m_vPositionX[i] - m_vPositionX[i-1];
        m_vDxSum[i] = dx1 + dx2;
        m_vInvDxSum[i] = 1.0 / m_vDxSum[i];
    }

    // Initialize hydraulic conditions (skipped if continuing from NetCDF)
    if (!m_bContinueSimulation) {
        if (m_nInitialEstuarineCondition == 1) {
            // Mode 1: Given discharge - calculate area using Manning equation
            for (int i = 0; i < m_nCrossSectionsNumber; i++) {
                // Calculate Manning factor: Q*n/sqrt(|S0|)
                double dManningFactor = 0.0;
                if (m_vCrossSectionBedSlope[i] == 0.0) {
                    dManningFactor = m_vCrossSectionQ[i] * estuary[i].dGetManningNumber() / sqrt(1e-3);
                } else {
                    dManningFactor = m_vCrossSectionQ[i] * estuary[i].dGetManningNumber() / sqrt(fabs(m_vCrossSectionBedSlope[i]));
                }
                vector<double> vCrossSectionAreaTmp = estuary[i].vGetArea();
                vector<double> vCrossSectionHydraulicRadiusTmp = estuary[i].vGetHydraulicRadius();
                vector<double> vCrossSectionElevationTmp = estuary[i].vGetWaterDepth();
                // Compute A*R^(2/3) for each geometry node
                vector<double> dSecondTerm;
                dSecondTerm.resize(vCrossSectionAreaTmp.size());
                for (size_t j = 0; j < vCrossSectionAreaTmp.size(); j++) {
                    dSecondTerm[j] = vCrossSectionAreaTmp[j] * pow(vCrossSectionHydraulicRadiusTmp[j], 2.0/3.0);
                }
                // Interpolate to find cross-sectional area from Manning equation
                m_vCrossSectionArea[i] = linearInterpolation1d(dManningFactor, dSecondTerm, vCrossSectionAreaTmp);
                m_vCrossSectionWaterDepth[i] = linearInterpolation1d(m_vCrossSectionArea[i], vCrossSectionAreaTmp, vCrossSectionElevationTmp);
                m_vCrossSectionWidth[i] = m_vCrossSectionArea[i] / m_vCrossSectionWaterDepth[i];
                m_vCrossSectionWaterElevation[i] = m_vCrossSectionWaterDepth[i] + estuary[i].dGetZ();
            }
        }
        else if (m_nInitialEstuarineCondition == 2) {
            // Mode 2: Given elevation - calculate discharge using Manning equation
            for (int i = 0; i < m_nCrossSectionsNumber; i++) {
                vector<double> vCrossSectionAreaTmp = estuary[i].vGetArea();
                vector<double> vCrossSectionElevationTmp = estuary[i].vGetWaterDepth();
                vector<double> vCrossSectionHydraulicRadiusTmp = estuary[i].vGetHydraulicRadius();
                m_vCrossSectionArea[i] = linearInterpolation1d(m_vCrossSectionWaterElevation[i] - estuary[i].dGetZ(),vCrossSectionElevationTmp, vCrossSectionAreaTmp);
                m_vCrossSectionWidth[i] = m_vCrossSectionArea[i] / m_vCrossSectionWaterDepth[i];
                m_vCrossSectionHydraulicRadius[i] = linearInterpolation1d(m_vCrossSectionArea[i], vCrossSectionAreaTmp, vCrossSectionHydraulicRadiusTmp);
                // Manning equation: Q = (A*R^(2/3)*sqrt(|S0|)) / n * sign(S0)
                double sign_S0 = (m_vCrossSectionBedSlope[i] >= 0) ? 1.0 : -1.0;
                m_vCrossSectionQ[i] = m_vCrossSectionArea[i] * 
                                      pow(m_vCrossSectionHydraulicRadius[i], 2.0/3.0) * 
                                      sqrt(fabs(m_vCrossSectionBedSlope[i]) + 1e-10) * 
                                      sign_S0 / estuary[i].dGetManningNumber();
            }
        }
        else {
            // Mode 0: Calm water - set water level at z=0 (sea level), zero discharge
            for (int i = 0; i < m_nCrossSectionsNumber; i++) {
                vector<double> vCrossSectionAreaTmp = estuary[i].vGetArea();
                vector<double> vCrossSectionElevationTmp = estuary[i].vGetWaterDepth();
                vector<double> vCrossSectionHydraulicRadiusTmp = estuary[i].vGetHydraulicRadius();
                m_vCrossSectionWaterDepth[i] = -estuary[i].dGetZ();
                m_vCrossSectionWaterElevation[i] = 0.0;
                m_vCrossSectionArea[i] = linearInterpolation1d(m_vCrossSectionWaterDepth[i],vCrossSectionElevationTmp, vCrossSectionAreaTmp);
                m_vCrossSectionWidth[i] = m_vCrossSectionArea[i] / m_vCrossSectionWaterDepth[i];
                m_vCrossSectionHydraulicRadius[i] = linearInterpolation1d(m_vCrossSectionArea[i], vCrossSectionAreaTmp, vCrossSectionHydraulicRadiusTmp);
                m_vCrossSectionQ[i] = 0.0;
            }
        }
    }
 
    // Copy Manning coefficients from geometry data
    for (int i = 0; i < m_nCrossSectionsNumber; i++) {
        m_vCrossSectionManningNumber[i] = estuary[i].dGetManningNumber();
    }
    
    // Precompute Manning² for friction term calculations
    for (int i = 0; i < m_nCrossSectionsNumber; i++) {
        m_vManningNumberSquared[i] = pow(m_vCrossSectionManningNumber[i], 2.0);
    }

    // Calculate dimensionless sediment diameter if transport is enabled
    if (m_bDoWaterDensity) {
        for (int i = 0; i < m_nCrossSectionsNumber; i++) {
            m_vCrossSectionDiamX[i] = m_vCrossSectionD50[i] * pow((m_vCrossSectionRhos[i] - 1.0) * G / (NU * NU), 1.0/3.0);
        }
    }
}


//======================================================================================================================
//! Perform 1D linear interpolation with extrapolation at boundaries
//! Uses binary search for O(log n) lookup efficiency
//! @param dValue X-value to interpolate at
//! @param vX X-coordinates (must be sorted in ascending order)
//! @param vY Y-values corresponding to vX
//! @return Interpolated/extrapolated Y-value
//======================================================================================================================
double CSimulation::linearInterpolation1d(const double dValue, const vector<double> &vX, const vector<double> &vY) {
    auto it = std::lower_bound(vX.begin(), vX.end(), dValue);
    
    if (it == vX.begin()) {
        // Extrapolate below minimum using linear trend from first two points
        double slope = (vY[1] - vY[0]) / (vX[1] - vX[0]);
        return vY[0] + slope * (dValue - vX[0]);
    } else if (it == vX.end()) {
        // Extrapolate above maximum using linear trend from last two points
        size_t n = vX.size();
        double slope = (vY[n-1] - vY[n-2]) / (vX[n-1] - vX[n-2]);
        return vY[n-1] + slope * (dValue - vX[n-1]);
    } else {
        // Interpolate between two points
        size_t i = std::distance(vX.begin(), it) - 1;
        double slope = (vY[i+1] - vY[i]) / (vX[i+1] - vX[i]);
        return vY[i] + slope * (dValue - vX[i]);
    }
}
//======================================================================================================================
//! Calculate hydraulic parameters (width, depth, hydraulic radius) from cross-sectional area
//! Uses precomputed geometry tables for fast lookup via interpolation
//! Updates both current state and predictor/corrector states depending on m_nPredictor flag
//======================================================================================================================
/**
 * @brief Compute hydraulic parameters (width, depth, hydraulic radius, etc.) by interpolating from geometry tables
 * 
 * OPTIMIZATIONS:
 * - Binary search cache exploits spatial coherence (~30% faster)
 * - Reduced conditional branches with pointer selection
 * - Compiler hints [[likely]] for typical flow path
 * - Single interpolation pass for all variables
 * 
 * For each cross-section, interpolates hydraulic properties based on current cross-sectional area.
 * Uses cached binary search for efficient lookup in pre-computed geometry tables.
 * Handles three cases: area below minimum, area above maximum, and interpolation between table entries.
 */
//======================================================================================================================
void CSimulation::calculateHydraulicParameters() {
    const int nCrossSections = m_nCrossSectionsNumber;
    
    // Select appropriate area vector based on predictor/corrector phase
    const auto& dArea = (m_nPredictor == 1) ? m_vCrossSectionArea : m_vPredictedCrossSectionArea;
    
    // Loop through all cross-sections and interpolate hydraulic properties
    for (int i = 0; i < nCrossSections; i++) {
        const int nElevationSectionsNumber = m_vElevationSectionsCount[i];
        const double currentArea = dArea[i];
        const auto& areas = m_vEstuaryAreas[i];
        
        // Case 1: Area below minimum table value (extrapolate using first entry)
        if (areas[0] > currentArea) [[unlikely]] {
            m_vCrossSectionHydraulicRadius[i] = m_vEstuaryHydraulicRadius[i][0];
            m_vCrossSectionWidth[i] = m_vWidth[i][0];
            m_vCrossSectionWaterDepth[i] = m_vEstuaryWaterDepths[i][0];
            m_vCrossSectionLeftRBLocation[i] = m_vLeftY[i][0];
            m_vCrossSectionRightRBLocation[i] = m_vRightY[i][0];
        }
        // Case 2: Area above maximum table value (extrapolate using last entry)
        else if (areas[nElevationSectionsNumber-1] < currentArea) [[unlikely]] {
            const int lastNode = nElevationSectionsNumber - 1;
            m_vCrossSectionHydraulicRadius[i] = m_vEstuaryHydraulicRadius[i][lastNode];
            m_vCrossSectionWidth[i] = m_vWidth[i][lastNode];
            m_vCrossSectionWaterDepth[i] = m_vEstuaryWaterDepths[i][lastNode];
            m_vCrossSectionLeftRBLocation[i] = m_vLeftY[i][lastNode];
            m_vCrossSectionRightRBLocation[i] = m_vRightY[i][lastNode];
        }
        // Case 3: Area within table range (linear interpolation) - MOST COMMON PATH
        else [[likely]] {
            // OPTIMIZATION: Use cached search start position (spatial coherence)
            // Water levels change slowly between timesteps, so last index is good hint
            int j = m_vLastInterpolationIndex[i];
            
            // Validate cache: ensure j is still valid bracket
            if (j < 0 || j >= nElevationSectionsNumber-1 || 
                currentArea < areas[j] || currentArea > areas[j+1]) {
                // Cache miss: perform binary search
                auto it = std::lower_bound(areas.begin(), areas.end(), currentArea);
                j = std::distance(areas.begin(), it) - 1;
                j = std::max(0, std::min(j, nElevationSectionsNumber-2));
                m_vLastInterpolationIndex[i] = j;  // Update cache
            }
            
            // Compute interpolation factor once for all variables
            const double denom = areas[j+1] - areas[j];
            const double factor = (currentArea - areas[j]) / denom;
            const double inv_factor = 1.0 - factor;
            
            // Linear interpolation: property = (1-f)*value[j] + f*value[j+1]
            m_vCrossSectionHydraulicRadius[i] = inv_factor * m_vEstuaryHydraulicRadius[i][j] + 
                                                factor * m_vEstuaryHydraulicRadius[i][j+1];
            m_vCrossSectionWaterDepth[i] = inv_factor * m_vEstuaryWaterDepths[i][j] + 
                                           factor * m_vEstuaryWaterDepths[i][j+1];
            m_vCrossSectionWidth[i] = inv_factor * m_vWidth[i][j] + 
                                      factor * m_vWidth[i][j+1];
            m_vCrossSectionLeftRBLocation[i] = inv_factor * m_vLeftY[i][j] + 
                                               factor * m_vLeftY[i][j+1];
            m_vCrossSectionRightRBLocation[i] = inv_factor * m_vRightY[i][j] + 
                                                factor * m_vRightY[i][j+1];
        }
        
        m_vCrossSectionWaterElevation[i] = m_vCrossSectionWaterDepth[i] + m_vBedZ[i];
    }
}

//======================================================================================================================
/**
 * @brief Calculate UTM coordinates of left and right river banks
 * 
 * Computes bank positions in UTM coordinate system by projecting from thalweg (channel centerline)
 * perpendicular to the channel direction using specified angles and distances.
 * 
 * For each cross-section:
 * - Starts at thalweg UTM position (x_thalweg, y_thalweg)
 * - Projects perpendicular distance to left bank using angle_left
 * - Projects perpendicular distance to right bank using angle_right
 * 
 * Angles are in radians, measured from reference direction (typically north or channel axis).
 */
//======================================================================================================================
void CSimulation::calculateRiverBankUTMCoordinates() {
    
    for (int i = 0; i < m_nCrossSectionsNumber; i++) {
        // Get thalweg (centerline) UTM coordinates from along_channel_data.csv
        const double x_thalweg = estuary[i].dGetX_UTM();
        const double y_thalweg = estuary[i].dGetY_UTM();
        
        // Get angles in radians (AngMd for left, AngMi for right)
        const double angle_left = estuary[i].dGetLeftRBAngle();   // Already in radians
        const double angle_right = estuary[i].dGetRightRBAngle(); // Already in radians
        
        // Get distances from thalweg to river banks
        const double dist_left = fabs(m_vCrossSectionLeftRBLocation[i]);
        const double dist_right = fabs(m_vCrossSectionRightRBLocation[i]);
        
        // Calculate UTM coordinates for left bank
        // Starting from thalweg position, moving perpendicular to channel direction
        m_vCrossSectionLeftRBLocation_UTM_X[i] = x_thalweg + dist_left * cos(angle_left);
        m_vCrossSectionLeftRBLocation_UTM_Y[i] = y_thalweg + dist_left * sin(angle_left);
        
        // Calculate UTM coordinates for right bank
        m_vCrossSectionRightRBLocation_UTM_X[i] = x_thalweg + dist_right * cos(angle_right);
        m_vCrossSectionRightRBLocation_UTM_Y[i] = y_thalweg + dist_right * sin(angle_right);
    }
}

//======================================================================================================================
/**
 * @brief Linearly interpolate hydraulic parameters between two elevation table entries
 * 
 * @param dArea Current cross-sectional area to interpolate for
 * @param nCrossSection Cross-section index
 * @param nElevationNode Lower elevation table entry index (interpolates between [j] and [j+1])
 * 
 * Interpolates: hydraulic radius, width, water depth, left bank position, right bank position
 */
//======================================================================================================================
void CSimulation::interpolateHydraulicParameters(double dArea, int nCrossSection, int nElevationNode) {
    // Compute interpolation factor: f = (A - A[j]) / (A[j+1] - A[j])
    const double dInterpolationFactor = (dArea - estuary[nCrossSection].dGetArea(nElevationNode)) / 
                                        (estuary[nCrossSection].dGetArea(nElevationNode+1) - estuary[nCrossSection].dGetArea(nElevationNode));
    
    // Interpolate all hydraulic properties: property = f*(P[j+1] - P[j]) + P[j]
    m_vCrossSectionHydraulicRadius[nCrossSection] = dInterpolationFactor * (estuary[nCrossSection].dGetHydraulicRadius(nElevationNode+1) - 
                                                                            estuary[nCrossSection].dGetHydraulicRadius(nElevationNode)) + 
                                                    estuary[nCrossSection].dGetHydraulicRadius(nElevationNode);
    
    m_vCrossSectionWidth[nCrossSection] = dInterpolationFactor * (estuary[nCrossSection].dGetWidth(nElevationNode+1) - 
                                                                  estuary[nCrossSection].dGetWidth(nElevationNode)) + 
                                          estuary[nCrossSection].dGetWidth(nElevationNode);
    
    m_vCrossSectionWaterDepth[nCrossSection] = dInterpolationFactor * (estuary[nCrossSection].dGetWaterDepth(nElevationNode+1) - 
                                                                       estuary[nCrossSection].dGetWaterDepth(nElevationNode)) + 
                                               estuary[nCrossSection].dGetWaterDepth(nElevationNode);
    
    m_vCrossSectionLeftRBLocation[nCrossSection] = dInterpolationFactor * (estuary[nCrossSection].dGetLeftY(nElevationNode+1) - 
                                                                           estuary[nCrossSection].dGetLeftY(nElevationNode)) + 
                                                   estuary[nCrossSection].dGetLeftY(nElevationNode);
    
    m_vCrossSectionRightRBLocation[nCrossSection] = dInterpolationFactor * (estuary[nCrossSection].dGetRightY(nElevationNode+1) - 
                                                                            estuary[nCrossSection].dGetRightY(nElevationNode)) + 
                                                    estuary[nCrossSection].dGetRightY(nElevationNode);
}

//======================================================================================================================
/**
 * @brief Assign hydraulic parameters from first (minimum) elevation table entry
 * @param nCrossSection Cross-section index
 * 
 * Used when cross-sectional area is below the minimum table value (extrapolation).
 */
//======================================================================================================================
void CSimulation::getFirstHydraulicParameters(const int nCrossSection) {
    m_vCrossSectionHydraulicRadius[nCrossSection] = estuary[nCrossSection].dGetHydraulicRadius(0);
    m_vCrossSectionWidth[nCrossSection] = estuary[nCrossSection].dGetWidth(0);
    m_vCrossSectionWaterDepth[nCrossSection] = estuary[nCrossSection].dGetWaterDepth(0);
    m_vCrossSectionLeftRBLocation[nCrossSection] = estuary[nCrossSection].dGetLeftY(0);
    m_vCrossSectionRightRBLocation[nCrossSection] = estuary[nCrossSection].dGetRightY(0);
}

//======================================================================================================================
/**
 * @brief Assign hydraulic parameters from last (maximum) elevation table entry
 * @param nCrossSection Cross-section index
 * 
 * Used when cross-sectional area exceeds the maximum table value (extrapolation).
 */
//======================================================================================================================
void CSimulation::getLastHydraulicParameters(const int nCrossSection) {
    const int nLastNode = estuary[nCrossSection].nGetElevationSectionsNumber() - 1;
    m_vCrossSectionHydraulicRadius[nCrossSection] = estuary[nCrossSection].dGetHydraulicRadius(nLastNode);
    m_vCrossSectionWidth[nCrossSection] = estuary[nCrossSection].dGetWidth(nLastNode);
    m_vCrossSectionWaterDepth[nCrossSection] = estuary[nCrossSection].dGetWaterDepth(nLastNode);
    m_vCrossSectionLeftRBLocation[nCrossSection] = estuary[nCrossSection].dGetLeftY(nLastNode);
    m_vCrossSectionRightRBLocation[nCrossSection] = estuary[nCrossSection].dGetRightY(nLastNode);
}


//======================================================================================================================
/**
 * @brief Compute adaptive timestep satisfying CFL stability condition
 * 
 * Calculates timestep based on:
 * 1. Courant-Friedrichs-Lewy (CFL) condition: Δt ≤ C·Δx / (|u| + c)
 *    where C is Courant number, u is flow velocity, c is wave celerity
 * 2. Diffusion stability (if salinity dispersion active): Δt ≤ 0.5·Δx² / Kh
 * 3. Density-driven flow constraints (if water density active)
 * 4. Output time synchronization (reduce timestep if approaching save time)
 * 
 * Also computes flow velocity (u) and wave celerity (c) at all nodes for use in numerical scheme.
 */
//======================================================================================================================
void CSimulation::calculateTimestep() {
    double dMinTimestep = 10000000.0;
    double dWaterDensityFactor = 1.0;
  
    // ⚡ OPTIMIZATION: Cache sqrt(G) to avoid recomputing in every iteration
    static const double sqrt_G = sqrt(G);
    
    // Compute timestep based on interior nodes only (exclude boundaries i=0 and i=n-1)
    for (int i=1; i < m_nCrossSectionsNumber-1; i++) {
        if (m_vCrossSectionArea[i] > DRY_AREA) {
            const double dX = m_vCrossSectionDX[i];
            
            // Mean flow velocity: u = Q / A
            const double u = m_vCrossSectionQ[i] / m_vCrossSectionArea[i];
            m_vCrossSectionU[i] = u;

            // Shallow water wave celerity: c = sqrt(g·A/B) = sqrt(g·h) where h=A/B is mean depth
            const double c = sqrt_G * sqrt(m_vCrossSectionArea[i] / m_vCrossSectionWidth[i]);
            m_vCrossSectionC[i] = c;

            // Additional stability constraint for density-driven flows (baroclinic adjustment)
            if (m_bDoWaterDensity) {
                dWaterDensityFactor = m_dCourantNumber * dX / (m_vCrossSectionWidth[i] * m_dBetaSalinityConstant);
            }

            // ⚡ OPTIMIZATION: Precalculate characteristic speeds to avoid recomputing in max()
            // CFL timestep: Δt = C·Δx / max(|u+c|, |u-c|)
            const double lambda_plus = fabs(u + c);   // Downstream characteristic
            const double lambda_minus = fabs(u - c);  // Upstream characteristic
            const double max_lambda = std::max(lambda_plus, lambda_minus);
            
            if (const double dTimestepTmp = m_dCourantNumber * dX / max_lambda; 
                dTimestepTmp < dMinTimestep) {
                dMinTimestep = dTimestepTmp;
            }

            // Apply density factor constraint if active
            if (m_bDoWaterDensity && (dWaterDensityFactor < dMinTimestep)) {
                dMinTimestep = dWaterDensityFactor;
            }
        }
        else {
            // Handle dry/nearly-dry nodes using minimum threshold values
            const double dX = m_vCrossSectionDX[i];
            const double u_dry = DRY_Q / DRY_AREA;
            m_vCrossSectionU[i] = u_dry;

            // Wave celerity for dry zone using minimum depth: h_dry = A_dry / B
            const double h_dry = DRY_AREA / m_vCrossSectionWidth[i];
            const double c_dry = sqrt_G * sqrt(h_dry);
            m_vCrossSectionC[i] = c_dry;

            // ⚡ OPTIMIZATION: Precalculate characteristic speeds
            const double lambda_plus = fabs(u_dry + c_dry);
            const double lambda_minus = fabs(u_dry - c_dry);
            const double dTimestepTmp = m_dCourantNumber * dX / std::max(lambda_plus, lambda_minus);
            
            if (dTimestepTmp < dMinTimestep) {
                dMinTimestep = dTimestepTmp;
        }
    }}
    
    // Compute u and c at boundary nodes (for flux calculations, but not used in Δt constraint)
    for (int i : {0, m_nCrossSectionsNumber-1}) {
        if (m_vCrossSectionArea[i] != DRY_AREA) {
            m_vCrossSectionU[i] = m_vCrossSectionQ[i] / m_vCrossSectionArea[i];
            m_vCrossSectionC[i] = sqrt_G * sqrt(m_vCrossSectionArea[i] / m_vCrossSectionWidth[i]);
        } else {
            m_vCrossSectionU[i] = DRY_Q / DRY_AREA;
            const double h_dry = DRY_AREA / m_vCrossSectionWidth[i];
            m_vCrossSectionC[i] = sqrt_G * sqrt(h_dry);
        }
    }
    
    m_dTimestep = dMinTimestep;

    // Additional stability constraint for salinity diffusion (if active)
    // Criterion: Δt ≤ 0.5·Δx² / Kh to prevent numerical instabilities in diffusion term
    if (m_bDoWaterSalinity && m_dLongitudinalDispersion > 0.0) {
        for (int i = 1; i < m_nCrossSectionsNumber-1; i++) {
            double dX = m_vCrossSectionDX[i];
            double dt_diffusion = 0.5 * dX * dX / m_dLongitudinalDispersion;  // Factor 0.5 for safety margin
            if (dt_diffusion < m_dTimestep) {
                m_dTimestep = dt_diffusion;
            }
        }
    }

    // Enforce minimum timestep threshold (0.1 seconds)
    if (m_dTimestep < 0.1) {
        m_dTimestep = 0.1;
    }
    
    // Enforce maximum timestep (user-specified output timestep)
    if (m_dTimestep > m_dSimTimestep) {
        m_dTimestep = m_dSimTimestep;
    }

    // Synchronize timestep with output times if approaching scheduled save point
    if (m_bSaveAllTimesteps) {
        // Save at every computational timestep
        m_bSaveTime = true;
    }
    else if ((m_nTimeId < static_cast<int>(m_vOutputTimes.size())) && (m_dCurrentTime + m_dTimestep > m_vOutputTimes[m_nTimeId])) {
        // Reduce timestep to hit exact output time
        double dt_to_save = m_vOutputTimes[m_nTimeId] - m_dCurrentTime;
        if (dt_to_save > 1e-6) {  // Avoid dt=0 when already at save time
            m_dTimestep = dt_to_save;
        }
        m_bSaveTime = true;
    }
    else {
        m_bSaveTime = false;
    }
    
    // Compute dimensionless timestep ratio λ = Δt / min(Δx) for use in numerical scheme
    m_dLambda = m_dTimestep / dMinVectorValue(m_vCrossSectionDX);
}


//===============================================================================================================================
/**
 * @brief Calculate boundary condition values at current time by interpolating from time series
 * 
 * Performs two-step interpolation:
 * 1. Temporal interpolation: get BC value at current simulation time from input time series
 * 2. Geometric conversion: convert elevation BC to cross-sectional area using geometry tables
 * 
 * Handles both upstream (upward) and downstream (downward) boundaries.
 * For elevation BCs, also computes value at adjacent interior node for smooth transitions.
 */
//===============================================================================================================================
void CSimulation::calculateBoundaryConditions() {
    // Upstream boundary (condition type 2 = water surface elevation prescribed)
    if (nGetUpwardEstuarineCondition() == 2) {
        // Step 1: Interpolate elevation from time series at current time
        double dUpwardBoundaryValue = linearInterpolation1d(m_dCurrentTime, 
                                                           m_vUpwardBoundaryConditionTime, 
                                                           m_vUpwardBoundaryConditionValue);
        
        // Step 2: Convert water surface elevation to cross-sectional area at boundary node (i=0)
        // Uses water depth h = η - z_bed as lookup in area-depth geometry table
        m_dUpwardBoundaryValue = linearInterpolation1d(-m_vBedZ[0] + dUpwardBoundaryValue, 
                                                       m_vEstuaryWaterDepths[0],
                                                       m_vEstuaryAreas[0]);
        
        // Also compute for adjacent interior node (i=1) to ensure smooth boundary transition
        m_dNextUpwardBoundaryValue = linearInterpolation1d(-m_vBedZ[1] + dUpwardBoundaryValue, 
                                                          m_vEstuaryWaterDepths[1],
                                                          m_vEstuaryAreas[1]);
    }

    // Downstream boundary (condition type 2 = water surface elevation prescribed)
    if (nGetDownwardEstuarineCondition() == 2) {
        // Step 1: Interpolate elevation from time series at current time
        double dDownwardBoundaryValue = linearInterpolation1d(m_dCurrentTime, 
                                                             m_vDownwardBoundaryConditionTime, 
                                                             m_vDownwardBoundaryConditionValue);
        
        // Step 2: Convert water surface elevation to cross-sectional area at boundary node (i=n-1)
        m_dDownwardBoundaryValue = linearInterpolation1d(-m_vBedZ[m_nCrossSectionsNumber-1] + dDownwardBoundaryValue, 
                                                         m_vEstuaryWaterDepths[m_nCrossSectionsNumber-1], 
                                                         m_vEstuaryAreas[m_nCrossSectionsNumber-1]);
        
        // Also compute for adjacent interior node (i=n-2) to ensure smooth boundary transition
        m_dNextDownwardBoundaryValue = linearInterpolation1d(-m_vBedZ[m_nCrossSectionsNumber-2] + dDownwardBoundaryValue, 
                                                            m_vEstuaryWaterDepths[m_nCrossSectionsNumber-2], 
                                                            m_vEstuaryAreas[m_nCrossSectionsNumber-2]);
    }

    // Lateral inflows (hydrographs) at specified cross-sections
    if (m_nHydrographsNumber > 0) {
        for (int i = 0; i < m_nHydrographsNumber; i++) {
            m_vLateralSourcesAtT[hydrographs[i].m_nNearestCrossSectionNo] = linearInterpolation1d(m_dCurrentTime, 
                                                                                                  hydrographs[i].vGetTime(), 
                                                                                                  hydrographs[i].vGetQ());
        }
    }
}

//===============================================================================================================================
/**
 * @brief Enforce minimum thresholds for cross-sectional area and flow to prevent numerical issues
 * 
 * When computed area falls below physical/numerical threshold (DRY_AREA), sets minimum values:
 * - Area = DRY_AREA (small but non-zero to avoid division by zero)
 * - Flow = DRY_Q (typically zero or very small)
 * 
 * Applied to appropriate state vector depending on predictor/corrector phase:
 * - Predictor phase (m_nPredictor=1): applies to predicted variables
 * - Corrector phase (m_nPredictor=2): applies to corrected variables
 * - Normal state (m_nPredictor=0): applies to current variables
 */
//===============================================================================================================================
void CSimulation::dryArea() {
    if (m_nPredictor == 1) {
        for (int i = 0; i < m_nCrossSectionsNumber; i++) {
            if (m_vPredictedCrossSectionArea[i] <= DRY_AREA) {
                m_vPredictedCrossSectionArea[i] = DRY_AREA;
                m_vPredictedCrossSectionQ[i] = DRY_Q;
            }
        }
    }
    else if (m_nPredictor == 2) {
        for (int i = 0; i < m_nCrossSectionsNumber; i++) {
            if (m_vCorrectedCrossSectionArea[i] <= DRY_AREA) {
                m_vCorrectedCrossSectionArea[i] = DRY_AREA;
                m_vCorrectedCrossSectionQ[i] = DRY_Q;
            }
        }
    }
    else {
        for (int i = 0; i < m_nCrossSectionsNumber; i++) {
            if (m_vCrossSectionArea[i] <= DRY_AREA) {
                m_vCrossSectionArea[i] = DRY_AREA;
                m_vCrossSectionQ[i] = DRY_Q;
            }
        }
    }
}

//===============================================================================================================================
/**
 * @brief Zero out flux terms in dry/nearly-dry cross-sections to prevent spurious fluxes
 * 
 * When cross-sectional area falls below threshold, sets flux terms to zero:
 * - F0 = 0 (continuity flux)
 * - F1 = 0 (momentum flux)
 * - Gv1 = 0 (source terms)
 * 
 * This prevents unphysical fluxes in dry zones where numerical solution may be unreliable.
 * Applied to current state (predictor phase) or predicted state (corrector phase).
 */
//===============================================================================================================================
void CSimulation::dryTerms() {
    if (m_nPredictor == 1) {
        for (int i = 0; i < m_nCrossSectionsNumber-1; i++) {
            if (m_vCrossSectionArea[i] <= DRY_AREA) {
                m_vCrossSectionF0[i] = 0.0;
                m_vCrossSectionF1[i] = 0.0;
                m_vCrossSectionGv1[i] = 0.0;
            }
        }
    }
    else {
        for (int i = 0; i < m_nCrossSectionsNumber-1; i++) {
            if (m_vPredictedCrossSectionArea[i] <= DRY_AREA) {
                m_vCrossSectionF0[i] = 0.0;
                m_vCrossSectionF1[i] = 0.0;
                m_vCrossSectionGv1[i] = 0.0;
            }
        }
    }
}

//===============================================================================================================================
/**
 * @brief Apply Murillo stability condition to limit Manning coefficient in shallow/steep sections
 * 
 * Murillo et al. (2010) stability criterion prevents numerical instabilities when bed friction
 * becomes too large relative to grid resolution and flow depth.
 * 
 * Computes critical Manning value: n_crit = C_dx * sqrt(2*R^(2/3) / (g*dx))
 * If actual Manning coefficient n > n_crit, reduces n to n_crit and stores scaling factor.
 * 
 * Only applied during corrector phase (m_nPredictor != 1).
 * 
 * Reference: Murillo et al. (2010), "The influence of source terms on stability..."
 */
//===============================================================================================================================
void CSimulation::doMurilloCondition() {
    if (m_nPredictor != 1) {
        for (int i = 0; i < m_nCrossSectionsNumber; i++) {
            const double dValue = CDX_MURILLO * sqrt(
                                      2 * pow(m_vCrossSectionHydraulicRadius[i], 2.0 / 3.0) / (G * m_vCrossSectionDX[i]));

            if (dValue < estuary[i].dGetManningNumber()) {
                m_vCrossSectionMurilloFactor[i] = dValue / m_vCrossSectionManningNumber[i];
                m_vCrossSectionManningNumber[i] = dValue;
            }
            else {
                m_vCrossSectionMurilloFactor[i] = 1.0;
            }
        }
    }
}

//======================================================================================================================
/**
 * @brief Calculate bed slope and friction source terms for momentum equation
 * 
 * Computes:
 * 1. Bed slope term: gAS0 = g * A * S0
 * 2. Friction slope term: gASf = g * A * Sf, where Sf = n²|Q|Q / (A²R^(4/3))
 * 
 * Two methods:
 * - Surface term balance (bGetDoSurfaceTermBalance=true): 
 *   Uses averaged hydraulic variables between adjacent nodes and upwind bed slopes
 *   Predictor: averages with i+1 (forward), uses forward difference S0
 *   Corrector: averages with i-1 (backward), uses backward difference S0
 *   
 * - Direct computation (bGetDoSurfaceTermBalance=false):
 *   Uses local values without averaging (simpler but less accurate)
 * 
 * Optionally applies water-level-dependent Manning coefficient if m_bManningDependsOnLevel=true.
 */
//======================================================================================================================
void CSimulation::calculate_GS_A_terms() {
    
    if (bGetDoSurfaceTermBalance()) {
        // Average hydraulic variables between adjacent cross-sections for balanced source term discretization
        // Predictor: average with i+1 (forward), Corrector: average with i-1 (backward)
        for (int i = 0; i < m_nCrossSectionsNumber; i++) {
            double dMeanArea, dMeanQ, dMeanHydraulicRadius;
            
            if (m_nPredictor == 1) {
                if (i < m_nCrossSectionsNumber - 1) {
                    dMeanArea = (m_vCrossSectionArea[i] + m_vCrossSectionArea[i+1]) / 2.0;
                    dMeanQ = (m_vCrossSectionQ[i] + m_vCrossSectionQ[i+1]) / 2.0;
                    dMeanHydraulicRadius = (m_vCrossSectionHydraulicRadius[i] + m_vCrossSectionHydraulicRadius[i+1]) / 2.0;
                } else {
                    dMeanArea = m_vCrossSectionArea[i];
                    dMeanQ = m_vCrossSectionQ[i];
                    dMeanHydraulicRadius = m_vCrossSectionHydraulicRadius[i];
                }
            } else {
                if (i > 0) {
                    dMeanArea = (m_vPredictedCrossSectionArea[i] + m_vPredictedCrossSectionArea[i-1]) / 2.0;
                    dMeanQ = (m_vPredictedCrossSectionQ[i] + m_vPredictedCrossSectionQ[i-1]) / 2.0;
                    dMeanHydraulicRadius = (m_vCrossSectionHydraulicRadius[i] + m_vCrossSectionHydraulicRadius[i-1]) / 2.0;
                } else {
                    dMeanArea = m_vPredictedCrossSectionArea[i];
                    dMeanQ = m_vPredictedCrossSectionQ[i];
                    dMeanHydraulicRadius = m_vCrossSectionHydraulicRadius[i];
                }
            }
            
            // Use upwind bed slopes for consistent source term balance
            // Predictor: forward difference S0, Corrector: backward difference S0
            double S0_to_use;
            if (m_nPredictor == 1) {
                S0_to_use = m_vCrossSectionBedSlopePredictor[i];  // Forward difference (zmedp)
            } else {
                S0_to_use = m_vCrossSectionBedSlopeCorrector[i];  // Backward difference (zmedc)
            }
            
            // Bed slope source term: gAS0 = g * A * S0
            m_vCrossSection_gAS0[i] = G * dMeanArea * S0_to_use;
            
            // Water-level-dependent Manning coefficient correction (if enabled)
            double neta = 1.0;
            if (m_bManningDependsOnLevel) {
                neta = pow(n_eta(m_vPredictedCrossSectionWaterDepth[i], m_dMaxAstronomicalTide, m_vEtaWidthGradientThreshold[i]), 2.0);
            }
            
            // Friction slope: Sf = (n²*η*|Q|*Q) / (A²*R^(4/3))
            // Friction term: gASf = g*A*Sf
            if (dMeanArea > DRY_AREA && dMeanHydraulicRadius > 1e-6) {
                // ⚡ OPTIMIZATION: Precalculate repeated terms
                const double A_squared = dMeanArea * dMeanArea;
                const double Rh_power = pow(dMeanHydraulicRadius, 4.0/3.0);
                const double Q_abs_Q = dMeanQ * fabs(dMeanQ);
                
                m_vCrossSectionFrictionSlope[i] = m_vManningNumberSquared[i] * neta *
                                                 Q_abs_Q / (A_squared * Rh_power);
                m_vCrossSection_gASf[i] = G * dMeanArea * m_vCrossSectionFrictionSlope[i];
            } else {
                m_vCrossSectionFrictionSlope[i] = 0.0;
                m_vCrossSection_gASf[i] = 0.0;
            }
        }
    }
    else {
        // Direct computation without averaging (simpler but less accurate)
        for (int i = 0; i < m_nCrossSectionsNumber; i++) {
            // Bed slope term: S0 already has correct sign, no direction factor needed
            m_vCrossSection_gAS0[i] = G * m_vCrossSectionArea[i] * m_vCrossSectionBedSlope[i];
            
            if (m_nPredictor == 1) {
                double neta = 1.0;
                if (m_bManningDependsOnLevel) {
                    neta = pow(n_eta(m_vPredictedCrossSectionWaterDepth[i], m_dMaxAstronomicalTide, m_vEtaWidthGradientThreshold[i]), 2.0);
                }
                if (m_vCrossSectionArea[i] > DRY_AREA && m_vCrossSectionHydraulicRadius[i] > 1e-6) {
                    // ⚡ OPTIMIZATION: Precalculate repeated terms
                    const double A = m_vCrossSectionArea[i];
                    const double A_squared = A * A;
                    const double Rh_power = pow(m_vCrossSectionHydraulicRadius[i], 4.0/3.0);
                    const double Q_abs_Q = m_vCrossSectionQ[i] * fabs(m_vCrossSectionQ[i]);
                    
                    m_vCrossSectionFrictionSlope[i] = m_vManningNumberSquared[i] * neta *
                                                     Q_abs_Q / (A_squared * Rh_power);
                    m_vCrossSection_gASf[i] = G * A * m_vCrossSectionFrictionSlope[i];
                } else {
                    m_vCrossSectionFrictionSlope[i] = 0.0;
                    m_vCrossSection_gASf[i] = 0.0;
                }
            }
            else {
                double neta = 1.0;
                if (m_bManningDependsOnLevel) {
                    neta = pow(n_eta(m_vPredictedCrossSectionWaterDepth[i], m_dMaxAstronomicalTide, m_vEtaWidthGradientThreshold[i]), 2.0);
                }
                if (m_vPredictedCrossSectionArea[i] > DRY_AREA && m_vCrossSectionHydraulicRadius[i] > 1e-6) {
                    // ⚡ OPTIMIZATION: Precalculate repeated terms
                    const double A = m_vPredictedCrossSectionArea[i];
                    const double A_squared = A * A;
                    const double Rh_power = pow(m_vCrossSectionHydraulicRadius[i], 4.0/3.0);
                    const double Q_abs_Q = m_vPredictedCrossSectionQ[i] * fabs(m_vPredictedCrossSectionQ[i]);
                    
                    m_vCrossSectionFrictionSlope[i] = m_vManningNumberSquared[i] * neta *
                                                     Q_abs_Q / (A_squared * Rh_power);
                    m_vCrossSection_gASf[i] = G * A * m_vCrossSectionFrictionSlope[i];
                } else {
                    m_vCrossSectionFrictionSlope[i] = 0.0;
                    m_vCrossSection_gASf[i] = 0.0;
                }
            }
        }
    }
}

//======================================================================================================================
/**
 * @brief Calculate flux terms F0 and F1 for conservative form of Saint-Venant equations
 * 
 * Computes:
 * - F0 = Q (continuity equation flux: simply the discharge)
 * - F1 = β*Q²/A (momentum equation flux: momentum advection with momentum correction factor β)
 * 
 * Pressure gradient is computed more rigorously in source terms (Gv1) using -gA∂η/∂x,
 * which provides better accuracy and stability for non-prismatic channels than integral methods.
 * 
 * Applied to current state (predictor) or predicted state (corrector) depending on m_nPredictor.
 */
//======================================================================================================================
void CSimulation::calculateFlowTerms() {
    
    if (m_nPredictor == 1) {
        for (int i = 0; i < m_nCrossSectionsNumber; i++) {
            // F0 = Q (discharge, continuity flux)
            m_vCrossSectionF0[i] = m_vCrossSectionQ[i];
            
            // F1 = β*Q²/A (momentum flux)
            // Note: Pressure integral I1 term disabled for stability (see function header)
            if (m_vCrossSectionArea[i] > DRY_AREA) {
                m_vCrossSectionF1[i] = m_vBeta[i] * pow(m_vCrossSectionQ[i], 2.0) / m_vCrossSectionArea[i];
            } else {
                m_vCrossSectionF1[i] = 0.0;
            }
        }
    }
    else {
        for (int i = 0; i < m_nCrossSectionsNumber; i++) {
            // F0 = Q (discharge from predicted state)
            m_vCrossSectionF0[i] = m_vPredictedCrossSectionQ[i];
            
            // F1 = β*Q²/A (momentum flux from predicted state, I1 disabled)
            if (m_vPredictedCrossSectionArea[i] > DRY_AREA) {
                m_vCrossSectionF1[i] = m_vBeta[i] * pow(m_vPredictedCrossSectionQ[i], 2.0) / m_vPredictedCrossSectionArea[i];
            } else {
                m_vCrossSectionF1[i] = 0.0;
            }
        }
    }
}

//======================================================================================================================
/**
 * @brief Calculate source terms for momentum and continuity equations
 * 
 * Computes:
 * 1. Lateral inflow source term: Gv0 = Q_lateral / Δx
 * 2. Momentum source terms: Gv1 = -gA∂η/∂x - gASf + S_baroclinic
 *    - Barotropic pressure gradient: -gA∂η/∂x (water surface slope)
 *    - Friction term: -gASf (bottom friction)
 *    - Baroclinic pressure gradient: -gA(h/ρ)∂ρ/∂x (density-driven flow, if salinity active)
 * 
 * Uses centered differences for interior nodes, one-sided differences at boundaries.
 * Automatically selects appropriate state variables (current or predicted) based on m_nPredictor.
 */
//======================================================================================================================
void CSimulation::calculateSourceTerms() {
    // Select appropriate state vectors based on predictor/corrector phase
    const auto& area_to_use_vec = (m_nPredictor == 1) ? m_vCrossSectionArea : m_vPredictedCrossSectionArea;
    const auto& depth_to_use_vec = (m_nPredictor == 1) ? m_vCrossSectionWaterDepth : m_vPredictedCrossSectionWaterDepth;

    // 1. Water surface gradient (barotropic pressure gradient): ∂η/∂x
    //    Interior nodes: centered difference ∂η/∂x ≈ (η[i+1] - η[i-1]) / (Δx[i] + Δx[i-1])
    //    Boundary nodes: forward/backward difference
    for (int i = 1; i < m_nCrossSectionsNumber-1; i++) {
        m_vCrossSectionDhDx[i] = (m_vCrossSectionWaterElevation[i+1] - m_vCrossSectionWaterElevation[i-1]) * m_vInvDxSum[i];
    }
    m_vCrossSectionDhDx[0] = (m_vCrossSectionWaterElevation[1] - m_vCrossSectionWaterElevation[0]) * m_vInvDX[0];
    m_vCrossSectionDhDx[m_nCrossSectionsNumber-1] = (m_vCrossSectionWaterElevation[m_nCrossSectionsNumber-1] - 
                                                     m_vCrossSectionWaterElevation[m_nCrossSectionsNumber-2]) * 
                                                     m_vInvDX[m_nCrossSectionsNumber-2];

    // 2. Density gradient (baroclinic pressure gradient): ∂ρ/∂x
    //    Only computed if salinity transport is active
    if (bGetDoWaterSalinity()) {
        const auto& density_to_use = (m_nPredictor == 1) ? m_vCrossSectionDensity : m_vPredictedCrossSectionDensity;
        
        // Centered differences for interior, one-sided for boundaries
        for (int i = 1; i < m_nCrossSectionsNumber-1; i++) {
            m_vCrossSectionDRhoDx[i] = (density_to_use[i+1] - density_to_use[i-1]) * m_vInvDxSum[i];
        }
        m_vCrossSectionDRhoDx[0] = (density_to_use[1] - density_to_use[0]) * m_vInvDX[0];
        m_vCrossSectionDRhoDx[m_nCrossSectionsNumber-1] = (density_to_use[m_nCrossSectionsNumber-1] - 
                                                           density_to_use[m_nCrossSectionsNumber-2]) * 
                                                           m_vInvDX[m_nCrossSectionsNumber-2];
    }

    // 3. Assemble total source terms
    for (int i = 0; i < m_nCrossSectionsNumber; i++) {
        // Grid spacing (use appropriate Δx for each node)
        const double inv_dx = (i < m_nCrossSectionsNumber - 1) ? m_vInvDX[i] : m_vInvDX[i-1];
        
        // Continuity source: lateral inflow distributed over grid cell
        m_vCrossSectionGv0[i] = m_vLateralSourcesAtT[i] * inv_dx;
        
        // Get current state variables
        const double area_to_use = area_to_use_vec[i];
        const double depth_to_use = depth_to_use_vec[i];
        
        // ⚡ OPTIMIZATION: Cache G*A since used in both pressure and baroclinic terms
        const double G_times_area = G * area_to_use;
        
        // Momentum source: barotropic pressure gradient + friction
        // Gv1 = -gA∂η/∂x - gASf
        m_vCrossSectionGv1[i] = -G_times_area * m_vCrossSectionDhDx[i] - m_vCrossSection_gASf[i];
        
        // Add baroclinic pressure gradient if salinity is active and node is wet
        // S_baroclinic = -gA(h/ρ)∂ρ/∂x (density-driven circulation)
        if (bGetDoWaterSalinity() && depth_to_use > DRY_AREA) {
            const double current_density = (m_nPredictor == 1) ? m_vCrossSectionDensity[i] : m_vPredictedCrossSectionDensity[i];
            const double S_baroclinic = -G_times_area * (depth_to_use / current_density) * m_vCrossSectionDRhoDx[i];
            m_vCrossSectionGv1[i] += S_baroclinic;
        }
    }
}


//======================================================================================================================
/**
 * @brief Predictor step of McCormack two-step explicit scheme
 * 
 * Computes predicted values using forward spatial differences:
 * U^* = U^n - λ(F[i+1] - F[i]) + Δt·G[i]
 * 
 * where:
 * - U^n: current state vector (A, Q)
 * - U^*: predicted state
 * - λ = Δt/Δx: dimensionless timestep ratio
 * - F: flux vector
 * - G: source term vector
 * 
 * Note: Density (ρ) only affects source terms in Gv1 (baroclinic term),
 *       not the flux terms F0 and F1.
 */
//======================================================================================================================
void CSimulation::calculatePredictor() {

    const int n = m_nCrossSectionsNumber - 1;
    for (int i = 1; i < n; i++) {
        // Continuity equation: ∂A/∂t + ∂Q/∂x = G0
        // Predictor: A* = A - λ(F0[i+1] - F0[i]) + Δt·G0[i]
        m_vPredictedCrossSectionArea[i] = m_vCrossSectionArea[i] - 
            m_dLambda * (m_vCrossSectionF0[i+1] - m_vCrossSectionF0[i]) + 
            m_dTimestep * m_vCrossSectionGv0[i];  // G0 = lateral inflow (usually 0)

        // Momentum equation: ∂Q/∂t + ∂F1/∂x = G1
        // Predictor: Q* = Q - λ(F1[i+1] - F1[i]) + Δt·G1[i]
        m_vPredictedCrossSectionQ[i] = m_vCrossSectionQ[i] - 
            m_dLambda * (m_vCrossSectionF1[i+1] - m_vCrossSectionF1[i]) + 
            m_dTimestep * m_vCrossSectionGv1[i];  // G1 = pressure gradient + friction + baroclinic
    }

    // Convert predicted cross-sectional area to water depth using geometry tables
    // This is needed for subsequent hydraulic parameter calculations
    for (int i = 0; i < m_nCrossSectionsNumber; ++i) {
        m_vPredictedCrossSectionWaterDepth[i] = linearInterpolation1d(
            m_vPredictedCrossSectionArea[i],
            m_vEstuaryAreas[i],
            m_vEstuaryWaterDepths[i]
        );
    }
}


//======================================================================================================================
/**
 * @brief Corrector step of McCormack two-step explicit scheme
 * 
 * Computes corrected values using backward spatial differences:
 * U^(n+1) = U^* - λ(F*[i] - F*[i-1]) + Δt·G*[i]
 * 
 * where:
 * - U^*: predicted state from predictor step
 * - U^(n+1): corrected (final) state
 * - F*: flux computed from predicted state
 * - G*: source terms computed from predicted state
 * 
 * The final solution is typically averaged: U^(n+1) = 0.5(predictor + corrector)
 * or TVD flux limiters are applied to reduce spurious oscillations.
 */
//======================================================================================================================
void CSimulation::calculateCorrector() {

    const int n = m_nCrossSectionsNumber - 1;
    for (int i = 1; i < n; i++) {
        // Continuity equation corrector: A^(n+1) = A* - λ(F0*[i] - F0*[i-1]) + Δt·G0*[i]
        m_vCorrectedCrossSectionArea[i] = m_vPredictedCrossSectionArea[i] - 
            m_dLambda * (m_vCrossSectionF0[i] - m_vCrossSectionF0[i-1]) + 
            m_dTimestep * m_vCrossSectionGv0[i];  // G0 = lateral inflow

        // Momentum equation corrector: Q^(n+1) = Q* - λ(F1*[i] - F1*[i-1]) + Δt·G1*[i]
        // G1 includes: friction, bed slope, pressure gradients, and baroclinic term
        m_vCorrectedCrossSectionQ[i] = m_vPredictedCrossSectionQ[i] - 
            m_dLambda * (m_vCrossSectionF1[i] - m_vCrossSectionF1[i-1]) + 
            m_dTimestep * m_vCrossSectionGv1[i];
    }
}

//======================================================================================================================
/**
 * @brief Apply boundary conditions after predictor step
 * 
 * Handles upstream and downstream boundary conditions with various types:
 * - Type 0: Open/free boundary (gradient extrapolation)
 * - Type 1: Reflective/wall (Q=0) or prescribed discharge
 * - Type 2: Prescribed water surface elevation (with Manning for discharge)
 * - Type 3: Prescribed discharge (with Manning for area)
 * 
 * Also applies lateral inflows from hydrographs at interior nodes.
 * Includes selective smoothing to prevent spurious oscillations at boundaries.
 */
//======================================================================================================================
void CSimulation::updatePredictorBoundaries() {
    //==============================================================================================================
    // UPSTREAM BOUNDARY CONDITION
    //==============================================================================================================
    
    if (nGetUpwardEstuarineCondition() == 0) {
        // Type 0: Open/free boundary - linear gradient extrapolation
        // Assumes zero second derivative at boundary (smooth outflow)
        const double dQ = m_vCrossSectionQ[2] - m_vCrossSectionQ[1];
        const double dA = m_vCrossSectionArea[2] - m_vCrossSectionArea[1];
        
        // Linear extrapolation: Q[0] = Q[1] - (Q[2] - Q[1]) = 2Q[1] - Q[2]
        m_vPredictedCrossSectionQ[0] = m_vCrossSectionQ[1] - dQ;
        m_vPredictedCrossSectionArea[0] = m_vCrossSectionArea[1] - dA;
        
        // Enforce minimum area threshold for stability
        if (m_vPredictedCrossSectionArea[0] < DRY_AREA) {
            m_vPredictedCrossSectionArea[0] = m_vCrossSectionArea[1];
        }
    }
    else if (nGetUpwardEstuarineCondition() == 1) {
        // Type 1: Reflective/wall boundary - no-flow condition (solid wall)
        // Imposes zero velocity (Q=0), extrapolates water level from interior
        m_vPredictedCrossSectionQ[0] = 0.0;
        
        // Extrapolate area maintaining interior gradient
        const double dA = m_vCrossSectionArea[2] - m_vCrossSectionArea[1];
        m_vPredictedCrossSectionArea[0] = m_vCrossSectionArea[1] - dA;
        
        // Enforce minimum area threshold
        if (m_vPredictedCrossSectionArea[0] < DRY_AREA) {
            m_vPredictedCrossSectionArea[0] = m_vCrossSectionArea[1];
        }
    }
    else if (nGetUpwardEstuarineCondition() == 2) {
        // Type 2: Prescribed water surface elevation
        // Impose area (converted from elevation), compute discharge using Manning equation
        m_vPredictedCrossSectionArea[0] = m_dUpwardBoundaryValue;
        m_vCrossSectionHydraulicRadius[0] = linearInterpolation1d(
            m_vPredictedCrossSectionArea[0], 
            m_vEstuaryAreas[0], 
            m_vEstuaryHydraulicRadius[0]);
        
        // Manning equation: Q = A·R^(2/3)·√S0 / n (with correct sign for flow direction)
        double sign_S0 = (m_vCrossSectionBedSlope[0] >= 0) ? 1.0 : -1.0;
        m_vPredictedCrossSectionQ[0] = m_vPredictedCrossSectionArea[0] * 
                                      sqrt(fabs(m_vCrossSectionBedSlope[0]) + 1e-10) * 
                                      pow(m_vCrossSectionHydraulicRadius[0], 2.0/3.0) * 
                                      sign_S0 / (m_vManningN[0] + 1e-10);
    }
    else {
        // Type 3: Prescribed discharge
        // Impose discharge, compute area using inverse Manning equation
        m_vPredictedCrossSectionQ[0] = m_dUpwardBoundaryValue;
        
        // Compute Manning factor: Qn/√S0 = AR^(2/3)
        double dManningFactor = m_vPredictedCrossSectionQ[0] * m_vManningN[0] / 
                              (sqrt(fabs(m_vCrossSectionBedSlope[0]) + 1e-10));
        
        // Invert Manning equation to get area from pre-computed AR^(2/3) table
        m_vPredictedCrossSectionArea[0] = linearInterpolation1d(
            dManningFactor, 
            m_vPrecalculatedSecondTerm[0], 
            m_vEstuaryAreas[0]);
    }

    //==============================================================================================================
    // LATERAL INFLOWS (HYDROGRAPHS)
    //==============================================================================================================
    // Add lateral discharge from tributaries/sources at specified cross-sections
    for (int i = 0; i < nGetHydrographsNumber(); i++) {
        int node = hydrographs[i].m_nNearestCrossSectionNo;
        if (node >= 0 && node < m_nCrossSectionsNumber) {
            m_vPredictedCrossSectionQ[node] += m_vLateralSourcesAtT[node];
        }
    }

    //==============================================================================================================
    // DOWNSTREAM BOUNDARY CONDITION
    //==============================================================================================================
    int n = m_nCrossSectionsNumber - 1;
    
    if (nGetDownwardEstuarineCondition() == 0) {
        // Type 0: Open/free boundary - linear gradient extrapolation
        const double dQ = m_vCrossSectionQ[n-1] - m_vCrossSectionQ[n-2];
        const double dA = m_vCrossSectionArea[n-1] - m_vCrossSectionArea[n-2];
        
        m_vPredictedCrossSectionQ[n] = m_vCrossSectionQ[n-1] + dQ;
        m_vPredictedCrossSectionArea[n] = m_vCrossSectionArea[n-1] + dA;
        
        // Enforce minimum area threshold
        if (m_vPredictedCrossSectionArea[n] < DRY_AREA) {
            m_vPredictedCrossSectionArea[n] = m_vCrossSectionArea[n-1];
        }
    }
    else if (nGetDownwardEstuarineCondition() == 1) {
        // Type 1: Prescribed discharge
        m_vPredictedCrossSectionQ[n] = m_dDownwardBoundaryValue;
        
        // Compute area using inverse Manning equation
        double dManningFactor = m_vPredictedCrossSectionQ[n] * m_vManningN[n] / 
                              (sqrt(fabs(m_vCrossSectionBedSlope[n]) + 1e-10));
        
        m_vPredictedCrossSectionArea[n] = linearInterpolation1d(
            dManningFactor, 
            m_vPrecalculatedSecondTerm[n],
            m_vEstuaryAreas[n]);
    }
    else {
        // Type 2: Prescribed tidal elevation
        // Impose area (converted from elevation), compute discharge from tide rate of change
        m_vPredictedCrossSectionArea[n] = m_dDownwardBoundaryValue;
        
        // Compute time rate of change of water surface elevation: dη/dt
        double dhdt = 0.0;
        if (m_vDownwardBoundaryConditionTime.size() > 1) {
            // Find current time interval in tidal time series
            int idx = 0;
            for (size_t i = 0; i < m_vDownwardBoundaryConditionTime.size() - 1; i++) {
                if (m_dCurrentTime >= m_vDownwardBoundaryConditionTime[i] && 
                    m_dCurrentTime <= m_vDownwardBoundaryConditionTime[i+1]) {
                    idx = i;
                    break;
                }
            }
            
            // Compute temporal derivative: dη/dt ≈ Δη/Δt
            double dt = m_vDownwardBoundaryConditionTime[idx+1] - m_vDownwardBoundaryConditionTime[idx];
            if (dt > 0) {
                double dh = m_vDownwardBoundaryConditionValue[idx+1] - m_vDownwardBoundaryConditionValue[idx];
                dhdt = dh / dt;
            }
        }
        
        // Mass conservation at boundary: Q = -A·dη/dt
        // Sign convention: dη/dt > 0 (rising tide) → Q < 0 (inflow)
        //                  dη/dt < 0 (falling tide) → Q > 0 (outflow)
        double Q_tide = -m_vPredictedCrossSectionArea[n] * dhdt;
        
        // Blend with interior extrapolation (50/50 weight) for smoothness
        const double dQ = m_vCrossSectionQ[n-1] - m_vCrossSectionQ[n-2];
        double Q_extrap = m_vCrossSectionQ[n-1] + dQ;
        
        m_vPredictedCrossSectionQ[n] = 0.5 * Q_tide + 0.5 * Q_extrap;
    }
}

//===============================================================================================================================
/**
 * @brief Apply boundary conditions after corrector step
 * 
 * Similar to updatePredictorBoundaries but operates on corrected values.
 * Uses predicted state from predictor step as reference for extrapolation.
 * Includes optional selective sponge layer smoothing to dampen boundary discontinuities.
 * 
 * Boundary condition types same as predictor: 0=open, 1=reflective/discharge, 2=elevation, 3=discharge
 */
//===============================================================================================================================
void CSimulation::updateCorrectorBoundaries() {
    //==============================================================================================================
    // UPSTREAM BOUNDARY CONDITION - Corrector
    //==============================================================================================================
    
    if (nGetUpwardEstuarineCondition() == 0) {
        // Type 0: Open/free boundary - extrapolate from predicted state
        const double dQ = m_vPredictedCrossSectionQ[2] - m_vPredictedCrossSectionQ[1];
        const double dA = m_vPredictedCrossSectionArea[2] - m_vPredictedCrossSectionArea[1];
        
        m_vCorrectedCrossSectionQ[0] = m_vPredictedCrossSectionQ[1] - dQ;
        m_vCorrectedCrossSectionArea[0] = m_vPredictedCrossSectionArea[1] - dA;
        
        if (m_vCorrectedCrossSectionArea[0] < DRY_AREA) {
            m_vCorrectedCrossSectionArea[0] = m_vPredictedCrossSectionArea[1];
        }
    }
    else if (nGetUpwardEstuarineCondition() == 1) {
        // Type 1: Reflective/wall boundary (Q=0)
        m_vCorrectedCrossSectionQ[0] = 0.0;
        
        // Extrapolate area from predicted interior state
        const double dA = m_vPredictedCrossSectionArea[2] - m_vPredictedCrossSectionArea[1];
        m_vCorrectedCrossSectionArea[0] = m_vPredictedCrossSectionArea[1] - dA;
        
        if (m_vCorrectedCrossSectionArea[0] < DRY_AREA) {
            m_vCorrectedCrossSectionArea[0] = m_vPredictedCrossSectionArea[1];
        }
    }
    else if (nGetUpwardEstuarineCondition() == 2) {
        // Type 2: Prescribed elevation - impose area, compute Q with Manning
        m_vCorrectedCrossSectionArea[0] = m_dUpwardBoundaryValue;
        
        m_vCrossSectionHydraulicRadius[0] = linearInterpolation1d(
            m_vCorrectedCrossSectionArea[0], 
            m_vEstuaryAreas[0], 
            m_vEstuaryHydraulicRadius[0]);
        
        double sign_S0 = (m_vCrossSectionBedSlope[0] >= 0) ? 1.0 : -1.0;
        m_vCorrectedCrossSectionQ[0] = m_vCorrectedCrossSectionArea[0] * 
                                      sqrt(fabs(m_vCrossSectionBedSlope[0]) + 1e-10) * 
                                      pow(m_vCrossSectionHydraulicRadius[0], 2.0/3.0) * 
                                      sign_S0 / (m_vManningN[0] + 1e-10);
    }
    else {
        // Type 3: Prescribed discharge with optional selective sponge layer
        m_vCorrectedCrossSectionQ[0] = m_dUpwardBoundaryValue;
        
        double dManningFactor = m_vCorrectedCrossSectionQ[0] * m_vManningN[0] / 
                              (sqrt(fabs(m_vCrossSectionBedSlope[0]) + 1e-10));
        
        m_vCorrectedCrossSectionArea[0] = linearInterpolation1d(
            dManningFactor, 
            m_vPrecalculatedSecondTerm[0], 
            m_vEstuaryAreas[0]);
        
        // Selective sponge layer: only smooth if discontinuity > 5%
        // Damps spurious oscillations at inflow boundary
        if (m_nCrossSectionsNumber >= 3) {
            double Q_diff = fabs(m_vCorrectedCrossSectionQ[0] - m_vCrossSectionQ[1]);
            double Q_avg = 0.5 * fabs(m_vCorrectedCrossSectionQ[0]) + 0.5 * fabs(m_vCrossSectionQ[1]) + 1e-6;
            
            // Only apply smoothing if relative discontinuity exceeds 5%
            if (Q_diff / Q_avg > 0.05) {
                double alpha = 0.8;  // Smoothing weight (80% computed, 20% smoothed)
                double Q_computed = m_vCorrectedCrossSectionQ[1];
                double Q_smooth = 0.5 * (m_vCorrectedCrossSectionQ[0] + m_vCrossSectionQ[2]);
                m_vCorrectedCrossSectionQ[1] = alpha * Q_computed + (1 - alpha) * Q_smooth;
            }
        }
    }

    //==============================================================================================================
    // LATERAL INFLOWS (HYDROGRAPHS)
    //==============================================================================================================
    for (int i = 0; i < nGetHydrographsNumber(); i++) {
        int node = hydrographs[i].m_nNearestCrossSectionNo;
        if (node >= 0 && node < m_nCrossSectionsNumber) {
            m_vCorrectedCrossSectionQ[node] += m_vLateralSourcesAtT[node];
        }
    }

    //==============================================================================================================
    // DOWNSTREAM BOUNDARY CONDITION - Corrector
    //==============================================================================================================
    int n = m_nCrossSectionsNumber - 1;
    
    if (nGetDownwardEstuarineCondition() == 0) {
        // Type 0: Open/free boundary
        const double dQ = m_vPredictedCrossSectionQ[n-1] - m_vPredictedCrossSectionQ[n-2];
        const double dA = m_vPredictedCrossSectionArea[n-1] - m_vPredictedCrossSectionArea[n-2];
        
        m_vCorrectedCrossSectionQ[n] = m_vPredictedCrossSectionQ[n-1] + dQ;
        m_vCorrectedCrossSectionArea[n] = m_vPredictedCrossSectionArea[n-1] + dA;
        
        if (m_vCorrectedCrossSectionArea[n] < DRY_AREA) {
            m_vCorrectedCrossSectionArea[n] = m_vPredictedCrossSectionArea[n-1];
        }
    }
    else if (nGetDownwardEstuarineCondition() == 1) {
        // Type 1: Prescribed discharge
        m_vCorrectedCrossSectionQ[n] = m_dDownwardBoundaryValue;
        
        const double dManningFactor = m_vCorrectedCrossSectionQ[n] * m_vManningN[n] / 
                                     (sqrt(fabs(m_vCrossSectionBedSlope[n]) + 1e-10));
        
        m_vCorrectedCrossSectionArea[n] = linearInterpolation1d(
            dManningFactor, 
            m_vPrecalculatedSecondTerm[n],
            m_vEstuaryAreas[n]);
    }
    else {
        // Type 2: Prescribed tidal elevation with tide rate of change
        m_vCorrectedCrossSectionArea[n] = m_dDownwardBoundaryValue;
        
        // Compute time rate of change of tidal elevation
        double dhdt = 0.0;
        if (m_vDownwardBoundaryConditionTime.size() > 1) {
            // Find current time interval in boundary condition series
            int idx = 0;
            for (size_t i = 0; i < m_vDownwardBoundaryConditionTime.size() - 1; i++) {
                if (m_dCurrentTime >= m_vDownwardBoundaryConditionTime[i] && 
                    m_dCurrentTime <= m_vDownwardBoundaryConditionTime[i+1]) {
                    idx = i;
                    break;
                }
            }
            
            // Compute dη/dt
            double dt = m_vDownwardBoundaryConditionTime[idx+1] - m_vDownwardBoundaryConditionTime[idx];
            if (dt > 0) {
                double dh = m_vDownwardBoundaryConditionValue[idx+1] - m_vDownwardBoundaryConditionValue[idx];
                dhdt = dh / dt;
            }
        }
        
        // Discharge from mass conservation: Q = -A·dη/dt
        double Q_tide = -m_vCorrectedCrossSectionArea[n] * dhdt;
        
        // Blend with extrapolation from predicted interior
        const double dQ = m_vPredictedCrossSectionQ[n-1] - m_vPredictedCrossSectionQ[n-2];
        double Q_extrap = m_vPredictedCrossSectionQ[n-1] + dQ;
        
        m_vCorrectedCrossSectionQ[n] = 0.5 * Q_tide + 0.5 * Q_extrap;
    }
}



//======================================================================================================================
/**
 * @brief Predictor step for salinity transport equation
 * 
 * Solves advection-diffusion equation for salinity:
 * ∂(AS)/∂t + ∂(QS)/∂x = ∂/∂x(Kh·A·∂S/∂x)
 * 
 * Uses upwind advection scheme based on flow direction:
 * - If u > 0: backward difference (i - i-1)
 * - If u < 0: forward difference (i+1 - i)
 * 
 * Diffusion uses centered differences for stability.
 * Conserves salt mass: updates A·S, then divides by A to get S.
 * Applies prescribed boundary conditions at upstream/downstream.
 */
//======================================================================================================================
void CSimulation::calculate_salinity_predictor() {
    // Compute advection and diffusion terms (centered scheme for diffusion, upwind for advection)
    for (int i = 1; i < m_nCrossSectionsNumber-1; ++i) {
        double dx = m_vCrossSectionDX[i];
        double u = m_vCrossSectionU[i];
        
        // Upwind advection: choose direction based on flow
        double adv = -u * (m_vCrossSectionSalinity[i] - m_vCrossSectionSalinity[i-1]) / dx;
        if (u < 0) adv = -u * (m_vCrossSectionSalinity[i+1] - m_vCrossSectionSalinity[i]) / dx;
        
        // Centered diffusion: ∂²S/∂x² ≈ (S[i+1] - 2S[i] + S[i-1])/Δx²
        double diff = 0.0;
        if (m_dLongitudinalDispersion > 0.0) {
            diff = m_dLongitudinalDispersion * (m_vCrossSectionSalinity[i+1] - 2*m_vCrossSectionSalinity[i] + m_vCrossSectionSalinity[i-1]) / (dx*dx);
        }
        
        // Store Δ(A·S) for mass update
        m_vCrossSectionSalinityASt[i] = m_dTimestep * (adv + diff) * m_vCrossSectionArea[i];
    }
    
    // Zero flux at boundaries (can be adjusted based on BC type)
    m_vCrossSectionSalinityASt[0] = 0.0;
    m_vCrossSectionSalinityASt[m_nCrossSectionsNumber-1] = 0.0;

    // Update salinity conserving salt mass: S_new = (A·S_old + Δ(A·S))/A_new
    for (int i = 0; i < m_nCrossSectionsNumber; i++) {
        if (m_vCrossSectionArea[i] > DRY_AREA) {
            double salt_mass = m_vCrossSectionArea[i] * m_vCrossSectionSalinity[i];
            double delta_mass = m_vCrossSectionSalinityASt[i];
            salt_mass += delta_mass;
            
            // Enforce physical bounds: salinity >= 0
            if (salt_mass < 0.0) salt_mass = 0.0;
            
            m_vPredictedCrossSectionS[i] = salt_mass / m_vCrossSectionArea[i];
            
            // Clamp to physical range [0, S_ocean]
            if (m_vPredictedCrossSectionS[i] < 0.0) m_vPredictedCrossSectionS[i] = 0.0;
            if (m_vPredictedCrossSectionS[i] > m_dDownwardSalinityBoundaryValue) {
                m_vPredictedCrossSectionS[i] = m_dDownwardSalinityBoundaryValue;
            }
        } else {
            // Dry area: zero salinity
            m_vPredictedCrossSectionS[i] = 0.0;
        }
    }
    
    // Apply boundary conditions
    // Upstream salinity boundary condition
    double up_sal = m_dUpwardSalinityBoundaryValue;
    
    // If time series provided, interpolate at current time
    if (!m_strUpwardSalinityBoundaryConditionFilename.empty() && 
        !m_vUpwardSalinityBoundaryConditionTime.empty() && 
        !m_vUpwardSalinityBoundaryConditionValue.empty()) {
        double t = m_dCurrentTime;
        auto it = std::lower_bound(m_vUpwardSalinityBoundaryConditionTime.begin(), 
                                   m_vUpwardSalinityBoundaryConditionTime.end(), t);
        
        if (it == m_vUpwardSalinityBoundaryConditionTime.begin()) {
            up_sal = m_vUpwardSalinityBoundaryConditionValue.front();
        } else if (it == m_vUpwardSalinityBoundaryConditionTime.end()) {
            up_sal = m_vUpwardSalinityBoundaryConditionValue.back();
        } else {
            // Linear interpolation between time points
            size_t idx = std::distance(m_vUpwardSalinityBoundaryConditionTime.begin(), it);
            double t1 = m_vUpwardSalinityBoundaryConditionTime[idx-1];
            double t2 = m_vUpwardSalinityBoundaryConditionTime[idx];
            double s1 = m_vUpwardSalinityBoundaryConditionValue[idx-1];
            double s2 = m_vUpwardSalinityBoundaryConditionValue[idx];
            up_sal = s1 + (s2-s1)*(t-t1)/(t2-t1);
        }
    }
    
    // Apply upstream BC based on condition type
    if (nGetUpwardSalinityCondition() == 0) {
        // Type 0: Only apply if inflow (Q > 0)
        if (m_vCrossSectionQ[0] > 0.0) {
            m_vPredictedCrossSectionS[0] = up_sal;
        }
    }
    else if (nGetUpwardSalinityCondition() == 1 || nGetUpwardSalinityCondition() == 2) {
        // Type 1/2: Always impose prescribed value
        m_vPredictedCrossSectionS[0] = up_sal;
    }
    // Downstream salinity boundary condition
    double down_sal = m_dDownwardSalinityBoundaryValue;
    
    // If time series provided, interpolate at current time
    if (!m_strDownwardSalinityBoundaryConditionFilename.empty() && 
        !m_vDownwardSalinityBoundaryConditionTime.empty() && 
        !m_vDownwardSalinityBoundaryConditionValue.empty()) {
        double t = m_dCurrentTime;
        auto it = std::lower_bound(m_vDownwardSalinityBoundaryConditionTime.begin(), 
                                   m_vDownwardSalinityBoundaryConditionTime.end(), t);
        
        if (it == m_vDownwardSalinityBoundaryConditionTime.begin()) {
            down_sal = m_vDownwardSalinityBoundaryConditionValue.front();
        } else if (it == m_vDownwardSalinityBoundaryConditionTime.end()) {
            down_sal = m_vDownwardSalinityBoundaryConditionValue.back();
        } else {
            size_t idx = std::distance(m_vDownwardSalinityBoundaryConditionTime.begin(), it);
            double t1 = m_vDownwardSalinityBoundaryConditionTime[idx-1];
            double t2 = m_vDownwardSalinityBoundaryConditionTime[idx];
            double s1 = m_vDownwardSalinityBoundaryConditionValue[idx-1];
            double s2 = m_vDownwardSalinityBoundaryConditionValue[idx];
            down_sal = s1 + (s2-s1)*(t-t1)/(t2-t1);
        }
    }
    
    // Apply downstream BC based on condition type
    if (nGetDownwardSalinityCondition() == 1) {
        // Type 1: Use upstream value (unusual but allowed)
        m_vPredictedCrossSectionS[m_nCrossSectionsNumber-1] = up_sal;
    }
    else if (nGetDownwardSalinityCondition() == 2) {
        // Type 2: Use prescribed downstream value (typical for ocean boundary)
        m_vPredictedCrossSectionS[m_nCrossSectionsNumber-1] = down_sal;
    }
}

//======================================================================================================================
/**
 * @brief Corrector step for salinity transport equation
 * 
 * Uses predicted salinity from predictor step to compute corrected values.
 * Applies same advection-diffusion scheme as predictor but with predicted state.
 * Part of McCormack predictor-corrector method for second-order accuracy.
 */
//======================================================================================================================
void CSimulation::calculate_salinity_corrector() {
    // Compute advection and diffusion using predicted state
    for (int i = 1; i < m_nCrossSectionsNumber-1; ++i) {
        double dx = m_vCrossSectionDX[i];
        double u = m_vCrossSectionU[i];
        
        // Upwind advection with predicted salinity
        double adv = -u * (m_vPredictedCrossSectionS[i] - m_vPredictedCrossSectionS[i-1]) / dx;
        if (u < 0) adv = -u * (m_vPredictedCrossSectionS[i+1] - m_vPredictedCrossSectionS[i]) / dx;
        
        // Centered diffusion with predicted salinity
        double diff = 0.0;
        if (m_dLongitudinalDispersion > 0.0) {
            diff = m_dLongitudinalDispersion * (m_vPredictedCrossSectionS[i+1] - 2*m_vPredictedCrossSectionS[i] + m_vPredictedCrossSectionS[i-1]) / (dx*dx);
        }
        
        m_vCrossSectionSalinityASt[i] = m_dTimestep * (adv + diff) * m_vPredictedCrossSectionArea[i];
    }
    
    m_vCrossSectionSalinityASt[0] = 0.0;
    m_vCrossSectionSalinityASt[m_nCrossSectionsNumber-1] = 0.0;

    for (int i = 0; i < m_nCrossSectionsNumber; i++)
    {
        if (m_vCrossSectionArea[i] > DRY_AREA) {
            double salt_mass = m_vPredictedCrossSectionArea[i] * m_vPredictedCrossSectionS[i];
            double delta_mass = m_vCrossSectionSalinityASt[i];
            salt_mass += delta_mass;
            if (salt_mass < 0.0) salt_mass = 0.0;
            m_vCorrectedCrossSectionS[i] = salt_mass / m_vCrossSectionArea[i];
            if (m_vCorrectedCrossSectionS[i] < 0.0) m_vCorrectedCrossSectionS[i] = 0.0;
            if (m_vCorrectedCrossSectionS[i] > m_dDownwardSalinityBoundaryValue) m_vCorrectedCrossSectionS[i] = m_dDownwardSalinityBoundaryValue;
        } else {
            m_vCorrectedCrossSectionS[i] = 0.0;
        }
    }
    
    // Condiciones de frontera
    double up_sal = m_dUpwardSalinityBoundaryValue;
    if (!m_strUpwardSalinityBoundaryConditionFilename.empty() && !m_vUpwardSalinityBoundaryConditionTime.empty() && !m_vUpwardSalinityBoundaryConditionValue.empty()) {
        double t = m_dCurrentTime;
        auto it = std::lower_bound(m_vUpwardSalinityBoundaryConditionTime.begin(), m_vUpwardSalinityBoundaryConditionTime.end(), t);
        if (it == m_vUpwardSalinityBoundaryConditionTime.begin()) {
            up_sal = m_vUpwardSalinityBoundaryConditionValue.front();
        } else if (it == m_vUpwardSalinityBoundaryConditionTime.end()) {
            up_sal = m_vUpwardSalinityBoundaryConditionValue.back();
        } else {
            size_t idx = std::distance(m_vUpwardSalinityBoundaryConditionTime.begin(), it);
            double t1 = m_vUpwardSalinityBoundaryConditionTime[idx-1];
            double t2 = m_vUpwardSalinityBoundaryConditionTime[idx];
            double s1 = m_vUpwardSalinityBoundaryConditionValue[idx-1];
            double s2 = m_vUpwardSalinityBoundaryConditionValue[idx];
            up_sal = s1 + (s2-s1)*(t-t1)/(t2-t1);
        }
    }
    if (nGetUpwardSalinityCondition() == 0) {
        if (m_vCrossSectionQ[0] > 0.0) {
            m_vCorrectedCrossSectionS[0] = up_sal;
        }
    }
    else if (nGetUpwardSalinityCondition() == 1) {
        m_vCorrectedCrossSectionS[0] = up_sal;
    }
    else if (nGetUpwardSalinityCondition() == 2) {
        m_vCorrectedCrossSectionS[0] = up_sal;
    }
    double down_sal = m_dDownwardSalinityBoundaryValue;
    if (!m_strDownwardSalinityBoundaryConditionFilename.empty() && !m_vDownwardSalinityBoundaryConditionTime.empty() && !m_vDownwardSalinityBoundaryConditionValue.empty()) {
        double t = m_dCurrentTime;
        auto it = std::lower_bound(m_vDownwardSalinityBoundaryConditionTime.begin(), m_vDownwardSalinityBoundaryConditionTime.end(), t);
        if (it == m_vDownwardSalinityBoundaryConditionTime.begin()) {
            down_sal = m_vDownwardSalinityBoundaryConditionValue.front();
        } else if (it == m_vDownwardSalinityBoundaryConditionTime.end()) {
            down_sal = m_vDownwardSalinityBoundaryConditionValue.back();
        } else {
            size_t idx = std::distance(m_vDownwardSalinityBoundaryConditionTime.begin(), it);
            double t1 = m_vDownwardSalinityBoundaryConditionTime[idx-1];
            double t2 = m_vDownwardSalinityBoundaryConditionTime[idx];
            double s1 = m_vDownwardSalinityBoundaryConditionValue[idx-1];
            double s2 = m_vDownwardSalinityBoundaryConditionValue[idx];
            down_sal = s1 + (s2-s1)*(t-t1)/(t2-t1);
        }
    }
    if (nGetDownwardSalinityCondition() == 1) {
        m_vCorrectedCrossSectionS[m_nCrossSectionsNumber-1] = up_sal;
    }
    else if (nGetDownwardSalinityCondition() == 2) {
        m_vCorrectedCrossSectionS[m_nCrossSectionsNumber-1] = down_sal;
    }
}

//======================================================================================================================
/**
 * @brief Predictor step for temperature transport equation
 * 
 * Solves advection-diffusion equation with heat flux source:
 * ∂T/∂t + u·∂T/∂x = Kh·∂²T/∂x² + Q_net/(ρ·Cp·h)
 * 
 * where:
 * - ρ = 1000 kg/m³ (water density)
 * - Cp = 4186 J/(kg·°C) (specific heat capacity of water)
 * - Q_net = net surface heat flux (W/m²) stored in m_vCrossSectionTemperatureASt
 * 
 * Uses upwind advection and centered diffusion.
 * Applies prescribed boundary conditions or extrapolation at upstream/downstream.
 */
//======================================================================================================================
void CSimulation::calculate_temperature_predictor() {
    const int N = m_nCrossSectionsNumber;
    if (m_vPredictedCrossSectionT.size() != static_cast<size_t>(N))
        m_vPredictedCrossSectionT.resize(N, 0.0);
    
    for (int i = 1; i < N-1; ++i) {
        double adv = 0.0;
        double u = m_vCrossSectionU[i];
        double dx = m_vCrossSectionDX[i];
        
        // Upwind advection
        if (u > 0) {
            adv = -u * (m_vCrossSectionTemperature[i] - m_vCrossSectionTemperature[i-1]) / dx;
        } else {
            adv = -u * (m_vCrossSectionTemperature[i+1] - m_vCrossSectionTemperature[i]) / dx;
        }
        
        // Centered diffusion
        double diff = m_dThermalDispersion * (m_vCrossSectionTemperature[i+1] - 2*m_vCrossSectionTemperature[i] + m_vCrossSectionTemperature[i-1]) / (dx*dx);
        
        // Surface heat flux contribution: Q_net/(ρ·Cp·h)
        double Qnet = m_vCrossSectionTemperatureASt[i];
        double dTdt = adv + diff + Qnet / (WATER_DENSITY * WATER_SPECIFIC_HEAT * m_vCrossSectionWaterDepth[i]);
        
        m_vPredictedCrossSectionT[i] = m_vCrossSectionTemperature[i] + m_dTimestep * dTdt;
    }
    
    // Boundary conditions
    // Upstream temperature boundary
    if (m_nUpwardTemperatureCondition == 1) {
        // Type 1: Constant prescribed value
        m_vPredictedCrossSectionT[0] = m_dUpwardTemperatureBoundaryValue;
    } else if (m_nUpwardTemperatureCondition == 2 && !m_vUpwardTemperatureBoundaryConditionTime.empty()) {
        // Type 2: Time-varying prescribed value
        m_vPredictedCrossSectionT[0] = linearInterpolation1d(m_dCurrentTime, 
                                                             m_vUpwardTemperatureBoundaryConditionTime, 
                                                             m_vUpwardTemperatureBoundaryConditionValue);
    } else {
        // Default: Extrapolate from interior
        m_vPredictedCrossSectionT[0] = m_vPredictedCrossSectionT[1];
    }
    
    // Downstream temperature boundary
    if (m_nDownwardTemperatureCondition == 1) {
        m_vPredictedCrossSectionT[N-1] = m_dDownwardTemperatureBoundaryValue;
    } else if (m_nDownwardTemperatureCondition == 2 && !m_vDownwardTemperatureBoundaryConditionTime.empty()) {
        m_vPredictedCrossSectionT[N-1] = linearInterpolation1d(m_dCurrentTime, 
                                                               m_vDownwardTemperatureBoundaryConditionTime, 
                                                               m_vDownwardTemperatureBoundaryConditionValue);
    } else {
        m_vPredictedCrossSectionT[N-1] = m_vPredictedCrossSectionT[N-2];
    }
}

//======================================================================================================================
/**
 * @brief Corrector step for temperature transport equation
 * 
 * Uses predicted temperature from predictor step to compute corrected values.
 * Applies same advection-diffusion-heating scheme as predictor but with predicted state.
 */
//======================================================================================================================
void CSimulation::calculate_temperature_corrector() {
    const int N = m_nCrossSectionsNumber;
    if (m_vCorrectedCrossSectionT.size() != static_cast<size_t>(N))
        m_vCorrectedCrossSectionT.resize(N, 0.0);
    
    for (int i = 1; i < N-1; ++i) {
        double adv = 0.0;
        double u = m_vCrossSectionU[i];
        double dx = m_vCrossSectionDX[i];
        if (u > 0) {
            adv = -u * (m_vPredictedCrossSectionT[i] - m_vPredictedCrossSectionT[i-1]) / dx;
        } else {
            adv = -u * (m_vPredictedCrossSectionT[i+1] - m_vPredictedCrossSectionT[i]) / dx;
        }
        double diff = m_dThermalDispersion * (m_vPredictedCrossSectionT[i+1] - 2*m_vPredictedCrossSectionT[i] + m_vPredictedCrossSectionT[i-1]) / (dx*dx);
        double Qnet = m_vCrossSectionTemperatureASt[i];
        double dTdt = adv + diff + Qnet / (WATER_DENSITY * WATER_SPECIFIC_HEAT * m_vCrossSectionWaterDepth[i]);
        m_vCorrectedCrossSectionT[i] = m_vPredictedCrossSectionT[i] + m_dTimestep * dTdt;
    }
    // Fronteras
    if (m_nUpwardTemperatureCondition == 1) {
        m_vCorrectedCrossSectionT[0] = m_dUpwardTemperatureBoundaryValue;
    } else if (m_nUpwardTemperatureCondition == 2 && !m_vUpwardTemperatureBoundaryConditionTime.empty()) {
        m_vCorrectedCrossSectionT[0] = linearInterpolation1d(m_dCurrentTime, m_vUpwardTemperatureBoundaryConditionTime, m_vUpwardTemperatureBoundaryConditionValue);
    } else {
        m_vCorrectedCrossSectionT[0] = m_vCorrectedCrossSectionT[1];
    }
    if (m_nDownwardTemperatureCondition == 1) {
        m_vCorrectedCrossSectionT[N-1] = m_dDownwardTemperatureBoundaryValue;
    } else if (m_nDownwardTemperatureCondition == 2 && !m_vDownwardTemperatureBoundaryConditionTime.empty()) {
        m_vCorrectedCrossSectionT[N-1] = linearInterpolation1d(m_dCurrentTime, m_vDownwardTemperatureBoundaryConditionTime, m_vDownwardTemperatureBoundaryConditionValue);
    } else {
        m_vCorrectedCrossSectionT[N-1] = m_vCorrectedCrossSectionT[N-2];
    }
}

//======================================================================================================================
/**
 * @brief Update reservoir temperature using 0D (tank) energy balance model with full radiative physics
 * 
 * Well-mixed tank model with comprehensive surface heat flux:
 * ρ·Cp·V·dT/dt = A_surface·Q_net + ρ·Cp·Q_in·(T_in - T_water)
 * 
 * Where Q_net includes:
 * - Solar radiation with dynamic Fresnel albedo
 * - Longwave radiation (atmospheric emission and water back-radiation)
 * - Sensible heat flux (turbulent heat transfer)
 * - Latent heat flux (evaporation/condensation)
 * 
 * Uses same physics as calculateRadiativeFluxes() for consistency.
 */
//======================================================================================================================
void CSimulation::updateReservoirTemperature0D() {
    // Get calibration parameters from YAML configuration
    double latitude = m_dHeatFluxLatitude;
    double cloud_cover = m_dHeatFluxCloudCover;
    double C_h = m_dHeatFlux_CS;
    double C_e = m_dHeatFlux_CL;
    
    for (size_t i = 0; i < m_vUpwardTemperatureBoundaryConditionValue.size(); ++i) {   
        double dt = ((i == 0) ? 0.0 : m_vUpwardTemperatureBoundaryConditionTime[i] - m_vUpwardTemperatureBoundaryConditionTime[i-1]);
        
        // === METEOROLOGICAL FORCING ===
        double Tair = (i < m_vHeatFluxAirTemp.size()) ? m_vHeatFluxAirTemp[i] : 15.0;
        double rh = (i < m_vHeatFluxRelHumidity.size()) ? m_vHeatFluxRelHumidity[i] : 70.0;
        double wind = (i < m_vHeatFluxWind.size()) ? m_vHeatFluxWind[i] : 1.0;
        double pressure = (i < m_vHeatFluxAtmosphericPressure.size()) ? 
            m_vHeatFluxAtmosphericPressure[i] : ATM_PRESSURE_DEFAULT;
        
        if (wind < 0.1) wind = 0.1;  // Prevent stagnation
        
        double Twater = m_vUpwardTemperatureBoundaryConditionValue[i];
        double Tin = Tair + m_dUpwardTemperatureOffsetBeta;  // Inflow temperature
        
        // === SOLAR GEOMETRY ===
        double sim_time_hours = m_vUpwardTemperatureBoundaryConditionTime[i] / 3600.0;
        double hour_of_day = std::fmod(sim_time_hours, 24.0);
        int day_of_year = 1 + static_cast<int>(sim_time_hours / 24.0);
        
        double delta = 23.45 * std::sin(DEG_TO_RAD * (360.0/365.0) * (284.0 + day_of_year));
        double delta_rad = delta * DEG_TO_RAD;
        double lat_rad = latitude * DEG_TO_RAD;
        double omega = 15.0 * (hour_of_day - 12.0);
        double omega_rad = omega * DEG_TO_RAD;
        
        double cos_theta = std::sin(lat_rad) * std::sin(delta_rad) + 
                          std::cos(lat_rad) * std::cos(delta_rad) * std::cos(omega_rad);
        cos_theta = std::max(0.0, cos_theta);
        
        double zenith_rad = std::acos(cos_theta);
        double albedo = calc_albedo_briegleb(zenith_rad, cloud_cover);
        
        // === SHORTWAVE RADIATION ===
        double cloud_factor = 1.0 - 0.65 * cloud_cover * cloud_cover;
        double R_clear = SOLAR_CONSTANT * cos_theta * ATM_TRANSMISSIVITY;
        double Qsw_in = R_clear * cloud_factor;
        double Qsw_net = Qsw_in * (1.0 - albedo);
        
        // === LONGWAVE RADIATION ===
        double Tair_K = Tair + 273.15;
        double Twater_K = Twater + 273.15;
        double epsilon_sky = 9.37e-6 * (Tair_K * Tair_K) * (1.0 + 0.17 * cloud_cover * cloud_cover);
        double Qlw_in = epsilon_sky * STEFAN_BOLTZMANN * std::pow(Tair_K, 4);
        double Qlw_out = WATER_EMISSIVITY * STEFAN_BOLTZMANN * std::pow(Twater_K, 4);
        double Qlw_net = Qlw_in - Qlw_out;
        
        // === TURBULENT FLUXES ===
        double rho_a = calc_rho_air(Tair, pressure);
        double L_v = calc_Lv(Twater);
        
        double e_sat_air = 611.2 * std::exp(17.67 * Tair / (Tair + 243.5));
        double e_air = (rh / 100.0) * e_sat_air;
        double q_air = 0.622 * e_air / (pressure - 0.378 * e_air);
        
        double e_sat_water = 611.2 * std::exp(17.67 * Twater / (Twater + 243.5));
        double qsat_water = 0.622 * e_sat_water / (pressure - 0.378 * e_sat_water);
        
        double Qsensible = rho_a * AIR_SPECIFIC_HEAT * C_h * wind * (Tair - Twater);
        double Qlatente = rho_a * L_v * C_e * wind * (q_air - qsat_water);
        
        // === NET SURFACE HEAT FLUX ===
        double Qnet = Qsw_net + Qlw_net + Qsensible + Qlatente;  // W/m²
        
        // === TEMPERATURE TENDENCY ===
        // Surface heat flux term: Q_net / (ρ·Cp·h_eff)
        // Inflow term: (Q_in/V) * (T_in - T_water)
        double h_effective = 1.0;  // Effective depth for 0D model (can be calibrated)
        double heat_flux_term = Qnet / (WATER_DENSITY * WATER_SPECIFIC_HEAT * h_effective);
        double inflow_term = (i < m_vCrossSectionQ.size()) ? 
            m_dUpwardInflowWaterEffectkQ * m_vCrossSectionQ[i] * (Tin - Twater) : 0.0;
        
        double dTdt = heat_flux_term + inflow_term;
        
        // Update temperature
        if (i == 0) {
            m_vUpwardTemperatureBoundaryConditionValue[i] += dt * dTdt;
        } else {
            m_vUpwardTemperatureBoundaryConditionValue[i] = 
                m_vUpwardTemperatureBoundaryConditionValue[i-1] + dt * dTdt;
        }
    }
}

//===============================================================================================================================
//! Merge Predictor and Corrector
//======================================================================================================================
void CSimulation::mergePredictorCorrector() {
    //! Merge predictor and corrector steps using McCormack averaging
    //! Applies TVD flux limiters to maintain monotonicity and avoid spurious oscillations
    if (bGetDoMcComarckLimiterFlux())
    {
        //! Include TVD-McComarck? - Usar vectores pre-alocados (miembros de clase)
        // Note: Vectors don't need zeroing - overwritten in loops below
        
        // References to maintain compatibility with existing code
        auto& a1_med = m_vTVD_a1_med;
        auto& a2_med = m_vTVD_a2_med;
        auto& alfa1_med = m_vTVD_alfa1_med;
        auto& alfa2_med = m_vTVD_alfa2_med;
        auto& psi1_med = m_vTVD_psi1_med;
        auto& psi2_med = m_vTVD_psi2_med;
        auto& r1_med = m_vTVD_r1_med;
        auto& r2_med = m_vTVD_r2_med;
        auto& fi1_med = m_vTVD_fi1_med;
        auto& fi2_med = m_vTVD_fi2_med;
        auto& vFactor1 = m_vTVD_Factor1;
        auto& vFactor2 = m_vTVD_Factor2;


        //! García-Navarro Psi formula
        double delta1 = m_dDeltaValue;
        double delta2 = m_dDeltaValue;
        for (int i = 0; i < m_nCrossSectionsNumber-1; i++)
        {
            const double u_med = (m_vCrossSectionQ[i+1] / sqrt(m_vCrossSectionArea[i+1]) + m_vCrossSectionQ[i] /
                                  sqrt(m_vCrossSectionArea[i])) / (
                                     sqrt(m_vCrossSectionArea[i+1]) + sqrt(m_vCrossSectionArea[i]));
            const double A_med = (m_vCrossSectionArea[i+1] + m_vCrossSectionArea[i]) / 2.0;
            const double c_med = (m_vCrossSectionC[i+1] + m_vCrossSectionC[i]) / 2.0;
            a1_med[i] = u_med + c_med;
            a2_med[i] = u_med - c_med;

            if (!bGetDoSurfaceGradientMethod()) {
                alfa1_med[i] = ((m_vCrossSectionQ[i+1] - m_vCrossSectionQ[i]) - a2_med[i]*(A_med - m_vCrossSectionArea[i]))/(2.0*c_med);
                alfa2_med[i] = -((m_vCrossSectionQ[i+1] - m_vCrossSectionQ[i]) - a1_med[i]*(A_med - m_vCrossSectionArea[i]))/(2.0*c_med);
            }
            else {
                alfa1_med[i] = m_vCrossSectionWidth[i]*((m_vCrossSectionQ[i+1]/m_vCrossSectionWidth[i+1] - m_vCrossSectionQ[i]/m_vCrossSectionWidth[i])-a2_med[i]*(m_vCrossSectionWaterDepth[i+1] - m_vCrossSectionWaterDepth[i]))/(2.0*c_med);
                alfa2_med[i] = -m_vCrossSectionWidth[i]*((m_vCrossSectionQ[i+1]/m_vCrossSectionWidth[i+1] - m_vCrossSectionQ[i]/m_vCrossSectionWidth[i])-a1_med[i]*(m_vCrossSectionWaterDepth[i+1] - m_vCrossSectionWaterDepth[i]))/(2.0*c_med);
            }

            if (nGetPsiFormula() != 1) {
                //! Tseng Psi formula
                vector<double> deltaValues = {0.0, a1_med[i] - (m_vCrossSectionU[i] + m_vCrossSectionC[i]), m_vCrossSectionU[i+1] + m_vCrossSectionC[i+1] - a1_med[i]};
                delta1 = dMaxVectorValue(deltaValues);

                deltaValues ={0.0, a2_med[i] - (m_vCrossSectionU[i] - m_vCrossSectionC[i]), m_vCrossSectionU[i+1] + m_vCrossSectionC[i+1] - a2_med[i]};
                delta2 = dMaxVectorValue(deltaValues);
            }

            //! Computing psi_i
            if (fabs(a1_med[i])>= delta1)
            {
                psi1_med[i] = fabs(a1_med[i]);
            }
            else
            {
                psi1_med[i] = delta1;
            }

            if (fabs(a2_med[i])>= delta2)
            {
                psi2_med[i] = fabs(a2_med[i]);
            }
            else
            {
                psi2_med[i] = delta2;
            }
        }


        // === STEP 3-4: Compute r-ratios and apply flux limiter (FUSED) ===
        // r = alfa[i±1]/alfa[i] measures solution smoothness
        // r ≈ 1: smooth → limiter inactive, r >> 1 or << 1: steep gradient → limiter active
        // Direction depends on characteristic speed sign (upwind stencil)
        
        // ⚡ OPTIMIZATION: Extract limiter type check outside hot loop
        const int limiter_type = nGetEquationLimiterFlux();
        const int n_nodes = m_nCrossSectionsNumber - 1;
        
        for (int i = 0; i < n_nodes; i++)
        {
            // === Compute r1 (upwind ratio for characteristic 1) ===
            double r1, r2;
            
            if (alfa1_med[i] == 0.0) {
                r1 = 1.0;  // No wave amplitude → assume smooth
            } else {
                // Upwind stencil based on wave direction
                if (a1_med[i] < 0 && i != n_nodes - 1) {
                    r1 = alfa1_med[i+1] / alfa1_med[i];  // Downstream gradient
                } else if (a1_med[i] > 0 && i != 0) {
                    r1 = alfa1_med[i-1] / alfa1_med[i];  // Upstream gradient
                } else {
                    r1 = 1.0;  // Boundary or zero velocity
                }
            }
            
            // === Compute r2 (upwind ratio for characteristic 2) ===
            if (alfa2_med[i] == 0.0) {
                r2 = 1.0;
            } else {
                if (a2_med[i] < 0 && i != n_nodes) {
                    r2 = alfa2_med[i+1] / alfa2_med[i];
                } else if (a2_med[i] > 0 && i != 0) {
                    r2 = alfa2_med[i-1] / alfa2_med[i];
                } else {
                    r2 = 1.0;
                }
            }
            
            // ✅ PROTECTION: Clamp extreme values that cause instability
            r1 = (std::isfinite(r1) && fabs(r1) <= 1e6) ? r1 : 1.0;
            r2 = (std::isfinite(r2) && fabs(r2) <= 1e6) ? r2 : 1.0;
            
            // === Apply flux limiter function Ψ(r) ===
            // TVD region: Ψ(1) = 1, Ψ(r) ≤ min(2, 2r) (Sweby's TVD condition)
            switch (limiter_type) {
                case 1:  // MinMod: Ψ(r) = clamp(r, 0, 1) - BRANCHLESS
                    // Most diffusive, very stable, first-order near extrema
                    fi1_med[i] = std::clamp(r1, 0.0, 1.0);
                    fi2_med[i] = std::clamp(r2, 0.0, 1.0);
                    break;
                    
                case 2:  // Roe's Superbee: Ψ(r) = max(0, min(2r,1), min(r,2))
                    // Least diffusive, sharpest shock capturing
                    fi1_med[i] = dMaxVectorValue({0.0, dMinVectorValue({2*r1, 1.}), dMinVectorValue({r1, 2.})});
                    fi2_med[i] = dMaxVectorValue({0.0, dMinVectorValue({2*r2, 1.}), dMinVectorValue({r2, 2.})});
                    break;
                    
                case 3:  // Van Leer: Ψ(r) = (r + |r|)/(1 + |r|)
                    // Good balance between accuracy and stability
                    fi1_med[i] = (fabs(r1) + r1) / (1.0 + fabs(r1));
                    fi2_med[i] = (fabs(r2) + r2) / (1.0 + fabs(r2));
                    break;
                    
                case 4:  // Van Albada: Ψ(r) = (r² + r)/(1 + r²)
                    // Very smooth, slightly more diffusive than Van Leer
                    fi1_med[i] = (r1*r1 + r1) / (1.0 + r1*r1);
                    fi2_med[i] = (r2*r2 + r2) / (1.0 + r2*r2);
                    break;
                    
                default:
                    fi1_med[i] = 0.0;  // No limiter
                    fi2_med[i] = 0.0;
            }
        }
        // === STEP 5: Compute anti-diffusive flux correction ===
        // Factor = α·Ψ·(1 - λ|a|)·(1 - φ(r))
        // - α: wave amplitude (jump size)
        // - Ψ: García-Navarro limiter (bounds correction)
        // - (1 - λ|a|): CFL-based weight
        // - (1 - φ): anti-diffusive component (φ=0 for 1st order, φ=1 for TVD)
        for (int i = 0; i < m_nCrossSectionsNumber-1; i++) {
            // Anti-diffusive flux for each characteristic
            vFactor1[i] = alfa1_med[i]*psi1_med[i]*(1-m_dLambda*fabs(a1_med[i]))*(1-fi1_med[i]);
            vFactor2[i] = alfa2_med[i]*psi2_med[i]*(1-m_dLambda*fabs(a2_med[i]))*(1-fi2_med[i]);

            // Project back to physical variables: D1 for continuity, D2 for momentum
            m_vCrossSectionD1Factor[i+1] = 0.5*(vFactor1[i] + vFactor2[i]);
            m_vCrossSectionD2Factor[i+1] = 0.5*(vFactor1[i]*a1_med[i] + vFactor2[i]*a2_med[i]);
        }
        // Extrapolate correction to boundaries (usually small effect)
        m_vCrossSectionD1Factor[0] = m_vCrossSectionD1Factor[1];
        m_vCrossSectionD2Factor[0] = m_vCrossSectionD2Factor[1];
        m_vCrossSectionD1Factor[m_nCrossSectionsNumber] = m_vCrossSectionD1Factor[m_nCrossSectionsNumber-1];
        m_vCrossSectionD2Factor[m_nCrossSectionsNumber] = m_vCrossSectionD2Factor[m_nCrossSectionsNumber-1];

        // === STEP 6: Apply TVD correction to averaged solution ===
        // U^{n+1} = 0.5*(U_pred + U_corr) + TVD_correction
        // Only interior nodes - boundaries set by BC
        for (int i = 1; i < m_nCrossSectionsNumber - 1; i++) {
            m_vCrossSectionArea[i] = 0.5 * (m_vPredictedCrossSectionArea[i] + m_vCorrectedCrossSectionArea[i]) + m_dLambda *
                (m_vCrossSectionD1Factor[i + 1] - m_vCrossSectionD1Factor[i]);
            m_vCrossSectionQ[i] = 0.5*(m_vPredictedCrossSectionQ[i] + m_vCorrectedCrossSectionQ[i]) +m_dLambda*(m_vCrossSectionD2Factor[i+1] - m_vCrossSectionD2Factor[i]);
        }
        
        //! ✅ Aplicar BC finales después del merge
        //! Las fronteras deben usar los valores corrected, no el promedio
        m_vCrossSectionArea[0] = m_vCorrectedCrossSectionArea[0];
        m_vCrossSectionQ[0] = m_vCorrectedCrossSectionQ[0];
        m_vCrossSectionArea[m_nCrossSectionsNumber-1] = m_vCorrectedCrossSectionArea[m_nCrossSectionsNumber-1];
        m_vCrossSectionQ[m_nCrossSectionsNumber-1] = m_vCorrectedCrossSectionQ[m_nCrossSectionsNumber-1];
    }
}

/**
 * @brief Merge predictor and corrector steps for passive tracer transport (salinity, temperature)
 * 
 * Uses simple arithmetic averaging for passive scalars:
 *   S^{n+1} = 0.5 * (S_pred + S_corr)
 * 
 * Unlike momentum/continuity, passive tracers don't require TVD correction because:
 * - No wave propagation (only advection-diffusion)
 * - TVD is already applied in the velocity field that transports them
 * - Overshoots prevented by explicit flux limiters in predictor/corrector
 * 
 * @note Boundaries already set in calculate_salinity/temperature_corrector
 */
void CSimulation::mergeTracerPredictorCorrector() {
    // Salinity averaging
    if (bGetDoWaterSalinity() && m_vPredictedCrossSectionS.size() == static_cast<size_t>(m_nCrossSectionsNumber) && m_vCorrectedCrossSectionS.size() == static_cast<size_t>(m_nCrossSectionsNumber)) {
        for (int i = 0; i < m_nCrossSectionsNumber; ++i) {
            m_vCrossSectionSalinity[i] = 0.5 * (m_vPredictedCrossSectionS[i] + m_vCorrectedCrossSectionS[i]);
        }
    }
    // Temperature averaging
    if (bGetDoWaterTemperature() && m_vPredictedCrossSectionT.size() == static_cast<size_t>(m_nCrossSectionsNumber) && m_vCorrectedCrossSectionT.size() == static_cast<size_t>(m_nCrossSectionsNumber)) {
        for (int i = 0; i < m_nCrossSectionsNumber; ++i) {
            m_vCrossSectionTemperature[i] = 0.5 * (m_vPredictedCrossSectionT[i] + m_vCorrectedCrossSectionT[i]);
        }
    }
}

/**
 * @brief Pre-simulation bathymetry smoothing using non-uniform 1D Laplacian diffusion
 * 
 * Applies explicit diffusion to remove high-frequency noise from bed elevation data:
 *   Z^{k+1}_i = Z^k_i + α·∇²Z
 * 
 * where the non-uniform Laplacian is:
 *   ∇²Z_i = (2/(dx_L + dx_R)) * [(Z_{i+1} - Z_i)/dx_R - (Z_i - Z_{i-1})/dx_L]
 * 
 * @note Physics:
 * - Multiple passes act as multi-step diffusion (like heat equation)
 * - Alpha controls diffusion strength: larger α → more smoothing
 * - Non-uniform stencil respects variable grid spacing
 * - Boundaries kept fixed to preserve domain extent
 * 
 * @warning Too much smoothing can:
 * - Remove important morphological features (sills, steps)
 * - Violate mass conservation if not symmetric
 * - Create artificial slopes in transition zones
 */
void CSimulation::smoothBathymetry() {

    const int n = m_nCrossSectionsNumber;
    const int num_passes = m_nBathymetrySmoothingPasses;      // Number of smoothing passes from config
    const double alpha = m_dBathymetrySmoothingAlpha;         // Diffusion coefficient from config

    std::cout << "      - Smoothing bathymetry (non-uniform Laplacian, "
              << num_passes << " passes, alpha=" << alpha << ")" << std::endl;

    std::vector<double> vZSmooth(n);

    for (int pass = 0; pass < num_passes; pass++) {
        // Keep boundary values unchanged
        vZSmooth[0]   = m_vBedZ[0];
        vZSmooth[n-1] = m_vBedZ[n-1];

        // Apply non-uniform Laplacian smoothing to interior nodes
        for (int i = 1; i < n-1; i++) {
            double dxL = m_vCrossSectionX[i] - m_vCrossSectionX[i-1];
            double dxR = m_vCrossSectionX[i+1] - m_vCrossSectionX[i];

            // Non-uniform 1D Laplacian
            double lap =
                (2.0 / (dxL + dxR)) *
                ( (m_vBedZ[i+1] - m_vBedZ[i]) / dxR
                - (m_vBedZ[i]   - m_vBedZ[i-1]) / dxL );

            // Explicit smoothing step
            vZSmooth[i] = m_vBedZ[i] + alpha * lap;
        }

        // Update values for the next pass
        for (int i = 1; i < n-1; i++) {
            m_vBedZ[i] = vZSmooth[i];
        }
    }
}

/**
 * @brief Runtime solution smoothing with adaptive 2-pass algorithm
 * 
 * Pass 1: Global Laplacian smoothing (similar to smoothBathymetry)
 * Pass 2: Targeted smoothing only at strong area jumps
 * 
 * Strong jump criteria: A[i]/A[i±1] > 1.5 or < 0.67
 * 
 * @note Purpose:
 * - Stabilize solution near sudden geometry changes (bridge piers, junctions)
 * - Remove Gibbs oscillations from sharp fronts
 * - Maintain CFL stability in regions with complex geometry
 * 
 * @warning Can introduce:
 * - Artificial diffusion (reduces effective Reynolds number)
 * - Phase errors in wave propagation
 * - Mass/momentum conservation errors if alpha too large
 * 
 * @see smoothBathymetry() for similar Laplacian formulation
 */
void CSimulation::smoothSolution() {
    const int n = m_nCrossSectionsNumber;
    const int num_passes = m_nSolutionSmoothingPasses;
    const double alpha = m_dSolutionSmoothingAlpha;

    for (int pass = 0; pass < num_passes; ++pass) {
        std::vector<double> vAreaSmooth(n);
        std::vector<double> vQSmooth(n);

        // Keep boundaries unchanged
        vAreaSmooth[0]   = m_vCrossSectionArea[0];
        vAreaSmooth[n-1] = m_vCrossSectionArea[n-1];
        vQSmooth[0]      = m_vCrossSectionQ[0];
        vQSmooth[n-1]    = m_vCrossSectionQ[n-1];

        // First pass: non-uniform Laplacian smoothing
        for (int i = 1; i < n-1; i++) {
            double dxL = m_vCrossSectionX[i] - m_vCrossSectionX[i-1];
            double dxR = m_vCrossSectionX[i+1] - m_vCrossSectionX[i];

            // Laplacian for Area
            double lapA =
                (2.0 / (dxL + dxR)) *
                ( (m_vCrossSectionArea[i+1] - m_vCrossSectionArea[i]) / dxR
                - (m_vCrossSectionArea[i]   - m_vCrossSectionArea[i-1]) / dxL );

            // Laplacian for Q
            double lapQ =
                (2.0 / (dxL + dxR)) *
                ( (m_vCrossSectionQ[i+1] - m_vCrossSectionQ[i]) / dxR
                - (m_vCrossSectionQ[i]   - m_vCrossSectionQ[i-1]) / dxL );

            vAreaSmooth[i] = m_vCrossSectionArea[i] + alpha * lapA;
            vQSmooth[i]    = m_vCrossSectionQ[i]    + alpha * lapQ;
        }

        // Second pass only in regions with strong area jumps
        for (int i = 2; i < n-2; i++) {
            double ratioL = m_vCrossSectionArea[i]   / (m_vCrossSectionArea[i-1] + 1e-10);
            double ratioR = m_vCrossSectionArea[i+1] / (m_vCrossSectionArea[i]   + 1e-10);

            bool strongJump =
                (ratioL > 1.5 || ratioL < 0.67 ||
                 ratioR > 1.5 || ratioR < 0.67);

            if (strongJump) {
                double dxL = m_vCrossSectionX[i]   - m_vCrossSectionX[i-1];
                double dxR = m_vCrossSectionX[i+1] - m_vCrossSectionX[i];

                double lapA =
                    (2.0 / (dxL + dxR)) *
                    ( (vAreaSmooth[i+1] - vAreaSmooth[i]) / dxR
                    - (vAreaSmooth[i]   - vAreaSmooth[i-1]) / dxL );

                double lapQ =
                    (2.0 / (dxL + dxR)) *
                    ( (vQSmooth[i+1] - vQSmooth[i]) / dxR
                    - (vQSmooth[i]   - vQSmooth[i-1]) / dxL );

                vAreaSmooth[i] += alpha * lapA;
                vQSmooth[i]    += alpha * lapQ;
            }
        }

        // Apply smoothed values and update velocity
        for (int i = 1; i < n-1; i++) {
            m_vCrossSectionArea[i] = vAreaSmooth[i];
            m_vCrossSectionQ[i]    = vQSmooth[i];
            m_vCrossSectionU[i]    = m_vCrossSectionQ[i] /
                                     (m_vCrossSectionArea[i] + 1e-10);
        }
    }
}

/**
 * @brief Read tides.csv and extract maximum astronomical tide level
 * 
 * Parses CSV file with format: time,elevation
 * 
 * @param tidesFile Path to tides.csv file
 * @return Maximum tide elevation (m), or 0.0 if file not found/invalid
 */
double CSimulation::getMaxAstronomicalTide(const std::string& tidesFile) {
    std::ifstream file(tidesFile);
    if (!file.is_open()) return 0.0;
    double maxTide = -1e9;
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        size_t comma = line.find(',');
        if (comma == std::string::npos) continue;
        try {
            double value = std::stod(line.substr(comma + 1));
            if (value > maxTide) maxTide = value;
        } catch (...) { continue; }
    }
    file.close();
    return maxTide;
}

/**
 * @brief Accessor for simulation state variables (for post-processing, NetCDF output, debugging)
 * 
 * Available variables:
 * - Flow: "A", "Q", "U", "c" (wave speed)
 * - Predictor/Corrector: "Ap", "Ac", "Qp", "Qc"
 * - Geometry: "B" (width), "Rh" (hydraulic radius), "level", "eta" (elevation)
 * - Transport: "S" (salinity), "T" (temperature), "rho" (density)
 * - Derivatives: "DhDx" (surface slope)
 * - Sediment: "Qb" (bedload), "Qs" (suspended), "Qt" (total)
 * - Bank locations: "xl", "xr", "xl_utm_x", "xl_utm_y", "xr_utm_x", "xr_utm_y"
 * 
 * @param strVariableName Name of variable to retrieve
 * @return Vector of values for all cross-sections, or empty vector if not found
 * 
 * @warning No bounds checking - caller must verify size matches m_nCrossSectionsNumber
 */
vector<double> CSimulation::vGetVariable(const string& strVariableName) const {

    if (strVariableName == "A") {
        return m_vCrossSectionArea;
    }
    else if (strVariableName == "Ap") {
        return m_vPredictedCrossSectionArea;
    }
    else if (strVariableName == "Ac") {
        return m_vCorrectedCrossSectionArea;
    }
    else if (strVariableName == "Q") {
        return m_vCrossSectionQ;
    }
    else if (strVariableName == "Qp") {
        return m_vPredictedCrossSectionQ;
    }
    else if (strVariableName == "Qc") {
        return m_vCorrectedCrossSectionQ;
    }
    else if (strVariableName == "Rh") {
        return m_vCrossSectionHydraulicRadius;
    }
    else if (strVariableName == "B") {
        return m_vCrossSectionWidth;
    }
    else if (strVariableName == "level") {
        return m_vCrossSectionWaterDepth;
    }
    else if (strVariableName == "eta") {
        return m_vCrossSectionWaterElevation;
    }
    else if (strVariableName == "DhDx") {
        return m_vCrossSectionDhDx;
    }
    else if (strVariableName == "rho") {
        return m_vCrossSectionDensity;
    }
    else if (strVariableName == "U") {
        return m_vCrossSectionU;
    }
    else if (strVariableName == "c") {
        return m_vCrossSectionC;
    }
    else if (strVariableName == "S") {
        return m_vCrossSectionSalinity;
    }
    else if (strVariableName == "T") {
        return m_vCrossSectionTemperature;
    }
    else if (strVariableName == "xl") {
        return m_vCrossSectionLeftRBLocation;
    }
    else if (strVariableName == "xr") {
        return m_vCrossSectionRightRBLocation;
    }
    else if (strVariableName == "xl_utm_x") {
        return m_vCrossSectionLeftRBLocation_UTM_X;
    }
    else if (strVariableName == "xl_utm_y") {
        return m_vCrossSectionLeftRBLocation_UTM_Y;
    }
    else if (strVariableName == "xr_utm_x") {
        return m_vCrossSectionRightRBLocation_UTM_X;
    }
    else if (strVariableName == "xr_utm_y") {
        return m_vCrossSectionRightRBLocation_UTM_Y;
    }
    else if (strVariableName == "Qb") {
        return m_vCrossSectionQb;
    }
    else if (strVariableName == "Qs") {
        return m_vCrossSectionQs;
    }
    else if (strVariableName == "Qt") {
        return m_vCrossSectionQt;
    }
    // else if (strVariableName == "q") {
    //     return m_vCrossSectionQt;
    // }
    else {
        return {};
    }
}

/**
 * @brief Display simulation progress with ETA and current simulation date/time
 * 
 * Prints: [Elapsed][Remaining] Time - Progress % - SimDate: YYYY-MM-DD HH:MM:SS
 * 
 * @note Performance optimization:
 * - Only updates every 100 timesteps to reduce I/O overhead
 * - I/O can be expensive (terminal flushing, string formatting)
 * - For dt=1s, 100 timesteps = ~1.7 minutes between updates
 * - Always updates when saving (m_bSaveTime=true)
 * 
 * @warning Static variables used for persistent state across calls:
 * - call_count: throttling counter
 * - t_start: cached simulation start time (computed once)
 */
void CSimulation::AnnounceProgress() {

    // ⚡ OPTIMIZATION: Reduce screen update frequency to minimize I/O overhead
    // Update every 1000 iterations (~2-3% of total time) or when saving output
    static int call_count = 0;
    if (++call_count % 3600 != 0 && !m_bSaveTime) return;

    // Stdout is connected to a tty, so not running as a background job
    static double sdElapsed = 0;
    static double sdToGo = 0;
    const time_t tNow = time(nullptr);

    // Calculate time elapsed and remaining
    sdElapsed = difftime(tNow, m_tSysStartLoopTime);
    sdToGo = (sdElapsed * m_dSimDuration / m_dCurrentTime) - sdElapsed;


    // Calcular t_start solo una vez (static)
    static time_t t_start = 0;
    if (t_start == 0) {
        std::tm tm_start = {};
        tm_start.tm_year = m_nSimStartYear - 1900;
        tm_start.tm_mon = m_nSimStartMonth - 1;
        tm_start.tm_mday = m_nSimStartDay;
        tm_start.tm_hour = m_nSimStartHour;
        tm_start.tm_min = m_nSimStartMin;
        tm_start.tm_sec = m_nSimStartSec;
        tm_start.tm_isdst = -1;
        t_start = std::mktime(&tm_start);
    }
    time_t t_current = t_start + static_cast<time_t>(m_dCurrentTime);
    std::tm* tm_current = std::localtime(&t_current);

    // Imprimir fecha y hora en formato YYYY-MM-DD HH:MM:SS
    char datetime_buf[32];
    std::strftime(datetime_buf, sizeof(datetime_buf), "%Y-%m-%d %H:%M:%S", tm_current);

    cout << "\r    - Elapsed[Remaining] Time: " << std::fixed << setprecision(3) << setw(6) << sdElapsed <<"[" << std::fixed << setprecision(3) << setw(6) << sdToGo << "] s - Progress: " << std::fixed << setprecision(3) << setw(6) << 100 * m_dCurrentTime / m_dSimDuration  << "% - SimDate: " << datetime_buf << std::flush;

}





/**
 * @brief Calculate comprehensive net surface heat flux for water temperature equation
 * 
 * Computes Q_net = Q_SW + Q_LW + Q_H + Q_E (W/m²) using physically-based parameterizations
 * 
 * RADIATION COMPONENTS:
 * --------------------
 * 1. Shortwave radiation (Q_SW): Solar radiation absorbed by water
 *    - Synthetic incident radiation from solar geometry (Kasten-Czeplak model)
 *    - Dynamic albedo from Briegleb formula (zenith-angle dependent, ~0.05-0.40)
 *    - Q_SW = SW↓ × (1 - α(θ_z))
 * 
 * 2. Longwave radiation (Q_LW): Net thermal radiation exchange
 *    - Incoming: εsky × σ × T_air⁴ (Swinbank formula with cloud correction)
 *    - Outgoing: εwater × σ × T_water⁴
 *    - Q_LW = LW↓ - LW↑ (typically negative, net cooling)
 * 
 * TURBULENT FLUXES (Bulk Aerodynamic Formulas):
 * ----------------------------------------------
 * 3. Sensible heat flux (Q_H): Turbulent heat conduction
 *    Q_H = ρ_air(T_air) × c_p × C_H × U × (T_air - T_water)
 *    - C_H ≈ 1.3×10⁻³: Stanton number (calibration parameter)
 *    - Positive when air warms water
 * 
 * 4. Latent heat flux (Q_E): Evaporation/condensation
 *    Q_E = ρ_air(T_air) × L_v(T_water) × C_E × U × (q_air - q_sat,water)
 *    - C_E ≈ 1.3×10⁻³: Dalton number (calibration parameter)
 *    - Negative during evaporation (typical), positive during condensation
 *    - q = specific humidity from Magnus-Tetens formula
 * 
 * PHYSICAL IMPROVEMENTS OVER PREVIOUS VERSION:
 * ---------------------------------------------
 * - Dynamic air density ρ_air(T) using ideal gas law (see calc_rho_air)
 * - Dynamic latent heat L_v(T) temperature-dependent (see calc_Lv)
 * - Solar geometry calculation for realistic diurnal cycle
 * - Zenith-angle dependent albedo (Briegleb et al. 1986 reflectance)
 * - Sky emissivity with cloud cover correction
 * - Proper sign conventions (Q_net > 0 → water warms)
 * 
 * CONFIGURATION PARAMETERS:
 * -------------------------
 * - Latitude: 36.5°N (modifiable for other locations)
 * - Cloud cover: 0.2 (20%, can be made time-dependent if data available)
 * - C_H, C_E: 1.3×10⁻³ (calibration parameters, adjust for local conditions)
 * 
 * @note Requires meteorological forcing data:
 * - m_vHeatFluxAirTemp: Air temperature (°C)
 * - m_vHeatFluxRelHumidity: Relative humidity (%)
 * - m_vHeatFluxWind: Wind speed (m/s)
 * 
 * @note Output stored in m_vCrossSectionTemperatureASt[i] for each cross-section
 * 
 * @warning Minimum wind speed clamped to 0.1 m/s to avoid numerical issues
 * 
 * @see calculateTemperatureSourceTerms() for integration of Q_net into dT/dt equation
 * @see calc_rho_air() for temperature-dependent air density
 * @see calc_Lv() for temperature-dependent latent heat
 * @see calc_albedo_briegleb() for zenith-angle dependent albedo with direct/diffuse separation
 * 
 * @ref Swinbank (1963): Long-wave radiation from clear skies
 * @ref Kasten & Czeplak (1980): Solar radiation parameterization
 * @ref Briegleb et al. (1986): Angle-dependent surface reflectance
 */
void CSimulation::calculateRadiativeFluxes() {
    // === SIMULATION CONFIGURATION ===
    // Derive current hour and day from simulation time
    double sim_time_hours = m_dCurrentTime / 3600.0; 
    double hour_of_day = std::fmod(sim_time_hours, 24.0);  // Hour 0-24
    int day_of_year = 1 + static_cast<int>(sim_time_hours / 24.0); // Day 1-365 (simplified)
    
    // Get calibration parameters from YAML configuration (or use defaults)
    double latitude = m_dHeatFluxLatitude;        // Geographic latitude (degrees, North positive)
    double cloud_cover = m_dHeatFluxCloudCover;  // Cloud cover fraction (0.0=clear, 1.0=overcast)
    double C_h = m_dHeatFlux_CS;                 // Sensible heat transfer coefficient (Stanton number)
    double C_e = m_dHeatFlux_CL;                 // Latent heat transfer coefficient (Dalton number)

    for (int i = 0; i < m_nCrossSectionsNumber; ++i) {
        // === 1. STATE VARIABLES ===
        // Meteorological forcing from input data (with defaults if unavailable)
        double Tair = (!m_vHeatFluxTime.empty()) ? 
            linearInterpolation1d(m_dCurrentTime, m_vHeatFluxTime, m_vHeatFluxAirTemp) : 15.0;
        double rh = (!m_vHeatFluxTime.empty()) ? 
            linearInterpolation1d(m_dCurrentTime, m_vHeatFluxTime, m_vHeatFluxRelHumidity) : 70.0;
        double wind = (!m_vHeatFluxTime.empty()) ? 
            linearInterpolation1d(m_dCurrentTime, m_vHeatFluxTime, m_vHeatFluxWind) : 1.0;
        double pressure = (!m_vHeatFluxAtmosphericPressure.empty()) ?
            linearInterpolation1d(m_dCurrentTime, m_vHeatFluxTime, m_vHeatFluxAtmosphericPressure) : ATM_PRESSURE_DEFAULT;
        double Twater = m_vCrossSectionTemperature[i];
        
        // Prevent zero wind speed (avoids unrealistic stagnation and division by zero)
        if (wind < 0.1) wind = 0.1;

        // Dynamic physical properties (temperature-dependent)
        double rho_a_dynamic = calc_rho_air(Tair, pressure);  // Air density (kg/m³) with actual pressure
        double L_v_dynamic = calc_Lv(Twater);                 // Latent heat of vaporization (J/kg)

        // === 2. SOLAR GEOMETRY (For Albedo and Synthetic Radiation) ===
        // Solar declination angle (degrees) - function of day of year
        // δ = 23.45° × sin[2π/365 × (284 + n)]
        double delta = 23.45 * std::sin(DEG_TO_RAD * (360.0/365.0) * (284.0 + day_of_year));
        double delta_rad = delta * DEG_TO_RAD;
        double lat_rad = latitude * DEG_TO_RAD;
        
        // Hour angle (degrees from solar noon)
        double omega = 15.0 * (hour_of_day - 12.0);
        double omega_rad = omega * DEG_TO_RAD;

        // Cosine of solar zenith angle (sun elevation)
        // cos(θ_z) = sin(φ)sin(δ) + cos(φ)cos(δ)cos(ω)
        double cos_theta = std::sin(lat_rad) * std::sin(delta_rad) + 
                           std::cos(lat_rad) * std::cos(delta_rad) * std::cos(omega_rad);
        
        double zenith_rad = std::acos(std::max(0.0, cos_theta));

        // === 3. Q_SW: NET SHORTWAVE RADIATION ===
        double albedo_dynamic = 0.06;  // Default albedo (Briegleb diffuse value for nighttime)
        double Qsw_in = 0.0;           // Incident shortwave radiation (W/m²)

        if (cos_theta > 0.0) {  // Daytime
            // Dynamic albedo from Briegleb et al. (1986) model
            albedo_dynamic = calc_albedo_briegleb(zenith_rad, cloud_cover);

            // Synthetic incident shortwave radiation (Kasten-Czeplak model)
            double R_clear = SOLAR_CONSTANT * cos_theta * ATM_TRANSMISSIVITY;
            double cloud_factor = 1.0 - 0.75 * std::pow(cloud_cover, 3.4);
            Qsw_in = R_clear * cloud_factor;
        }

        // Net shortwave radiation absorbed by water (positive = heat gain)
        double Qsw = Qsw_in * (1.0 - albedo_dynamic);

        // === 4. Q_LW: NET LONGWAVE RADIATION ===
        double Tair_K = Tair + 273.15;
        double Twater_K = Twater + 273.15;
        
        // Sky emissivity (Swinbank formula with cloud correction)
        // ε_sky = 9.37×10⁻⁶ T_air² (1 + 0.17 C²)
        double epsilon_sky = 9.37e-6 * (Tair_K * Tair_K) * 
                             (1.0 + 0.17 * cloud_cover * cloud_cover);
        
        // Incoming longwave radiation from atmosphere
        double Qlw_in = epsilon_sky * STEFAN_BOLTZMANN * std::pow(Tair_K, 4);
        
        // Outgoing longwave radiation from water surface
        double Qlw_out = WATER_EMISSIVITY * STEFAN_BOLTZMANN * std::pow(Twater_K, 4);

        // Net longwave radiation (typically negative/cooling)
        double Qlw = Qlw_in - Qlw_out;

        // === 5. TURBULENT HEAT FLUXES (Bulk Aerodynamic Formulas) ===
        
        // --- SENSIBLE HEAT FLUX (Q_H) ---
        // Q_H = ρ_air × c_p × C_h × U × (T_air - T_water)
        // Sign convention: Positive when air is warmer (water gains heat)
        double bulk_coeff_sensible = rho_a_dynamic * AIR_SPECIFIC_HEAT * C_h * wind;
        double Qsensible = bulk_coeff_sensible * (Tair - Twater); 

        // --- LATENT HEAT FLUX (Q_E) ---
        // Q_E = ρ_air × L_v × C_e × U × (q_air - q_sat,water)
        
        // Saturation vapor pressure (Magnus-Tetens formula, result in hPa)
        double esat_water = 6.112 * std::exp((17.62 * Twater) / (243.12 + Twater));
        double esat_air = 6.112 * std::exp((17.62 * Tair) / (243.12 + Tair));
        
        // Specific humidity (kg water vapor / kg air)
        double qsat_water = 0.622 * esat_water / 1013.25;  // At water surface (saturated)
        double q_air = (rh / 100.0) * 0.622 * esat_air / 1013.25;  // In air (from relative humidity)
        
        // Latent heat flux calculation
        // Sign convention: Positive when air is more humid (condensation, water gains heat)
        //                  Negative when evaporation occurs (water loses heat)
        double bulk_coeff_latent = rho_a_dynamic * L_v_dynamic * C_e * wind;
        double Qlatente = bulk_coeff_latent * (q_air - qsat_water); 

        // === 6. NET SURFACE HEAT FLUX ===
        // Q_net = Q_SW + Q_LW + Q_H + Q_E
        // Positive Q_net → water warms
        // Negative Q_net → water cools
        double Qnet = Qsw + Qlw + Qsensible + Qlatente;

        // Store result for temperature equation (used in calculateTemperatureSourceTerms)
        m_vCrossSectionTemperatureASt[i] = Qnet;
    }
}


/**
 * @brief Calculate dynamic albedo from latitude, day of year, and hour of day
 * 
 * Computes albedo based on solar position geometry:
 * 
 * 1. Solar declination (Cooper 1969):
 *    δ = 23.45°·sin[2π/365 × (284 + N)]
 *    where N = day of year (1-365)
 * 
 * 2. Hour angle (15° per hour from solar noon):
 *    ω = 15° × (hour - 12)
 * 
 * 3. Solar zenith angle:
 *    cosθ = sinφ·sinδ + cosφ·cosδ·cosω
 *    where φ = latitude
 * 
 * 4. Albedo from Briegleb et al. (1986) model (see calc_albedo_briegleb)
 * 
 * @param lat_deg Latitude in degrees (e.g., 36.5 for Cádiz)
 * @param day_of_year Day of year (1-365, 1=Jan 1)
 * @param hour_of_day Decimal hour (0-24, e.g., 14.5 = 14:30)
 * @return Albedo (0.03-1.0)
 * 
 * @note Simplifications:
 * - Assumes solar noon at 12:00 local time (neglects equation of time ~±15 min)
 * - Neglects longitude correction (use UTC + timezone offset for precision)
 * - Good accuracy for heat budget (~5% error acceptable)
 * 
 * @warning For high-precision applications:
 * - Use NOAA solar position algorithm
 * - Include atmospheric refraction
 * - Account for equation of time
 * 
 * @see calc_albedo_briegleb() for zenith-to-albedo conversion
 */
double calc_dynamic_albedo(double lat_deg, int day_of_year, double hour_of_day) {
    
    // 1. SOLAR POSITION CALCULATION (Basic geometry)
    
    // Solar declination (delta): Sun angle relative to equator as function of day of year
    // Cooper formula (1969) - Sufficient for engineering applications
    double delta = 23.45 * std::sin(CSimulation::DEG_TO_RAD * (360.0/365.0) * (284.0 + day_of_year));
    double delta_rad = delta * CSimulation::DEG_TO_RAD;

    // Hour angle (omega): 0 at solar noon, +/- for each hour
    // (15 degrees per hour)
    double time_offset = 12.0; // Assuming solar noon at 12:00 (simplification)
    double omega = 15.0 * (hour_of_day - time_offset); 
    double omega_rad = omega * CSimulation::DEG_TO_RAD;

    double lat_rad = lat_deg * CSimulation::DEG_TO_RAD;

    // Cosine of zenith angle (cos_theta)
    // theta = 0 (sun vertical), theta = 90 (horizon)
    double cos_theta = std::sin(lat_rad) * std::sin(delta_rad) + 
                       std::cos(lat_rad) * std::cos(delta_rad) * std::cos(omega_rad);

    // 2. CHECK IF NIGHTTIME
    // If sun is below horizon, albedo doesn't matter (Radiation = 0), 
    // but we return 0.07 as default to avoid errors
    if (cos_theta <= 0.0) {
        return 0.07; // Nighttime or default value
    }

    // 3. ALBEDO CALCULATION (Fresnel approximation for water)
    // Many approximations exist. This one is robust and widely used in oceanography.
    // Albedo A = 0.5 * (Parallel_reflection + Perpendicular_reflection)
    
    // We use empirical approximation (Payne 1972, Briegleb et al. 1986)
    // This formula is more efficient than complete Fresnel and very accurate:
    // Approx from "Briegleb et al. (1986)" for marine surface albedo:
    
    double num = 0.037;
    double den = (1.1 * std::pow(cos_theta, 1.4)) + 0.15;
    double albedo = num / den;

    // Limit to physical values (albedo never > 1.0)
    return std::min(std::max(albedo, 0.03), 1.0);
}

/**
 * @brief Calculate sediment transport using van Rijn (1984) formulation
 * 
 * Separates transport into:
 * 1. Bedload (Qb): Grains rolling/sliding along bed
 * 2. Suspended load (Qs): Grains carried in water column
 * 3. Total (Qt): Qb + Qs
 * 
 * BEDLOAD (van Rijn 1984):
 *   u* = (g^0.5 / C) |U|  (shear velocity, C = Chézy coefficient)
 *   θ = u*² / [(s-1)·g·D50]  (Shields parameter)
 *   T = (θ - θ_cr) / θ_cr  (transport stage parameter)
 *   
 *   For T < 3:  gb = 0.053·√[(s-1)·g·D50³]·T^2.1·D*^{-0.3}
 *   For T ≥ 3:  gb = 0.10·√[(s-1)·g·D50³]·T^1.5·D*^{-0.3}
 *   
 *   where D* = D50·[(s-1)·g/ν²]^{1/3} (dimensionless grain size)
 * 
 * SUSPENDED LOAD:
 *   Uses Rouse-Vanoni profile with reference concentration at z = Δb:
 *   C(z) = C_a·[(h-z)/z · Δb/(h-Δb)]^Z
 *   
 *   where:
 *   - Z = w_s/(β·κ·u*): Rouse number (w_s = settling velocity)
 *   - Δb = 0.3·D*^0.7·√T·D_ave: Reference height
 *   - C_a = 0.117·ρ_s·T/D*: Reference concentration
 * 
 * @note Valid for:
 * - Non-cohesive sediments (sand/gravel)
 * - D50: 0.0001 - 0.01 m
 * - Flow not dominated by vegetation drag
 * 
 * @warning Limitations:
 * - Assumes steady, uniform flow (not ideal for tides)
 * - No bed slope correction (can be important in estuaries)
 * - No turbulence damping by stratification
 * - Calibration needed for site-specific conditions
 * 
 * @see van Rijn, L.C. (1984). Sediment Transport, Part I: Bed Load Transport.
 *      Journal of Hydraulic Engineering, 110(10), 1431-1456.
 */
void CSimulation::calculate_sediment_transport()
{
    for (int i=0; i< m_nCrossSectionsNumber; i++)
    {
        // The sign should be included at the end of the calculation
        double vel_abs = fabs(m_vCrossSectionU[i]);

        // Sediment transport direction (+ve upward, -ne downward)
        double dSedimentDirection = 0.0;
        if (m_vCrossSectionU[i] > 0) dSedimentDirection = 1;
        else if (m_vCrossSectionU[i] < 0) dSedimentDirection = -1;
        else dSedimentDirection = 0;


        //====================================================================================================
        //! Bed load sediment transport
        //====================================================================================================
        double c_1 = 18*log(12*m_vCrossSectionHydraulicRadius[i]/(3*m_vCrossSectionD90[i]));
        if (c_1 < 1e-3) c_1 = 1e-3;

        double u_star = pow(G, 0.5)/c_1*vel_abs;
        double shields_crit = 0.0013*pow(m_vCrossSectionDiamX[i], 0.29);
        if (m_vCrossSectionDiamX[i] < 150) shields_crit = 0.055;

        double u_star_crit = sqrt(shields_crit*(m_vCrossSectionRhos[i] - 1)*G*m_vCrossSectionD50[i]);
        //! Sediment transport
        double transport = (u_star*u_star - u_star_crit*u_star_crit)/(u_star_crit*u_star_crit);

        // If sediment transport is negative means no sediment transport
        if (transport < 0) transport = 0;

        // Mass sediment transport
        double gb = 0.0;
        if (transport >= 3)
        {
            gb = 0.1*sqrt((m_vCrossSectionRhos[i]-1)*G*pow(m_vCrossSectionD50[i], 3.0))*pow(transport, 1.5)*pow(m_vCrossSectionDiamX[i], -0.3);
        }
        else
        {
            gb = 0.053*sqrt((m_vCrossSectionRhos[i]-1)*G*pow(m_vCrossSectionD50[i], 3.0))*pow(transport, 2.1)*pow(m_vCrossSectionDiamX[i], -0.3);
        }
        //! Volumetric sediment transport
        m_vCrossSectionQb[i] = dSedimentDirection*gb*m_vCrossSectionWidth[i]*FRESH_WATER_DENSITY*m_vCrossSectionRhos[i]*dSedimentDirection;

        //====================================================================================================
        //! Suspended sediment transport
        //====================================================================================================
        //! Jump height
        double dDeltaB = m_vCrossSectionDaveraged[i]*0.3*pow(m_vCrossSectionDiamX[i], 0.7)*sqrt(transport);
        // Correction of delta_b
        if (dDeltaB < 0.01*m_vCrossSectionWaterDepth[i]) dDeltaB = 0.01*m_vCrossSectionWaterDepth[i];
        if (dDeltaB > 0.5*m_vCrossSectionWaterDepth[i]) dDeltaB = 0.5*m_vCrossSectionWaterDepth[i];

        //! Shear velocity
        double dUx = sqrt(G*m_vCrossSectionHydraulicRadius[i]*fabs(m_vCrossSectionFrictionSlope[i]));

        //! Reference concentration at z=delta_b with plane bottom
        double c_a = 0.117*FRESH_WATER_DENSITY*m_vCrossSectionRhos[i]*transport/m_vCrossSectionDiamX[i];

        //! Representative diameter of suspended particle
        double dRepresentativeDiameter = m_vCrossSectionD50[i]*(1.0+0.11*(m_vCrossSectionSedimentSigma[i]-1.0));

        if (transport >= 25) dRepresentativeDiameter = m_vCrossSectionD50[i];

        double dSettlingVelocity = 0.0;
        //! Settling velocity of representative diameters
        if (dRepresentativeDiameter < 0.0001)
        {
            dSettlingVelocity = (m_vCrossSectionRhos[i] - 1.0)*G*pow(dRepresentativeDiameter, 2.0)/(18.0*NU);
        }
        else if ((dRepresentativeDiameter >= 0.0001) & (dRepresentativeDiameter < 0.001))
        {
            dSettlingVelocity = 10.0*NU*sqrt((1.0 + 0.01*G*(m_vCrossSectionRhos[i] - 1.00)*pow(dRepresentativeDiameter, 3.0)/(NU*NU)) - 1.0)/dRepresentativeDiameter;
        }
        else
        {
            dSettlingVelocity = 1.1*sqrt((m_vCrossSectionRhos[i] - 1.0)*G*dRepresentativeDiameter);
        }

        //! Factor that considers the different diffusion between fluid and sediment particles
        double dBeta = 2.0;

        if (dUx != 0.0) dBeta = 1.0 + 2.0*pow(dSettlingVelocity/dUx, 2.0);
        if (dBeta > 2.0) dBeta = 2.0;

        //! Rouse suspension parameter
        double dRouse = dSettlingVelocity/(dBeta*KAPPA*dUx);

        //! Maximum concentration of sediments at the reference level
        double C0 = 0.65*FRESH_WATER_DENSITY*m_vCrossSectionRhos[i];

        //! Global correction factor
        double dPsi = 2.5*pow(dSettlingVelocity/dUx, 0.8)*pow(c_a/C0, 0.4);

        //! Update Rouse number
        dRouse = dRouse + dPsi;

        //! Adimensional factor A and F
        double dA_Factor = dDeltaB/m_vCrossSectionWaterDepth[i];
        double dF_Factor = 0.0;
        if (fabs(pow(1.0 - dA_Factor, dRouse)*(1.2 + dRouse)) >= 1e-4)
        {
            dF_Factor = (pow(dA_Factor, dRouse) - pow(dA_Factor, 1.2))/(pow(1.0-dA_Factor, dRouse)*(1.2 - dRouse));
        }

        double gbs1 = c_a*m_vCrossSectionWaterDepth[i]*vel_abs*dF_Factor;
        //! Suspended sediment transport
        m_vCrossSectionQs[i] = gbs1/(m_vCrossSectionRhos[i]*FRESH_WATER_DENSITY)*m_vCrossSectionWidth[i]*dSedimentDirection;

        //! Total sediment transport
        m_vCrossSectionQt[i] = m_vCrossSectionQb[i] + m_vCrossSectionQs[i];
        if (m_vCrossSectionQt[i] < 1e-6) m_vCrossSectionQt[i] = 0.0;

        return;
    }

};


/**
 * @brief Calculate water density from salinity, temperature, and sediment concentration
 * 
 * Linear equation of state (EOS):
 *   ρ = ρ₀·[1 + β_S·S - β_T·(T - 4)]
 * 
 * where:
 * - ρ₀ = 1000 kg/m³: Reference density (freshwater at 4°C)
 * - β_S ≈ 0.00078: Haline contraction coefficient (kg/m³/psu)
 * - β_T ≈ 0.0002: Thermal expansion coefficient (1/°C)
 * - S: Salinity (psu = g/kg)
 * - T: Temperature (°C)
 * 
 * Sediment effect (if enabled):
 *   ρ_total = ρ_water + (1 - ρ_w/(ρ_s·1000))·(Qt / (A·dx))·dt
 * 
 * where Qt = sediment transport, ρ_s = sediment relative density (≈2.65)
 * 
 * @note Physical interpretation:
 * - Salinity: +1 psu → +0.78 kg/m³ (denser)
 * - Temperature: +1°C → -0.2 kg/m³ (lighter)
 * - Freshwater densest at 4°C (hence T-4 term)
 * 
 * @warning Linear EOS limitations:
 * - Valid for: 0-40 psu, 0-30°C, pressures < 100 bar
 * - For higher accuracy use UNESCO EOS-80 or TEOS-10
 * - Neglects pressure effect (OK for shallow water)
 * 
 * @see Gill (1982): Atmosphere-Ocean Dynamics, Appendix 3
 */
void CSimulation::calculate_density()
{
    // ⚡ OPTIMIZATION: Precalculate constants outside loop (checked once per timestep)
    // Parámetros betaS y betaT
    double betaS = bGetDoWaterSalinity() ? dGetBetaSalinityConstant() : 0.0;
    double betaT = bGetDoWaterTemperature() ? 
                   (bGetDoBetaCoefficient() && m_dBetaTemperatureConstant > 0.0 ? 
                    m_dBetaTemperatureConstant : dGetBetaTemperatureConstant()) : 0.0;
    
    // Cache constant factors
    const double rho_base = 1000.0;
    const double T_ref = 4.0;  // Reference temperature for density calculation
    
    for (int i = 0; i < m_nCrossSectionsNumber; i++)
    {
        const double S = m_vCrossSectionSalinity[i];
        const double T = m_vCrossSectionTemperature[i];
        const double rhow = rho_base * (1.0 + betaS * S - betaT * (T - T_ref));
        if (m_bDoSedimentTransport) {
            if (m_nPredictor == 1) {
                m_vPredictedCrossSectionDensity[i] =  rhow + (1 - rhow/1000.0/m_vCrossSectionRhos[i])*m_vCrossSectionQt[i]/(m_vCrossSectionArea[i]*m_vCrossSectionDX[i])*m_dTimestep;
            }
            else 
                m_vCrossSectionDensity[i] = rhow + (1 - rhow/1000.0/m_vCrossSectionRhos[i])*m_vCrossSectionQt[i]/(m_vCrossSectionArea[i]*m_vCrossSectionDX[i])*m_dTimestep;
            }
        else {
            if (m_nPredictor == 1) {
                m_vPredictedCrossSectionDensity[i] = rhow;
            }
            else {
                m_vCrossSectionDensity[i] = rhow;
            }
        }
    }
}

//===============================================================================================================================
//! Notifies the user that the simulation has ended, asks for keypress if necessary, and if compiled under GNU can send an email
//===============================================================================================================================
void CSimulation::bDoSimulationEnd(){
    presenter.EndingRun();

    const string strText = strGetErrorText(m_nStringError);
    CScreenPresenter::AnnounceEnding(strText);
};


/**
 * @brief Precompute and cache estuary geometry/hydraulics for O(1) runtime access
 * 
 * Optimization strategy:
 * - Extract data from CCrossSection objects into flat vectors
 * - Avoids repeated virtual function calls (e.g., estuary[i].dGetZ())
 * - Enables better CPU cache locality (contiguous memory)
 * - Precompute expensive terms: A·R_h^{2/3} for Manning friction
 * 
 * Cached data:
 * - Scalars: bed elevation, Manning's n, x-position, β coefficient
 * - Vectors: width(h), area(h), Rh(h), depth(h) - tabulated vs water level
 * - Bank locations: left/right riverbank positions
 * - Friction term: A·R_h^{2/3} (saves one pow() call per timestep)
 * 
 * @note Performance impact:
 * - ~10-15% speedup for large domains (N > 500)
 * - Most benefit in calculateHydraulicParameters() (called every timestep)
 * - Trade-off: ~50 MB extra memory for N=1000 domain
 * 
 * @warning Must be called after:
 * - estuary[] objects fully initialized
 * - All cross-section data loaded from CSV
 * 
 * @see calculateHydraulicParameters() for usage of cached data
 */
void CSimulation::precomputeEstuaryData() {
    
    // Pre-calcular datos escalares básicos
    m_vElevationSectionsCount.resize(m_nCrossSectionsNumber);
    m_vBedZ.resize(m_nCrossSectionsNumber);
    m_vManningN.resize(m_nCrossSectionsNumber);
    m_vPositionX.resize(m_nCrossSectionsNumber);
    m_vBeta.resize(m_nCrossSectionsNumber);
    
    // Initialize binary search cache (speeds up calculateHydraulicParameters ~30%)
    m_vLastInterpolationIndex.assign(m_nCrossSectionsNumber, 0);
    
    // ✅ MEJOR: Inicializar con dimensiones conocidas
    m_vWidth.assign(m_nCrossSectionsNumber, vector<double>());        // Vector vacío inicial       
    m_vLeftY.assign(m_nCrossSectionsNumber, vector<double>());        
    m_vRightY.assign(m_nCrossSectionsNumber, vector<double>());       
    m_vEstuaryAreas.assign(m_nCrossSectionsNumber, vector<double>());
    m_vEstuaryHydraulicRadius.assign(m_nCrossSectionsNumber, vector<double>());
    m_vEstuaryWaterDepths.assign(m_nCrossSectionsNumber, vector<double>());
    m_vPrecalculatedSecondTerm.assign(m_nCrossSectionsNumber, vector<double>());

    
    for (int i = 0; i < m_nCrossSectionsNumber; i++) {
        // Datos escalares básicos
        m_vBedZ[i] = estuary[i].dGetZ();
        m_vManningN[i] = estuary[i].dGetManningNumber();
        m_vPositionX[i] = estuary[i].dGetX();
        m_vElevationSectionsCount[i] = estuary[i].nGetElevationSectionsNumber();
        m_vBeta[i] = estuary[i].dGetBeta();
        
        // Vectores de datos hidráulicos
        m_vWidth[i] = estuary[i].vGetWidth();   
        m_vEstuaryAreas[i] = estuary[i].vGetArea();
        m_vEstuaryHydraulicRadius[i] = estuary[i].vGetHydraulicRadius();
        m_vEstuaryWaterDepths[i] = estuary[i].vGetWaterDepth();
        m_vLeftY[i] = estuary[i].vGetLeftRBLocation();
        m_vRightY[i] = estuary[i].vGetRightRBLocation();

        const auto& areas = m_vEstuaryAreas[i];
        const auto& hydraulicRadius = m_vEstuaryHydraulicRadius[i];
        
        m_vPrecalculatedSecondTerm[i].resize(areas.size());
        for (size_t j = 0; j < areas.size(); j++) {
            m_vPrecalculatedSecondTerm[i][j] = areas[j] * pow(hydraulicRadius[j], 2.0/3.0);
        }
    }
    }

/**
 * @brief Generate descriptive NetCDF output filename from simulation parameters
 * 
 * Format: barrier_sim_YYYYMMDD_HHMM_CSxxx_Txxx_dtxxx_BCxy_[OPTIONS]_CFLxx.nc
 * 
 * Components:
 * - Date/Time: Current wall-clock time when simulation starts
 * - CS: Number of cross-sections (zero-padded to 3 digits)
 * - T: Simulation duration (minutes < 60, hours < 24, else days)
 * - dt: Timestep in seconds
 * - BC: Boundary conditions (xy = upstream/downstream type)
 * - OPTIONS: _SED (sediment), _SAL (salinity), _TVD (limiters), _DRY (dry bed)
 * - CFL: Courant number × 100 (e.g., CFL15 = 0.15)
 * 
 * Example:
 *   barrier_sim_20240115_1430_CS100_T7d_dt3600_BC12_SAL_TVD_DRY_CFL15.nc
 * 
 * @return Formatted filename string
 * 
 * @note Purpose:
 * - Self-documenting filenames for parameter studies
 * - Easy sorting by date
 * - Unique names prevent accidental overwrites
 * - Quickly identify simulation configuration without opening file
 * 
 * @see CDataWriter for NetCDF attribute metadata (more detailed)
 */
std::string CSimulation::generateOutputFileName() const {
    std::ostringstream filename;
    
    // ✅ FECHA Y HORA actual
    auto now = std::time(nullptr);
    auto tm = *std::localtime(&now);
    
    // ✅ NOMBRE BASE del proyecto
    filename << "barrier_sim";
    
    // ✅ FECHA: YYYYMMDD_HHMM
    filename << "_" << std::put_time(&tm, "%Y%m%d_%H%M");
    
    // ✅ NÚMERO DE SECCIONES
    filename << "_CS" << std::setfill('0') << std::setw(3) << m_nCrossSectionsNumber;
    
    // ✅ DURACIÓN DE SIMULACIÓN (en horas)
    double hours = m_dSimDuration / 3600.0;
    if (hours < 1.0) {
        filename << "_T" << std::setfill('0') << std::setw(2) << static_cast<int>(m_dSimDuration / 60.0) << "min";
    } else if (hours < 24.0) {
        filename << "_T" << std::setfill('0') << std::setw(2) << static_cast<int>(hours) << "h";
    } else {
        filename << "_T" << static_cast<int>(hours / 24.0) << "d";
    }
    
    // ✅ TIMESTEP (en segundos)
    filename << "_dt" << std::setfill('0') << std::setw(3) << static_cast<int>(m_dSimTimestep);
    
    // ✅ CONDICIONES DE FRONTERA
    filename << "_BC" << m_nUpwardEstuarineCondition << m_nDownwardEstuarineCondition;
    
    // ✅ OPCIONES ESPECIALES
    if (m_bDoSedimentTransport) filename << "_SED";
    if (m_bDoWaterSalinity) filename << "_SAL";
    if (m_bDoMcCormackLimiterFlux) filename << "_TVD";
    if (m_bDoDryBed) filename << "_DRY";
    
    // ✅ COURANT NUMBER
    filename << "_CFL" << std::setfill('0') << std::setw(2) << static_cast<int>(m_dCourantNumber * 100);
    
    // ✅ EXTENSIÓN
    filename << ".nc";
    
    return filename.str();
}




// === INLINE HELPER FUNCTIONS ===
// calc_Lv(), calc_rho_air(), and calc_albedo_briegleb() are defined as inline functions
// in simulation.h for maximum performance (~10x faster than function calls for small functions).

