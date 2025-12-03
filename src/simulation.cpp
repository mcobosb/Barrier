#include <iostream>

using std::cerr;
using std::cout;
using std::endl;
using std::ios;

#include <iomanip>
using std::setprecision;
using std::setw;

#include <string>
using std::to_string;

#include <cmath>
using std::sqrt;
using std::fabs;
using std::pow;

#include "simulation.h"
#include "screen_presenter.h"
#include "hydrograph.h"
#include "cross_section.h"
#include "main.h"
#include "utils.h"
#include "error_handling.h"



//===============================================================================================================================
//! The CSimulation constructor
//===============================================================================================================================
CSimulation::CSimulation() {
    m_dSimDuration = 0.0;
    m_dSimTimestep = 0.0;
    m_dTimeFactor = 0.0;
    m_dTimestep = 0.0;
    m_dLambda = 0.0;
    m_dCurrentTime = 0.0;

    m_nPredictor = -1;
    m_nTimeLogId = 0;
    m_nStringError = 0;
    m_nCrossSectionsNumber = 0;        
    m_nHydrographsNumber = 0;            
    m_nTimeId = 0;                      

    m_nInitialEstuarineCondition = 0;
    m_nLogFileDetail = 0;
    m_nEquationSedimentTransport = 0;
    m_nUpwardSalinityCondition = 0;
    m_nDownwardSalinityCondition = 0;
    m_nUpwardEstuarineCondition = 0;     
    m_nDownwardEstuarineCondition = 0;  

    m_bSaveTime = false;
    m_bDoSedimentTransport = false;
    m_bHydroFile = false;
    m_bReturnError = false;
    m_bDoWaterSalinity = false;          
    m_bDoMcCormackLimiterFlux = false;  
    m_bDoSurfaceGradientMethod = false;  
    m_bDoSourceTermBalance = false;     
    m_bDoBetaCoefficient = false;        
    m_bDoDryBed = false;                
    m_bDoMurilloCondition = false;       
    m_bDoWaterDensity = false;          

    m_dBetaSalinityConstant = 0.0;      
    m_dLongitudinalDispersion = 0.0;     
    m_dUpwardBoundaryValue = 0.0;        
    m_dNextUpwardBoundaryValue = 0.0;    
    m_dDownwardBoundaryValue = 0.0;      
    m_dNextDownwardBoundaryValue = 0.0;  
    m_dCourantNumber = 0.0;              
    m_nEquationMcCormackLimiterFlux = 0; 
    m_nPsiFormula = 0;                   
    m_dDeltaValue = 0.0;                 

    // ✅ FECHA DE INICIO
    m_nSimStartSec = 0;
    m_nSimStartMin = 0;
    m_nSimStartHour = 0;
    m_nSimStartDay = 1;
    m_nSimStartMonth = 1;
    m_nSimStartYear = 2024;

    m_vOutputVariables.clear();
    estuary.clear();
    hydrographs.clear();

    m_tSysStartTime = 0;
    m_tSysStartLoopTime = 0;
}

//===============================================================================================================================
//! The CSimulation destructor
//===============================================================================================================================
CSimulation::~CSimulation () {
    // Close output files if open
    if (LogStream && LogStream.is_open())
    {
        LogStream.flush();
        LogStream.close();
    }
};


//! Method for getting the simulation duration
double CSimulation::dGetSimulationDuration() const {
    return m_dSimDuration;
}

//! Method for setting the simulation duration
void CSimulation::dSetSimulationDuration(double simDuration) {
    m_dSimDuration = simDuration;
};

//! Method for getting the simulation timestep
double CSimulation::dGetSimulationTimestep() const {
    return m_dSimTimestep;
}

//! Method for setting the simulation timestep
void CSimulation::dSetSimulationTimestep(double simTimestep) {
    m_dSimTimestep = simTimestep;
};

//! Method for getting the initial estuarine condition
int CSimulation::nGetInitialEstuarineCondition() const {
    return m_nInitialEstuarineCondition;
}

//! Method for setting the initial estuarine condition
void CSimulation::nSetInitialEstuarineCondition(int initialCondition) {
    m_nInitialEstuarineCondition = initialCondition;
};

//! Method for getting the upward estuarine condition
int CSimulation::nGetUpwardEstuarineCondition() const {
    return m_nUpwardEstuarineCondition;
}

//! Method for setting the upward estuarine condition
void CSimulation::nSetUpwardEstuarineCondition(int upwardEstuarineCondition) {
    m_nUpwardEstuarineCondition = upwardEstuarineCondition;
};

//! Method for getting the downward estuarine condition
int CSimulation::nGetDownwardEstuarineCondition() const {
    return m_nDownwardEstuarineCondition;
}

//! Method for setting the downward estuarine condition
void CSimulation::nSetDownwardEstuarineCondition(int downwardEstuarineCondition) {
    m_nDownwardEstuarineCondition = downwardEstuarineCondition;
};

//! Method for getting the courant number
double CSimulation::dGetCourantNumber() const {
    return m_dCourantNumber;
}

//! Method for setting the courant number
void CSimulation::dSetCourantNumber(double courantNumber) {
    m_dCourantNumber = courantNumber;
}

//! Method for getting if McComarck limiter flux is applied
bool CSimulation::bGetDoMcComarckLimiterFlux() const {
    return m_bDoMcCormackLimiterFlux;
}
//! Method for setting if McComarck limiter flux is applied
void CSimulation::bSetDoMcComarckLimiterFlux(bool doMcComarckLimiterFlux) {
    m_bDoMcCormackLimiterFlux = doMcComarckLimiterFlux;
}

//! Method for getting equation limiter flux
int CSimulation::nGetEquationLimiterFlux() const {
    return m_nEquationMcCormackLimiterFlux;
}
//! Method for setting equation limiter flux
void CSimulation::nSetEquationLimiterFlux(int equationLimiterFlux) {
    m_nEquationMcCormackLimiterFlux = equationLimiterFlux;
}

//! Method for getting Psi formula
int CSimulation::nGetPsiFormula() const {
    return m_nPsiFormula;
}
//! Method for setting Psi formula
void CSimulation::nSetPsiFormula(int psiFormula) {
    m_nPsiFormula = psiFormula;
}

//! Method for getting Delta Value
double CSimulation::dGetDeltaValue() const {
    return m_dDeltaValue;
}
//! Method for setting Delta Value
void CSimulation::dSetDeltaValue(double deltaValue) {
    m_dDeltaValue = deltaValue;
}

//! Method for getting if surface gradient method is applied
bool CSimulation::bGetDoSurfaceGradientMethod() const {
    return m_bDoSurfaceGradientMethod;
}
//! Method for setting if surface gradient method is applied
void CSimulation::bSetDoSurfaceGradientMethod(bool doSurfaceGradientMethod) {
    m_bDoSurfaceGradientMethod = doSurfaceGradientMethod;
}

//! Method for getting if source term balance is applied
bool CSimulation::bGetDoSurfaceTermBalance() const {
    return m_bDoSourceTermBalance;
}
//! Method for setting if source term balance is applied
void CSimulation::bSetDoSurfaceTermBalance(bool doSourceTermBalance) {
    m_bDoSourceTermBalance = doSourceTermBalance;
}

//! Method for getting if beta coefficient is applied
bool CSimulation::bGetDoBetaCoefficient() const {
    return m_bDoBetaCoefficient;
}
//! Method for setting if beta coefficient is applied
void CSimulation::bSetDoBetaCoefficient(bool doBetaCoefficient) {
    m_bDoBetaCoefficient = doBetaCoefficient;
}

//! Method for getting if dry bed is applied
bool CSimulation::bGetDoDryBed() const {
    return m_bDoDryBed;
}
//! Method for setting if dry bed is applied
void CSimulation::bSetDoDryBed(bool doDryBed) {
    m_bDoDryBed = doDryBed;
}

//! Method for getting if Murillo condition is applied
bool CSimulation::bGetDoMurilloCondition() const {
    return m_bDoMurilloCondition;
}

//! Method for setting if Murillo condition is applied
void CSimulation::bSetDoMurilloCondition(bool doMurilloCondition) {
    m_bDoMurilloCondition = doMurilloCondition;
}

//! Method for getting the computation of water density
bool CSimulation::bGetDoWaterSalinity() const {
    return m_bDoWaterSalinity;
}

//! Method for setting the computation of water density
void CSimulation::bSetDoWaterSalinity(bool computeWaterSalinity) {
    m_bDoWaterSalinity = computeWaterSalinity;
}

//! Method for getting the upward salinity condition
int CSimulation::nGetUpwardSalinityCondition() const {
    return m_nUpwardSalinityCondition;
}

//! Method for setting the upward salinity condition
void CSimulation::nSetUpwardSalinityCondition(const int nUpwardCondition) {
    m_nUpwardSalinityCondition = nUpwardCondition;
}

//! Method for getting the downward salinity condition
int CSimulation::nGetDownwardSalinityCondition() const {
    return m_nDownwardSalinityCondition;
}

//! Method for setting the downward salinity condition
void CSimulation::nSetDownwardSalinityCondition(const int nDownwardCondition) {
    m_nDownwardSalinityCondition = nDownwardCondition;
}

//! Method for getting the along channel constant elevation
double CSimulation::dGetBetaSalinityConstant() const {
    return m_dBetaSalinityConstant;
}

//! Method for setting the along channel constant elevation
void CSimulation::dSetBetaSalinityConstant(double betaSalinity) {
    m_dBetaSalinityConstant = betaSalinity;
}

//! Method for getting the along channel constant elevation
double CSimulation::dGetLongitudinalDispersionConstant() const {
    return m_dLongitudinalDispersion;
}

//! Method for setting the along channel constant elevation
void CSimulation::dSetLongitudinalDispersionConstant(double longitudinalDispersion) {
    m_dLongitudinalDispersion = longitudinalDispersion;
};


//! Method for getting the computation of sediment transport
bool CSimulation::bGetDoSedimentTransport() const {
    return m_bDoSedimentTransport;
}

//! Method for setting the computation of sediment transport
void CSimulation::bSetDoSedimentTransport(bool computeSedimentTransport) {
    m_bDoSedimentTransport = computeSedimentTransport;
}

//! Method for getting equation for sediment transport
int CSimulation::nGetEquationSedimentTransport() const {
    return m_nEquationSedimentTransport;
}

//! Method for setting equation limiter flux
void CSimulation::nSetEquationSedimentTransport(int equationSedimentTransport) {
    m_nEquationSedimentTransport = equationSedimentTransport;
}

//! Method for getting the computation of water density
bool CSimulation::bGetDoWaterDensity() const {
    return m_bDoWaterDensity;
}

//! Method for setting the computation of water density
void CSimulation::bSetDoWaterDensity(const bool doWaterDensity) {
    m_bDoWaterDensity = doWaterDensity;
}

//! Method for getting the save all timesteps flag
bool CSimulation::bGetSaveAllTimesteps() const {
    return m_bSaveAllTimesteps;
}

//! Method for setting the save all timesteps flag
void CSimulation::bSetSaveAllTimesteps(const bool saveAllTimesteps) {
    m_bSaveAllTimesteps = saveAllTimesteps;
}

//===============================================================================================================================
//! Appends output variables
//===============================================================================================================================
void CSimulation::strAddOutputVariable(const string& strItem){
    m_vOutputVariables.push_back(strItem);
}

//===============================================================================================================================
//! Appends a CHydrograph objet to the CSimulation
//===============================================================================================================================
void CSimulation::AddHydrograph(){
    // Create the
    CHydrograph hydrograph;
    hydrographs.push_back(hydrograph);
}


//! Method for getting the number of hydrographs
int CSimulation::nGetHydrographsNumber() const{
    return m_nHydrographsNumber;
}
//! Method for setting the number of hydrographs
void CSimulation::nSetHydrographsNumber(int nValue) {
    m_nHydrographsNumber = nValue;
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

    //! Read the configuration file
    reader.bReadConfigurationFile(this);
    if (m_bReturnError)
    {
        return;
    }

    CDataReader::bOpenLogFile(this);
    reader.bReadAlongChannelDataFile(this);
    reader.bReadCrossSectionGeometryFile(this);
    if (m_nUpwardEstuarineCondition > 1) {
        CDataReader::bReadUpwardBoundaryConditionFile(this);
    }
    if (m_nDownwardEstuarineCondition > 1) {
        CDataReader::bReadDownwardBoundaryConditionFile(this);
    }
    if (m_bDoSedimentTransport)
    {
        reader.bReadAlongChannelSedimentsFile(this);
    }

    if (m_bHydroFile)
    {
        reader.bReadHydrographsFile(this);
    }


    initializeVectors();
    precomputeEstuaryData();
    calculateBedSlope();
    calculateAlongEstuaryInitialConditions();


    // ✅ CON esto:
    std::string m_strOutFileName = generateOutputFileName();
    m_strOutFile += m_strOutFileName;
    
    std::cout << "📁 Output file: " << m_strOutFile << std::endl;
    //! Create the NetCDF file with dimensions, coordinates and variables
    writer.nDefineNetCDFFile(this);

    //! Check if any part of the estuary is dry
    if (bGetDoDryBed())
        dryArea();

    m_nTimeId = 0;
    //! Save initial timestep
    writer.nSetOutputData(this);
    m_nTimeId++;

    int m_nStep = 1;
    cout << "    - Running" << endl;
    m_tSysStartLoopTime = time(nullptr);

    //! Main time-stepping loop: advance simulation until final time
    while (m_dCurrentTime <= m_dSimDuration)
    {
        calculateBoundaryConditions();

        m_nPredictor = 1;
        //! Calculate the hydraulic parameters for prediction
        calculateHydraulicParameters();

        //! Check if Murillo condition is applied
        if (bGetDoMurilloCondition())
            doMurilloCondition();

        ///! Calculate the timestep using the Courant Number
        calculateTimestep();

        //==============================================================================================================
        //! STEP 01: Compute the predictor
        //==============================================================================================================

        // Compute the sediment transport
        if (bGetDoSedimentTransport()) calculate_sediment_transport();

        // Calculate density due to salinity and sediment concentration
        if (bGetDoWaterSalinity()) calculate_density();

        //! Compute source and sink terms
        calculate_GS_A_terms();
        calculateFlowTerms();
        calculateSourceTerms();

        //! Change the values of the source terms if any part of the estuary is dry
        if (bGetDoDryBed())
            dryTerms();

        calculatePredictor();

        //! Boundary conditions
        updatePredictorBoundaries();

        //! Check if any part of the estuary is dry (predictor)
        if (bGetDoDryBed())
            dryArea();

        //! Save predicted variables
        // if (m_bSaveTime) {
        //     int nRtn = writer.vSetOutputData()
        //
        // }

        //! Calculate the hydraulic parameters for prediction
        calculateHydraulicParameters();

        //! Check if Murillo condition is applied
        if (bGetDoMurilloCondition())
            doMurilloCondition();

        //==============================================================================================================
        //! STEP 02: Compute the corrector
        //==============================================================================================================
        m_nPredictor = 2;
        calculate_GS_A_terms();
        calculateFlowTerms();
        calculateSourceTerms();

        //! Change the values of the source terms if any part of the estuary is dry
        if (bGetDoDryBed())
            dryTerms();

        calculateCorrector();

        //! Boundary conditions
        updateCorrectorBoundaries();
        //! Check if any part of the estuary is dry (predictor)
        if (bGetDoDryBed())
            dryArea();

        //==============================================================================================================
        //! STEP 03: Merge predictor-corrector and apply TVD limiters
        //==============================================================================================================
        mergePredictorCorrector();

        m_nPredictor = 0;
        //! Check if any part of the estuary is dry after merge
        if (bGetDoDryBed())
            dryArea();

        //! Compute salinity gradient for next timestep
        if (bGetDoWaterSalinity())
        {
            calculate_salinity_gradient();
        }

        //! Compute salinity transport
        if (bGetDoWaterSalinity())
        {
            calculate_salinity();
        }

        //! Calculate the hydraulic parameters for saving after the correction
        // calculateHydraulicParameters();

        AnnounceProgress();

        if (m_bSaveTime || (m_nLogFileDetail == 2)) {
            // 🎨 Suavizar solución solo para output (no afecta cálculo)
            smoothSolution();
            writer.nSetOutputData(this);
        }

        //! Increase counter
        m_dCurrentTime = m_dCurrentTime + m_dTimestep;
        if (m_bSaveTime) {
            m_nTimeId++;
        }
        m_nStep++;
    }

    // Avoid overwriting the last line
    cout << endl;
    writer.nCloseNetCDFFile();

    // Cerrar log file si está abierto
    if (LogStream.is_open()) {
        LogStream.close();
    }
}


//======================================================================================================================
//! Initialize vector variables
//======================================================================================================================
void CSimulation::initializeVectors()
{

    const double nCrossSectionsNumber = m_nCrossSectionsNumber;
    const vector<double> vZeros(static_cast<size_t>(nCrossSectionsNumber), 0.0);

    // m_vCrossSectionArea =
    m_vPredictedCrossSectionArea =
    m_vCorrectedCrossSectionArea =
    m_vPredictedCrossSectionQ =
    m_vCorrectedCrossSectionQ =
    m_vCrossSectionHydraulicRadius =
    m_vCrossSectionDhDx =
    m_vCrossSectionWidth =
    m_vCrossSectionBeta =
    m_vCrossSectionI1 =
    m_vPredictedCrossSectionI1 =
    m_vCrossSectionI2 =
    m_vCrossSectionU =
    m_vCrossSectionC =
    m_vCrossSectionSalinity =
    m_vCrossSectionQb =
    m_vCrossSectionQs =
    m_vCrossSectionQt =
    m_vCrossSectionLeftRBLocation =
    m_vCrossSectionRightRBLocation =
    m_vCrossSectionBedSlope =
    m_vCrossSectionBedSlopePredictor =
    m_vCrossSectionBedSlopeCorrector =
    m_vCrossSectionFrictionSlope =
    m_vCrossSectionDX =
    m_vCrossSectionManningNumber =
    m_vCrossSection_gAS0 =
    m_vCrossSection_gASf =
    m_vCrossSectionF0 =
    m_vCrossSectionF1 =
    m_vCrossSectionGv0 =
    m_vCrossSectionGv1 =
    m_vCrossSectionMurilloFactor =
    m_vLateralSourcesAtT =
    m_vCrossSectionDiamX =
    m_vCrossSectionSalinityASt = vZeros;

    const vector<double> vRhos(static_cast<size_t>(nCrossSectionsNumber), FRESH_WATER_DENSITY);
    m_vCrossSectionRho = vRhos;

    const vector<double> vOnes(static_cast<size_t>(nCrossSectionsNumber+1), 1.0);
    m_vCrossSectionD1Factor =
    m_vCrossSectionD2Factor = vOnes;

    const int nTimestepsNumber = static_cast<int> (m_dSimDuration/m_dSimTimestep) + 1;
    
    // Reservar espacio antes del bucle paralelo
    m_vOutputTimesIds.resize(nTimestepsNumber);
    m_vOutputTimes.resize(nTimestepsNumber);

    for (int i = 0; i < nTimestepsNumber; i++)
    {
        m_vOutputTimesIds[i] = i;
        m_vOutputTimes[i] = static_cast<double>(i) * m_dSimTimestep;  // Usar [] en lugar de push_back
    }

    //! ✅ Ya no necesitamos m_vCrossSectionBedSlopeDirection (eliminado)
}

//======================================================================================================================
//! Calculate the initial bed slope
//======================================================================================================================
void CSimulation::calculateBedSlope() {
    double dX1;
    double dX2;
    
    //! Calculate bed slope S0 using central differences (interior nodes) and forward/backward (boundaries)
    for (int i = 0; i < m_nCrossSectionsNumber ; i++) {
        if (i == 0) {
            // Compute dx as forward difference
            dX1 = estuary[1].dGetX() - estuary[0].dGetX();
            m_vCrossSectionBedSlope[i] = (estuary[0].dGetZ() - estuary[1].dGetZ())/dX1;
            //! Save dX into a vector
            m_vCrossSectionDX[i] = dX1;
        }
        else if (i == m_nCrossSectionsNumber - 1){
            // Compute dx as backward difference
            dX2 = estuary[i].dGetX() - estuary[i-1].dGetX();
            m_vCrossSectionBedSlope[i] = (estuary[i-1].dGetZ() - estuary[i].dGetZ())/dX2;
            //! Save dX into a vector (last node backward)
            m_vCrossSectionDX[i] = dX2;
        }
        else {
            dX1 = estuary[i+1].dGetX() - estuary[i].dGetX();
            dX2 = estuary[i].dGetX() - estuary[i-1].dGetX();
            m_vCrossSectionBedSlope[i] = (estuary[i-1].dGetZ() - estuary[i+1].dGetZ())/(dX1 + dX2);

            //! Save dX into a vector
            m_vCrossSectionDX[i] = dX1;
        }
        
        //! ✅ Ya no necesitamos el factor de dirección confuso
        //! La ecuación de Manning usa directamente sign(S0) * sqrt(|S0|)
    }
    
    //! ✅ CRÍTICO: Calcular pendientes unilaterales para balance de términos fuente (como Fortran)
    //! zmedp(ii) = (z(ii) - z(ii+1))/dx   para predictor (diferencia hacia adelante)
    //! zmedc(ii) = (z(ii-1) - z(ii))/dx   para corrector (diferencia hacia atrás)
    for (int i = 0; i < m_nCrossSectionsNumber; i++) {
        if (i != m_nCrossSectionsNumber - 1) {
            // Predictor: diferencia hacia adelante (z[i] - z[i+1])/dx
            double dx_forward = estuary[i+1].dGetX() - estuary[i].dGetX();
            m_vCrossSectionBedSlopePredictor[i] = (estuary[i].dGetZ() - estuary[i+1].dGetZ()) / dx_forward;
        } else {
            // Último nodo: usar diferencia hacia atrás
            double dx_backward = estuary[i].dGetX() - estuary[i-1].dGetX();
            m_vCrossSectionBedSlopePredictor[i] = (estuary[i-1].dGetZ() - estuary[i].dGetZ()) / dx_backward;
        }
        
        if (i != 0) {
            // Corrector: diferencia hacia atrás (z[i-1] - z[i])/dx
            double dx_backward = estuary[i].dGetX() - estuary[i-1].dGetX();
            m_vCrossSectionBedSlopeCorrector[i] = (estuary[i-1].dGetZ() - estuary[i].dGetZ()) / dx_backward;
        } else {
            // Primer nodo: usar diferencia hacia adelante
            double dx_forward = estuary[i+1].dGetX() - estuary[i].dGetX();
            m_vCrossSectionBedSlopeCorrector[i] = (estuary[i].dGetZ() - estuary[i+1].dGetZ()) / dx_forward;
        }
    }
}



//======================================================================================================================
//! Calculate the initial conditions along the estuary
//======================================================================================================================
void CSimulation::calculateAlongEstuaryInitialConditions() {

    // Reservar espacio para todos los vectores antes de los bucles paralelos
    m_vCrossSectionWaterDepth.resize(m_nCrossSectionsNumber);
    m_vCrossSectionWaterElevation.resize(m_nCrossSectionsNumber);

    // ⚡ Pre-allocar vectores de trabajo TVD una sola vez
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

    // ⚡ Pre-allocar vectores de salinidad una sola vez
    m_vSalinity_KAS_forward.resize(m_nCrossSectionsNumber);
    m_vSalinity_KAS_backward.resize(m_nCrossSectionsNumber);
    m_vSalinity_AUS_diff.resize(m_nCrossSectionsNumber);

    // ⚡ PRECALCULAR CONSTANTES (valores que no cambian durante la simulación)
    m_vManningNumberSquared.resize(m_nCrossSectionsNumber);
    m_vInvDX.resize(m_nCrossSectionsNumber);
    m_vDxSum.resize(m_nCrossSectionsNumber);
    m_vInvDxSum.resize(m_nCrossSectionsNumber);
    m_vGtimesDX.resize(m_nCrossSectionsNumber);
    m_vCrossSectionDI1Dx.resize(m_nCrossSectionsNumber, 0.0);  // Gradiente de I1
    
    for (int i = 0; i < m_nCrossSectionsNumber; i++) {
        // Manning² se usa ~4000 veces por día de simulación
        m_vManningNumberSquared[i] = pow(m_vCrossSectionManningNumber[i], 2.0);
        
        // 1/ΔX se usa en múltiples gradientes
        if (m_vCrossSectionDX[i] > 1e-10) {
            m_vInvDX[i] = 1.0 / m_vCrossSectionDX[i];
        } else {
            m_vInvDX[i] = 0.0;
        }
        
        // g*ΔX término constante
        m_vGtimesDX[i] = G * m_vCrossSectionDX[i];
        
        // Calcular área mínima para cada cross-section (10% del área mínima de geometría)
        // Esto evita problemas cuando DRY_AREA es menor que la geometría mínima del canal
        DRY_AREA = 0.1 * m_vEstuaryAreas[i][0];
    }
    
    // Precalcular sumas de ΔX para diferencias centradas
    for (int i = 1; i < m_nCrossSectionsNumber - 1; i++) {
        double dx1 = m_vPositionX[i+1] - m_vPositionX[i];
        double dx2 = m_vPositionX[i] - m_vPositionX[i-1];
        m_vDxSum[i] = dx1 + dx2;
        m_vInvDxSum[i] = 1.0 / m_vDxSum[i];
    }

    if (m_nInitialEstuarineCondition == 1) {
        //! Along estuary water flow given
        for (int i = 0; i < m_nCrossSectionsNumber; i++) {
            double dManningFactor = 0.0;
            //! As initial condition, it is assumed the independency of A with +/- value of S0
            if (m_vCrossSectionBedSlope[i] == 0.0) {
                //! Only for initial conditions S0 = 1e-3
                dManningFactor = m_vCrossSectionQ[i]*estuary[i].dGetManningNumber()/(sqrt(1e-3));
            }
            else {
                dManningFactor = m_vCrossSectionQ[i]*estuary[i].dGetManningNumber()/(sqrt(fabs(m_vCrossSectionBedSlope[i])));
            }

            vector<double> vCrossSectionAreaTmp = estuary[i].vGetArea();
            vector<double> vCrossSectionHydraulicRadiusTmp = estuary[i].vGetHydraulicRadius();
            vector<double> vCrossSectionElevationTmp = estuary[i].vGetWaterDepth();
            
            // Reservar espacio para dSecondTerm antes del bucle interno
            vector<double> dSecondTerm;
            dSecondTerm.resize(vCrossSectionAreaTmp.size());

            //! Second term to obtain the area from slope equation in open channels
            for (size_t j = 0; j < vCrossSectionAreaTmp.size(); j++) {
                dSecondTerm[j] = vCrossSectionAreaTmp[j]*pow(vCrossSectionHydraulicRadiusTmp[j], 2.0/3.0);
            }
            
            m_vCrossSectionArea[i] = linearInterpolation1d(dManningFactor, dSecondTerm, vCrossSectionAreaTmp);
            m_vCrossSectionWaterDepth[i] = linearInterpolation1d(m_vCrossSectionArea[i], vCrossSectionAreaTmp, vCrossSectionElevationTmp);
            m_vCrossSectionWidth[i] = m_vCrossSectionArea[i] / m_vCrossSectionWaterDepth[i];
            m_vCrossSectionWaterElevation[i] = m_vCrossSectionWaterDepth[i] + estuary[i].dGetZ();
        }
    }
    else if (m_nInitialEstuarineCondition == 2) {
        //! Along estuary elevation given
        for (int i = 0; i < m_nCrossSectionsNumber; i++) {
            vector<double> vCrossSectionAreaTmp = estuary[i].vGetArea();
            vector<double> vCrossSectionElevationTmp = estuary[i].vGetWaterDepth();
            vector<double> vCrossSectionHydraulicRadiusTmp = estuary[i].vGetHydraulicRadius();

            //! m_vCrossSectionElevation is the elevation from the water depth m_dZ of every cross-section
            m_vCrossSectionArea[i] = linearInterpolation1d(m_vCrossSectionWaterElevation[i] - estuary[i].dGetZ(),vCrossSectionElevationTmp, vCrossSectionAreaTmp);
            m_vCrossSectionWidth[i] = m_vCrossSectionArea[i] / m_vCrossSectionWaterDepth[i];
            m_vCrossSectionHydraulicRadius[i] = linearInterpolation1d(m_vCrossSectionArea[i], vCrossSectionAreaTmp, vCrossSectionHydraulicRadiusTmp);

            //! ✅ Usar pendiente de la línea de energía (más robusta que pendiente del lecho)
            //! Q = (A * R^(2/3) / n) * sqrt(|S0|) * sign(S0)
            double sign_S0 = (m_vCrossSectionBedSlope[i] >= 0) ? 1.0 : -1.0;
            m_vCrossSectionQ[i] = m_vCrossSectionArea[i] * 
                                 pow(m_vCrossSectionHydraulicRadius[i], 2.0/3.0) * 
                                 sqrt(fabs(m_vCrossSectionBedSlope[i]) + 1e-10) * 
                                 sign_S0 / estuary[i].dGetManningNumber();
        }
    }
    else {
        //! Water level in calm - REPOSO HIDROSTÁTICO
        for (int i = 0; i < m_nCrossSectionsNumber; i++) {
            vector<double> vCrossSectionAreaTmp = estuary[i].vGetArea();
            vector<double> vCrossSectionElevationTmp = estuary[i].vGetWaterDepth();
            vector<double> vCrossSectionHydraulicRadiusTmp = estuary[i].vGetHydraulicRadius();

            // Calado desde el nivel de referencia (z=0)
            m_vCrossSectionWaterDepth[i] = -estuary[i].dGetZ();
            m_vCrossSectionWaterElevation[i] = 0.0;
            m_vCrossSectionArea[i] = linearInterpolation1d(m_vCrossSectionWaterDepth[i],vCrossSectionElevationTmp, vCrossSectionAreaTmp);
            m_vCrossSectionWidth[i] = m_vCrossSectionArea[i] / m_vCrossSectionWaterDepth[i];
            m_vCrossSectionHydraulicRadius[i] = linearInterpolation1d(m_vCrossSectionArea[i], vCrossSectionAreaTmp, vCrossSectionHydraulicRadiusTmp);

            //! ✅ CORRECTO: En reposo (calm), Q = 0 en todo el dominio
            m_vCrossSectionQ[i] = 0.0;
        }
    }
 
    for (int i = 0; i < m_nCrossSectionsNumber; i++) {
        //! Insert the Manning number onto the simulation object
        m_vCrossSectionManningNumber[i] = estuary[i].dGetManningNumber();
    }

    if (m_bDoWaterDensity)
    {
        //! Compute along channel sediment parameter    
        for (int i = 0; i < m_nCrossSectionsNumber; i++)
        {
            m_vCrossSectionDiamX[i] = m_vCrossSectionD50[i] * pow((m_vCrossSectionRhos[i] - 1.0) * G / (NU*NU), 1.0/3.0);
        }
    }
}


//======================================================================================================================
//! Interpolate function
//======================================================================================================================
double CSimulation::linearInterpolation1d(const double dValue, const vector<double> &vX, const vector<double> &vY) {
    // Búsqueda binaria
    auto it = std::lower_bound(vX.begin(), vX.end(), dValue);
    
    if (it == vX.begin()) {
        // Extrapolación hacia abajo
        double slope = (vY[1] - vY[0]) / (vX[1] - vX[0]);
        return vY[0] + slope * (dValue - vX[0]);
    } else if (it == vX.end()) {
        // Extrapolación hacia arriba
        size_t n = vX.size();
        double slope = (vY[n-1] - vY[n-2]) / (vX[n-1] - vX[n-2]);
        return vY[n-1] + slope * (dValue - vX[n-1]);
    } else {
        // Interpolación normal
        size_t i = std::distance(vX.begin(), it) - 1;
        double slope = (vY[i+1] - vY[i]) / (vX[i+1] - vX[i]);
        return vY[i] + slope * (dValue - vX[i]);
    }
}
//======================================================================================================================
//! Obtain the hydraulic parameters from the cross-sectional data
//======================================================================================================================
void CSimulation::calculateHydraulicParameters() {
    //! Number of estuarine cross-sections
    const int nCrossSections = m_nCrossSectionsNumber;
    
    // ✅ Referencia condicional (sin copia)
    const auto& dArea = (m_nPredictor == 1) ? m_vCrossSectionArea : m_vPredictedCrossSectionArea;
    
    // ⚡ OPTIMIZACIÓN: Loop optimizado con interpolación lineal rápida
    for (int i = 0; i < nCrossSections; i++) {
        const int nElevationSectionsNumber = m_vElevationSectionsCount[i];
        const double currentArea = dArea[i];
        
        // Caso 1: Área menor que mínimo
        if (m_vEstuaryAreas[i][0] > currentArea) {
            m_vCrossSectionHydraulicRadius[i] = m_vEstuaryHydraulicRadius[i][0];
            m_vCrossSectionWidth[i] = m_vWidth[i][0];
            m_vCrossSectionWaterDepth[i] = m_vEstuaryWaterDepths[i][0];
            m_vCrossSectionBeta[i] = m_vBeta[i][0];
            m_vCrossSectionLeftRBLocation[i] = m_vLeftY[i][0];
            m_vCrossSectionRightRBLocation[i] = m_vRightY[i][0];
            
            // Interpolate I1 for predictor/normal states
            if (m_nPredictor == 1) {
                m_vCrossSectionI1[i] = m_vEstuaryI1[i][0];
            } else {
                m_vPredictedCrossSectionI1[i] = m_vEstuaryI1[i][0];
            }
        }
        // Caso 2: Área mayor que máximo
        else if (m_vEstuaryAreas[i][nElevationSectionsNumber-1] < currentArea) {
            const int lastNode = nElevationSectionsNumber - 1;
            m_vCrossSectionHydraulicRadius[i] = m_vEstuaryHydraulicRadius[i][lastNode];
            m_vCrossSectionWidth[i] = m_vWidth[i][lastNode];
            m_vCrossSectionWaterDepth[i] = m_vEstuaryWaterDepths[i][lastNode];
            m_vCrossSectionBeta[i] = m_vBeta[i][lastNode];
            m_vCrossSectionLeftRBLocation[i] = m_vLeftY[i][lastNode];
            m_vCrossSectionRightRBLocation[i] = m_vRightY[i][lastNode];
            
            // Interpolate I1 for predictor/normal states
            if (m_nPredictor == 1) {
                m_vCrossSectionI1[i] = m_vEstuaryI1[i][lastNode];
            } else {
                m_vPredictedCrossSectionI1[i] = m_vEstuaryI1[i][lastNode];
            }
        }
        // Caso 3: Interpolación (optimizada)
        else {
            // ⚡ Búsqueda binaria optimizada con hint del último índice
            const auto& areas = m_vEstuaryAreas[i];
            auto it = std::lower_bound(areas.begin(), areas.end(), currentArea);
            const int j = std::distance(areas.begin(), it) - 1;
            
            if (j >= 0 && j < nElevationSectionsNumber-1) {
                // ⚡ Calcular factor una sola vez
                const double denom = areas[j+1] - areas[j];
                const double factor = (currentArea - areas[j]) / denom;
                
                // ⚡ Interpolación lineal vectorizada (compilador puede usar SIMD)
                const double inv_factor = 1.0 - factor;
                m_vCrossSectionHydraulicRadius[i] = inv_factor * m_vEstuaryHydraulicRadius[i][j] + 
                                                    factor * m_vEstuaryHydraulicRadius[i][j+1];
                m_vCrossSectionWaterDepth[i] = inv_factor * m_vEstuaryWaterDepths[i][j] + 
                                               factor * m_vEstuaryWaterDepths[i][j+1];
                m_vCrossSectionWidth[i] = inv_factor * m_vWidth[i][j] + 
                                          factor * m_vWidth[i][j+1];
                m_vCrossSectionBeta[i] = inv_factor * m_vBeta[i][j] + 
                                         factor * m_vBeta[i][j+1];
                
                // Interpolate I1 for predictor/normal states
                const double I1_interpolated = inv_factor * m_vEstuaryI1[i][j] + 
                                               factor * m_vEstuaryI1[i][j+1];
                if (m_nPredictor == 1) {
                    m_vCrossSectionI1[i] = I1_interpolated;
                } else {
                    m_vPredictedCrossSectionI1[i] = I1_interpolated;
                }
            }        
        }
        
        m_vCrossSectionWaterElevation[i] = m_vCrossSectionWaterDepth[i] + m_vBedZ[i];
    }
}
    // else {
    //     for (int i = 0; i < nCrossSections; i++) {
    //         const double dArea = m_vPredictedCrossSectionArea[i];
    //         const int nElevationSectionsNumber = estuary[i].nGetElevationSectionsNumber();
    //
    //         if (estuary[i].dGetArea(0) > dArea) {
    //             //! Take the first node if dArea < the Area of the first elevation node
    //             // vElevationSection.push_back(0);
    //             getFirstHydraulicParameters(i);
    //         }
    //         else if (estuary[i].dGetArea(nElevationSectionsNumber -1) < dArea) {
    //             //! Take the last node if dArea > the Area of the last elevation node
    //             // vElevationSection.push_back(nCrossSections - 1);
    //             getLastHydraulicParameters(i);
    //         }
    //         else
    //         {
    //             for (int j = 0; j < nElevationSectionsNumber; j++) {
    //                 //! Getting the elevations node which Area is below dArea and next node Area higher than dArea
    //                 if ((estuary[i].dGetArea(j) <= dArea) && (estuary[i].dGetArea(j+1) > dArea)) {
    //                     // vElevationSection.push_back(j);
    //                     interpolateHydraulicParameters(dArea, i, j);
    //                     break;
    //                 }
    //             }
    //         }
    //         m_vCrossSectionWaterElevation[i] =  m_vCrossSectionWaterDepth[i] + estuary[i].dGetZ();
    //     }
    // }

//======================================================================================================================
//! Interpolate the Hydraulic parameters between the elevation node given.
//======================================================================================================================
void CSimulation::interpolateHydraulicParameters(double dArea, int nCrossSection, int nElevationNode) {
    const double dInterpolationFactor = (dArea - estuary[nCrossSection].dGetArea(nElevationNode)) / (estuary[nCrossSection].dGetArea(nElevationNode+1) - estuary[nCrossSection].dGetArea(nElevationNode));
    m_vCrossSectionHydraulicRadius[nCrossSection] = dInterpolationFactor*(estuary[nCrossSection].dGetHydraulicRadius(nElevationNode+1) - estuary[nCrossSection].dGetHydraulicRadius(nElevationNode)) + estuary[nCrossSection].dGetHydraulicRadius(nElevationNode);
    m_vCrossSectionWidth[nCrossSection] =  dInterpolationFactor*(estuary[nCrossSection].dGetWidth(nElevationNode+1) - estuary[nCrossSection].dGetWidth(nElevationNode)) + estuary[nCrossSection].dGetWidth(nElevationNode);
    m_vCrossSectionWaterDepth[nCrossSection] = dInterpolationFactor*(estuary[nCrossSection].dGetWaterDepth(nElevationNode+1) - estuary[nCrossSection].dGetWaterDepth(nElevationNode)) + estuary[nCrossSection].dGetWaterDepth(nElevationNode);
    m_vCrossSectionBeta[nCrossSection] = dInterpolationFactor*(estuary[nCrossSection].dGetBeta(nElevationNode+1) - estuary[nCrossSection].dGetBeta(nElevationNode)) + estuary[nCrossSection].dGetBeta(nElevationNode);
    m_vCrossSectionLeftRBLocation[nCrossSection] = dInterpolationFactor*(estuary[nCrossSection].dGetLeftY(nElevationNode+1) - estuary[nCrossSection].dGetLeftY(nElevationNode)) + estuary[nCrossSection].dGetLeftY(nElevationNode);
    m_vCrossSectionRightRBLocation[nCrossSection] = dInterpolationFactor*(estuary[nCrossSection].dGetRightY(nElevationNode+1) - estuary[nCrossSection].dGetRightY(nElevationNode)) + estuary[nCrossSection].dGetRightY(nElevationNode);
}

//======================================================================================================================
//! Get first node hydraulic parameters
//======================================================================================================================
void CSimulation::getFirstHydraulicParameters(const int nCrossSection) {
    m_vCrossSectionHydraulicRadius[nCrossSection] = estuary[nCrossSection].dGetHydraulicRadius(0);
    m_vCrossSectionWidth[nCrossSection] =  estuary[nCrossSection].dGetWidth(0);
    m_vCrossSectionWaterDepth[nCrossSection] = estuary[nCrossSection].dGetWaterDepth(0);
    m_vCrossSectionBeta[nCrossSection] = estuary[nCrossSection].dGetBeta(0);
    m_vCrossSectionLeftRBLocation[nCrossSection] = estuary[nCrossSection].dGetLeftY(0);
    m_vCrossSectionRightRBLocation[nCrossSection] = estuary[nCrossSection].dGetRightY(0);
}


//======================================================================================================================
//! Get last node hydraulic parameters
//======================================================================================================================
void CSimulation::getLastHydraulicParameters(const int nCrossSection) {
    const int nLastNode = estuary[nCrossSection].nGetElevationSectionsNumber() - 1;
     m_vCrossSectionHydraulicRadius[nCrossSection] = estuary[nCrossSection].dGetHydraulicRadius(nLastNode);
     m_vCrossSectionWidth[nCrossSection] =  estuary[nCrossSection].dGetWidth(nLastNode);
     m_vCrossSectionWaterDepth[nCrossSection] = estuary[nCrossSection].dGetWaterDepth(nLastNode);
     m_vCrossSectionBeta[nCrossSection] = estuary[nCrossSection].dGetBeta(nLastNode);
     m_vCrossSectionLeftRBLocation[nCrossSection] = estuary[nCrossSection].dGetLeftY(nLastNode);
     m_vCrossSectionRightRBLocation[nCrossSection] = estuary[nCrossSection].dGetRightY(nLastNode);
}


//======================================================================================================================
//! Compute the Courant number
//======================================================================================================================
void CSimulation::calculateTimestep() {
    double dMinTimestep = 10000000.0;
    double dWaterDensityFactor = 1.0;
  
    // ⚡ Calcular timestep solo en nodos interiores (excluir fronteras i=0 y i=n-1)
    for (int i=1; i < m_nCrossSectionsNumber-1; i++) {
        if (m_vCrossSectionArea[i] > DRY_AREA) {
            double dX = m_vCrossSectionDX[i];
            //! Mean water flow
            m_vCrossSectionU[i] = m_vCrossSectionQ[i]/m_vCrossSectionArea[i];

            //! Mean wave celerity
            m_vCrossSectionC[i] = sqrt(G * m_vCrossSectionArea[i] / m_vCrossSectionWidth[i]);

            //! Courant number accounting for density effects (future: baroclinic wave celerity)
            if (m_bDoWaterDensity) {
                dWaterDensityFactor = m_dCourantNumber*dX/(m_vCrossSectionWidth[i]*m_dBetaSalinityConstant);
            }

            //! Compute the timestep given the Courant number
            if (const double dTimestepTmp =  m_dCourantNumber * dX / 
                     std::max(fabs(m_vCrossSectionU[i] + m_vCrossSectionC[i]), 
                             fabs(m_vCrossSectionU[i] - m_vCrossSectionC[i])); dTimestepTmp < dMinTimestep) {
                dMinTimestep = dTimestepTmp;
            }

            if (m_bDoWaterDensity && (dWaterDensityFactor < dMinTimestep)) {
                dMinTimestep = dWaterDensityFactor;
            }
        }
        else {
            double dX = m_vCrossSectionDX[i];
            //! Mean water flow
            m_vCrossSectionU[i] = DRY_Q/DRY_AREA;

            //! ✅ Celeridad para zona seca (profundidad mínima)
            double h_dry = DRY_AREA / m_vCrossSectionWidth[i];
            m_vCrossSectionC[i] = sqrt(G * h_dry);

            const double dTimestepTmp =  m_dCourantNumber * dX / 
                     std::max(fabs(m_vCrossSectionU[i] + m_vCrossSectionC[i]), 
                             fabs(m_vCrossSectionU[i] - m_vCrossSectionC[i]));
            if (dTimestepTmp < dMinTimestep) {
                dMinTimestep = dTimestepTmp;
        }
    }}
    
    // ⚡ Calcular u y c en nodos frontera (pero no los usamos para dt)
    for (int i : {0, m_nCrossSectionsNumber-1}) {
        if (m_vCrossSectionArea[i] != DRY_AREA) {
            m_vCrossSectionU[i] = m_vCrossSectionQ[i]/m_vCrossSectionArea[i];
            m_vCrossSectionC[i] = sqrt(G * m_vCrossSectionArea[i] / m_vCrossSectionWidth[i]);
        } else {
            m_vCrossSectionU[i] = DRY_Q/DRY_AREA;
            double h_dry = DRY_AREA / m_vCrossSectionWidth[i];
            m_vCrossSectionC[i] = sqrt(G * h_dry);
        }
    }
    
    m_dTimestep = dMinTimestep;

    //! Check if m_dTimestep is lower than 0.1 seconds
    if (m_dTimestep < 0.1) {
        m_dTimestep = 0.1;
    }
    
    //! Check if m_dTimestep is higher than save timestep
    if (m_dTimestep > m_dSimTimestep) {
        m_dTimestep = m_dSimTimestep;
    }

    //! Check if the timestep must be reduced to the following save step
    if (m_bSaveAllTimesteps) {
        // Save at every computational timestep
        m_bSaveTime = true;
    }
    else if ((m_nTimeId < static_cast<int>(m_vOutputTimes.size())) && (m_dCurrentTime + m_dTimestep > m_vOutputTimes[m_nTimeId])) {
        double dt_to_save = m_vOutputTimes[m_nTimeId] - m_dCurrentTime;
        // Evitar dt=0 cuando ya estamos en tiempo de guardado
        if (dt_to_save > 1e-6) {
            m_dTimestep = dt_to_save;
        }
        m_bSaveTime = true;
    }
    else {
            m_bSaveTime = false;
    }
    m_dLambda = m_dTimestep / dMinVectorValue(m_vCrossSectionDX);
}


//===============================================================================================================================
//! Calculate boundary condition values at current time t by interpolating from time series
//! Converts elevation boundary conditions to cross-sectional areas using geometry tables
//===============================================================================================================================
void CSimulation::calculateBoundaryConditions() {
    //! Upstream boundary: interpolate elevation and convert to area
    if (nGetUpwardEstuarineCondition() == 2) {
        double dUpwardBoundaryValue = linearInterpolation1d(m_dCurrentTime, m_vUpwardBoundaryConditionTime, m_vUpwardBoundaryConditionValue);
        //! Convert water surface elevation to cross-sectional area at node 0
        m_dUpwardBoundaryValue = linearInterpolation1d(-m_vBedZ[0] + dUpwardBoundaryValue, 
                                                       m_vEstuaryWaterDepths[0],
                                                       m_vEstuaryAreas[0]);
        //! Convert for adjacent node 1 (for smooth transition)
        m_dNextUpwardBoundaryValue = linearInterpolation1d(-m_vBedZ[1] + dUpwardBoundaryValue, 
                                                          m_vEstuaryWaterDepths[1],
                                                          m_vEstuaryAreas[1]);
    }

    //! Downstream boundary: interpolate elevation and convert to area
    // if (nGetDownwardEstuarineCondition() == 1) {
    //     //! Water flow given upstream
    //     m_dDownwardBoundaryValue = linearInterpolation1d(m_dCurrentTime, m_vDownwardBoundaryConditionTime, m_vDownwardBoundaryConditionValue);
    //     m_dNextDownwardBoundaryValue = m_dDownwardBoundaryValue;
    // }
    if (nGetDownwardEstuarineCondition() == 2) {
        //! Elevation given upstream
        double dDownwardBoundaryValue = linearInterpolation1d(m_dCurrentTime, m_vDownwardBoundaryConditionTime, m_vDownwardBoundaryConditionValue);
        m_dDownwardBoundaryValue = linearInterpolation1d(-m_vBedZ[m_nCrossSectionsNumber-1] + dDownwardBoundaryValue, 
                                                         m_vEstuaryWaterDepths[m_nCrossSectionsNumber-1], 
                                                         m_vEstuaryAreas[m_nCrossSectionsNumber-1]);
        m_dNextDownwardBoundaryValue = linearInterpolation1d(-m_vBedZ[m_nCrossSectionsNumber-2] + dDownwardBoundaryValue, 
                                                            m_vEstuaryWaterDepths[m_nCrossSectionsNumber-2], 
                                                            m_vEstuaryAreas[m_nCrossSectionsNumber-2]);
    }
    //     double dDownwardBoundaryValue = linearInterpolation1d(m_dCurrentTime, m_vDownwardBoundaryConditionTime, m_vDownwardBoundaryConditionValue);
    //     m_dDownwardBoundaryValue = linearInterpolation1d(-estuary[m_nCrossSectionsNumber-1].dGetZ() + dDownwardBoundaryValue, estuary[m_nCrossSectionsNumber-1].vGetWaterDepth(), estuary[m_nCrossSectionsNumber-1].vGetArea());
    //     m_dNextDownwardBoundaryValue = linearInterpolation1d(-estuary[m_nCrossSectionsNumber-2].dGetZ() + dDownwardBoundaryValue, estuary[m_nCrossSectionsNumber-2].vGetWaterDepth(), estuary[m_nCrossSectionsNumber-2].vGetArea());
    // }
    // else {};

    //! Hydrographs
    if (m_nHydrographsNumber > 0) {
        for (int i = 0; i < m_nHydrographsNumber; i++) {
            m_vLateralSourcesAtT[hydrographs[i].m_nNearestCrossSectionNo] = linearInterpolation1d(m_dCurrentTime, hydrographs[i].vGetTime(), hydrographs[i].vGetQ());
        }
    }

}

//===============================================================================================================================
//! Change Area and Water flow if soil is dry
//===============================================================================================================================
void CSimulation::dryArea()
{
    if (m_nPredictor == 1) {     
    for (int i = 0; i < m_nCrossSectionsNumber; i++)
        {
            //! Define the dry area based on cross-section geometry
            if (m_vPredictedCrossSectionArea[i] <= DRY_AREA) {
                m_vPredictedCrossSectionArea[i] = DRY_AREA;
                m_vPredictedCrossSectionQ[i] = DRY_Q;
            }
        }
    }
    else if (m_nPredictor == 2) {
        for (int i = 0; i < m_nCrossSectionsNumber; i++)
        {
            //! Define the dry area based on cross-section geometry
            if (m_vCorrectedCrossSectionArea[i] <= DRY_AREA) {
                m_vCorrectedCrossSectionArea[i] = DRY_AREA;
                m_vCorrectedCrossSectionQ[i] = DRY_Q;
            }
        }
    }
    else {
        for (int i = 0; i < m_nCrossSectionsNumber; i++)
        {
            //! Define the dry area based on cross-section geometry
            if (m_vCrossSectionArea[i] <= DRY_AREA) {
                m_vCrossSectionArea[i] = DRY_AREA;
                m_vCrossSectionQ[i] = DRY_Q;
            }
        }
    }

}

//===============================================================================================================================
//! Change equation terms if soil is dry
//===============================================================================================================================
void CSimulation::dryTerms()
{
    if (m_nPredictor == 1) {
    for (int i = 0; i < m_nCrossSectionsNumber-1; i++)
        {
            //! Define the dry area based on cross-section geometry
            if (m_vCrossSectionArea[i] <= DRY_AREA) {
                m_vCrossSectionF0[i] = 0.0;
                m_vCrossSectionF1[i] = 0.0;
                m_vCrossSectionGv1[i] = 0.0;
            }
        }
    }
    else {      
        for (int i = 0; i < m_nCrossSectionsNumber-1; i++)
        {
            //! Define the dry area based on cross-section geometry
            if (m_vPredictedCrossSectionArea[i] <= DRY_AREA) {
                m_vCrossSectionF0[i] = 0.0;
                m_vCrossSectionF1[i] = 0.0;
                m_vCrossSectionGv1[i] = 0.0;
            }
        }
    }

}

//===============================================================================================================================
//! Do Murillo condition
//===============================================================================================================================
void CSimulation::doMurilloCondition()
{
    if (m_nPredictor != 1) {
        for (int i = 0; i < m_nCrossSectionsNumber; i++)
        {
            const double dValue = CDX_MURILLO * sqrt(
                                      2 * pow(m_vCrossSectionHydraulicRadius[i], 2.0 / 3.0) / (G * m_vCrossSectionDX[i]));

            if (dValue < estuary[i].dGetManningNumber()) {
                m_vCrossSectionMurilloFactor[i] = dValue/m_vCrossSectionManningNumber[i];
                m_vCrossSectionManningNumber[i] = dValue;
            }
            else {
                m_vCrossSectionMurilloFactor[i] = 1.0;
            }
        }
    }
}

//======================================================================================================================
//! Calculate bed slope source terms (gAS0) and friction terms (gASf)
//! Uses upwind bed slopes to ensure consistency with predictor-corrector discretization
//======================================================================================================================
void CSimulation::calculate_GS_A_terms() {
    
    if (bGetDoSurfaceTermBalance()) {
        //! Average hydraulic variables between adjacent cross-sections for source term balance
        //! Predictor: average with i+1 (forward), Corrector: average with i-1 (backward)      
        for (int i = 0; i < m_nCrossSectionsNumber; i++) {
            double dMeanArea, dMeanQ, dMeanHydraulicRadius;
            
            if (m_nPredictor == 1) {
                if (i < m_nCrossSectionsNumber - 1) {
                    dMeanArea = (m_vCrossSectionArea[i] + m_vCrossSectionArea[i+1])/2.0;
                    dMeanQ = (m_vCrossSectionQ[i] + m_vCrossSectionQ[i+1])/2.0;
                    dMeanHydraulicRadius = (m_vCrossSectionHydraulicRadius[i] + m_vCrossSectionHydraulicRadius[i+1])/2.0;
                } else {
                    dMeanArea = m_vCrossSectionArea[i];
                    dMeanQ = m_vCrossSectionQ[i];
                    dMeanHydraulicRadius = m_vCrossSectionHydraulicRadius[i];
                }
            } else {
                if (i > 0) {
                    dMeanArea = (m_vPredictedCrossSectionArea[i] + m_vPredictedCrossSectionArea[i-1])/2.0;
                    dMeanQ = (m_vPredictedCrossSectionQ[i] + m_vPredictedCrossSectionQ[i-1])/2.0;
                    dMeanHydraulicRadius = (m_vCrossSectionHydraulicRadius[i] + m_vCrossSectionHydraulicRadius[i-1])/2.0;
                } else {
                    dMeanArea = m_vPredictedCrossSectionArea[i];
                    dMeanQ = m_vPredictedCrossSectionQ[i];
                    dMeanHydraulicRadius = m_vCrossSectionHydraulicRadius[i];
                }
            }
            
            // ✅ CRÍTICO: Usar pendientes unilaterales para balance de términos fuente (como Fortran)
            // Predictor: diferencia hacia adelante, Corrector: diferencia hacia atrás
            double S0_to_use;
            if (m_nPredictor == 1) {
                S0_to_use = m_vCrossSectionBedSlopePredictor[i];  // zmedp en Fortran
            } else {
                S0_to_use = m_vCrossSectionBedSlopeCorrector[i];  // zmedc en Fortran
            }
            
            // gAS0 = g * A * S0 (con pendiente unilateral apropiada)
            m_vCrossSection_gAS0[i] = G * dMeanArea * S0_to_use;
            
            // ✅ CORRECTO: Pendiente de fricción según Manning
            if (dMeanArea > DRY_AREA && dMeanHydraulicRadius > 1e-6) {
                m_vCrossSectionFrictionSlope[i] = m_vManningNumberSquared[i] * 
                                                 dMeanQ * fabs(dMeanQ) / 
                                                 (dMeanArea * dMeanArea * pow(dMeanHydraulicRadius, 4.0/3.0));
                m_vCrossSection_gASf[i] = G * dMeanArea * m_vCrossSectionFrictionSlope[i];
            } else {
                m_vCrossSectionFrictionSlope[i] = 0.0;
                m_vCrossSection_gASf[i] = 0.0;
            }
        }
    }
    else {
        // Sin promediado (más simple pero menos preciso)
        for (int i = 0; i < m_nCrossSectionsNumber; i++) {
            // ✅ CORRECTO: S0 ya tiene el signo correcto, no necesita factor de dirección
            m_vCrossSection_gAS0[i] = G * m_vCrossSectionArea[i] * m_vCrossSectionBedSlope[i];
            
            if (m_nPredictor == 1) {
                if (m_vCrossSectionArea[i] > DRY_AREA && m_vCrossSectionHydraulicRadius[i] > 1e-6) {
                    m_vCrossSectionFrictionSlope[i] = m_vManningNumberSquared[i] * 
                                                     m_vCrossSectionQ[i] * fabs(m_vCrossSectionQ[i]) / 
                                                     (m_vCrossSectionArea[i] * m_vCrossSectionArea[i] * pow(m_vCrossSectionHydraulicRadius[i], 4.0/3.0));
                    m_vCrossSection_gASf[i] = G * m_vCrossSectionArea[i] * m_vCrossSectionFrictionSlope[i];
                } else {
                    m_vCrossSectionFrictionSlope[i] = 0.0;
                    m_vCrossSection_gASf[i] = 0.0;
                }
            }
            else {
                if (m_vPredictedCrossSectionArea[i] > DRY_AREA && m_vCrossSectionHydraulicRadius[i] > 1e-6) {
                    m_vCrossSectionFrictionSlope[i] = m_vManningNumberSquared[i] * 
                                                     m_vPredictedCrossSectionQ[i] * fabs(m_vPredictedCrossSectionQ[i]) / 
                                                     (m_vPredictedCrossSectionArea[i] * m_vPredictedCrossSectionArea[i] * pow(m_vCrossSectionHydraulicRadius[i], 4.0/3.0));
                    m_vCrossSection_gASf[i] = G * m_vPredictedCrossSectionArea[i] * m_vCrossSectionFrictionSlope[i];
                } else {
                    m_vCrossSectionFrictionSlope[i] = 0.0;
                    m_vCrossSection_gASf[i] = 0.0;
                }
            }
        }
    }
}

//======================================================================================================================
//! Calculate F terms
//======================================================================================================================
void CSimulation::calculateFlowTerms() {
    
    if (m_nPredictor == 1) {      
        for (int i = 0; i < m_nCrossSectionsNumber; i++) {
            m_vCrossSectionF0[i] = m_vCrossSectionQ[i];
            
            // F1 = β*Q²/A (I1 DISABLED - no stable scaling factor found)
            // Fortran: F1 = β*(Q²/A + g*I1)
            // Problem: g*I1 >> Q²/A causes instability, but scaling factors alter physics
            if (m_vCrossSectionArea[i] > DRY_AREA) {
                m_vCrossSectionF1[i] = m_vCrossSectionBeta[i] * pow(m_vCrossSectionQ[i], 2.0) / m_vCrossSectionArea[i];
            } else {
                m_vCrossSectionF1[i] = 0.0;
            }
        }
    }
    else {
        for (int i = 0; i < m_nCrossSectionsNumber; i++) {
            m_vCrossSectionF0[i] = m_vPredictedCrossSectionQ[i];
            
            // F1 = β*Q²/A para corrector (I1 DISABLED)
            if (m_vPredictedCrossSectionArea[i] > DRY_AREA) {
                m_vCrossSectionF1[i] = m_vCrossSectionBeta[i] * pow(m_vPredictedCrossSectionQ[i], 2.0) / m_vPredictedCrossSectionArea[i];
            } else {
                m_vCrossSectionF1[i] = 0.0;
            }
        }

    }
}

//======================================================================================================================
//! Calculate source terms for momentum equation: gravity, friction, and pressure gradients
//! Uses adaptive upwind scheme for water surface gradients at locations with steep bed slope changes
//======================================================================================================================
void CSimulation::calculateSourceTerms() {
    //! Calculate water surface gradient dh/dx using central differences
    //! NOTA: El esquema upwind adaptativo causaba inestabilidad en zonas con depresiones/elevaciones del lecho
    //! donde hay cambios de signo en la pendiente. Usar siempre diferencias centrales es más robusto.
    for (int i = 1; i < m_nCrossSectionsNumber-1; i++) {
        //! Use central differences with precomputed inverse for all nodes
        m_vCrossSectionDhDx[i] = (m_vCrossSectionWaterElevation[i+1] - m_vCrossSectionWaterElevation[i-1]) * m_vInvDxSum[i];
        
        /* ESQUEMA UPWIND DESACTIVADO - Causaba inestabilidad con pendiente variable
        double dx1 = m_vPositionX[i+1] - m_vPositionX[i];
        double dx2 = m_vPositionX[i] - m_vPositionX[i-1];
        
        //! Detect abrupt bed slope changes to apply upwind differencing
        double S0_left = (estuary[i].dGetZ() - estuary[i-1].dGetZ()) / dx2;
        double S0_right = (estuary[i+1].dGetZ() - estuary[i].dGetZ()) / dx1;
        double slope_change = fabs(S0_right - S0_left);
        
        //! If bed slope change > 0.001 (0.1%), use upwind scheme based on flow direction
        if (slope_change > 0.001) {
            double Q_current = (m_nPredictor == 1) ? m_vCrossSectionQ[i] : m_vPredictedCrossSectionQ[i];
            
            if (Q_current > 0) {
                //! Positive flow (downstream): use backward difference
                m_vCrossSectionDhDx[i] = (m_vCrossSectionWaterElevation[i] - m_vCrossSectionWaterElevation[i-1]) / dx2;
            } else {
                //! Negative flow (upstream): use forward difference
                m_vCrossSectionDhDx[i] = (m_vCrossSectionWaterElevation[i+1] - m_vCrossSectionWaterElevation[i]) / dx1;
            }
        } else {
            //! Smooth bed region: use central differences with precomputed inverse
            m_vCrossSectionDhDx[i] = (m_vCrossSectionWaterElevation[i+1] - m_vCrossSectionWaterElevation[i-1]) * m_vInvDxSum[i];
        }
        */
    }
    
    // Condiciones de frontera con valores precalculados
    m_vCrossSectionDhDx[0] = (m_vCrossSectionWaterElevation[1] - m_vCrossSectionWaterElevation[0]) * m_vInvDX[0];
    m_vCrossSectionDhDx[m_nCrossSectionsNumber-1] = (m_vCrossSectionWaterElevation[m_nCrossSectionsNumber-1] - m_vCrossSectionWaterElevation[m_nCrossSectionsNumber-2]) * m_vInvDX[m_nCrossSectionsNumber-2];

    //! NOTE: I1 gradient term (dI1/dx) for non-prismatic channels is DISABLED
    //! Causes numerical instability even with proper calculation. The large I1 gradients
    //! in highly non-prismatic sections overwhelm the scheme stability.
    //! For future: consider flux-limiting or implicit treatment of I1 term.
    /*
    // Seleccionar vectores de I1 según estado predictor/corrector
    const auto& I1_to_use = (m_nPredictor == 1) ? m_vCrossSectionI1 : m_vPredictedCrossSectionI1;
    
    // Calcular dI1/dx con diferencias centradas (nodos interiores)
    for (int i = 1; i < m_nCrossSectionsNumber-1; i++) {
        double dx1 = m_vPositionX[i+1] - m_vPositionX[i];
        double dx2 = m_vPositionX[i] - m_vPositionX[i-1];
        
        // Para canales no prismáticos: usar dI1/dx directamente
        // (incluye variación de geometría en x, no solo de h)
        m_vCrossSectionDI1Dx[i] = (I1_to_use[i+1] - I1_to_use[i-1]) / (dx1 + dx2);
    }
    
    // Condiciones de frontera para dI1/dx
    m_vCrossSectionDI1Dx[0] = (I1_to_use[1] - I1_to_use[0]) * m_vInvDX[0];
    m_vCrossSectionDI1Dx[m_nCrossSectionsNumber-1] = 
        (I1_to_use[m_nCrossSectionsNumber-1] - I1_to_use[m_nCrossSectionsNumber-2]) * m_vInvDX[m_nCrossSectionsNumber-2];
    */

    // Calcular términos fuente
    for (int i = 0; i < m_nCrossSectionsNumber; i++) {
        double dx = 0.0;
        if (i < m_nCrossSectionsNumber - 1) {
            dx = m_vPositionX[i+1] - m_vPositionX[i];
        }
        else {
            dx = m_vPositionX[i] - m_vPositionX[i-1];
        }
        
        // ⚡ Usar 1/dx precalculado cuando sea el mismo índice
        double inv_dx = (i < m_nCrossSectionsNumber - 1) ? m_vInvDX[i] : m_vInvDX[i-1];
        m_vCrossSectionGv0[i] = m_vLateralSourcesAtT[i] * inv_dx;
        
        // Gv1 = -g*A*dh/dx + gAS0 - gASf (aproximación hidrostática)
        // NOTA: Para canales no prismáticos, lo correcto sería -g*dI1/dx
        // pero causa inestabilidad numérica. La aproximación -g*A*dh/dx es estable
        // y válida cuando la geometría varía suavemente en x
        double area_to_use = (m_nPredictor == 1) ? m_vCrossSectionArea[i] : m_vPredictedCrossSectionArea[i];
        m_vCrossSectionGv1[i] = -G * area_to_use * m_vCrossSectionDhDx[i] + 
                               m_vCrossSection_gAS0[i] - m_vCrossSection_gASf[i];
        
        // Versión con dI1/dx (DESACTIVADA - inestable):
        // m_vCrossSectionGv1[i] = -G * m_vCrossSectionDI1Dx[i] + 
        //                        m_vCrossSection_gAS0[i] - m_vCrossSection_gASf[i];
    }
}


//======================================================================================================================
//! Calculate the predictor
//======================================================================================================================
void CSimulation::calculatePredictor() {
    //! Predictor step of McCormack scheme: forward differences
    //! Computes only interior nodes (i=1 to i=n-2), boundaries handled separately
    //! Los nodos i=0 y i=n-1 se manejan exclusivamente por las BC
    
    // ⚡ OPTIMIZACIÓN: Loop optimizado por compilador con -O3 -march=native
    const int n = m_nCrossSectionsNumber - 1;
    for (int i = 1; i < n; i++) {
        const double rho = m_vCrossSectionRho[i];
        const double inv_rho = 1.0 / rho;
        
        m_vPredictedCrossSectionArea[i] = (m_vCrossSectionArea[i]*rho - 
            m_dLambda*(m_vCrossSectionF0[i+1]*m_vCrossSectionRho[i+1] - m_vCrossSectionF0[i]*rho) + 
            m_dTimestep*m_vCrossSectionGv0[i]*rho) * inv_rho;
            
        m_vPredictedCrossSectionQ[i] = (m_vCrossSectionQ[i]*rho - 
            m_dLambda*(m_vCrossSectionF1[i+1]*m_vCrossSectionRho[i+1] - m_vCrossSectionF1[i]*rho) + 
            m_dTimestep*m_vCrossSectionGv1[i]*rho) * inv_rho;
    }
}


//======================================================================================================================
//! Calculate the corrector
//======================================================================================================================
void CSimulation::calculateCorrector() {
    //! Corrector step of McCormack scheme: backward differences
    //! Computes only interior nodes (i=1 to i=n-2), boundaries handled separately
    //! Los nodos i=0 y i=n-1 se manejan exclusivamente por las BC
    
    // ⚡ OPTIMIZACIÓN: Loop optimizado por compilador con -O3 -march=native
    const int n = m_nCrossSectionsNumber - 1;
    for (int i = 1; i < n; i++) {
        const double rho = m_vCrossSectionRho[i];
        const double inv_rho = 1.0 / rho;
        
        m_vCorrectedCrossSectionArea[i] = (m_vPredictedCrossSectionArea[i]*rho - 
            m_dLambda*(m_vCrossSectionF0[i]*rho - m_vCrossSectionF0[i-1]*m_vCrossSectionRho[i-1]) + 
            m_dTimestep*m_vCrossSectionGv0[i]*rho) * inv_rho;
            
        m_vCorrectedCrossSectionQ[i] = (m_vPredictedCrossSectionQ[i]*rho - 
            m_dLambda*(m_vCrossSectionF1[i]*rho - m_vCrossSectionF1[i-1]*m_vCrossSectionRho[i-1]) + 
            m_dTimestep*m_vCrossSectionGv1[i]*rho) * inv_rho;
    }
}

//======================================================================================================================
//! Update Predictor boundaries - CONDICIONES CON SUAVIZADO SELECTIVO
//======================================================================================================================
void CSimulation::updatePredictorBoundaries() {
    //==============================================================================================================
    //! CONDICIÓN DE FRONTERA AGUAS ARRIBA (UPSTREAM)
    //==============================================================================================================
    
    if (nGetUpwardEstuarineCondition() == 0) {
        //! ✅ CONDICIÓN ABIERTA (FREE): Extrapolación de gradiente cero
        //! Extrapolar linealmente pero con limitador para estabilidad
        const double dQ = m_vCrossSectionQ[2] - m_vCrossSectionQ[1];
        const double dA = m_vCrossSectionArea[2] - m_vCrossSectionArea[1];
        
        // Extrapolación lineal: Q[0] = Q[1] + (Q[1]-Q[2])
        m_vPredictedCrossSectionQ[0] = m_vCrossSectionQ[1] - dQ;
        m_vPredictedCrossSectionArea[0] = m_vCrossSectionArea[1] - dA;
        
        // Limitar área mínima
        if (m_vPredictedCrossSectionArea[0] < DRY_AREA) {
            m_vPredictedCrossSectionArea[0] = m_vCrossSectionArea[1];
        }
    }
    else if (nGetUpwardEstuarineCondition() == 1) {
        //! ✅ CONDICIÓN REFLECTANTE (WALL): Q = 0 (pared sólida)
        //! Imponer velocidad cero, extrapolar área (nivel de agua) desde interior
        m_vPredictedCrossSectionQ[0] = 0.0;
        
        // Extrapolar área manteniendo gradiente desde interior
        const double dA = m_vCrossSectionArea[2] - m_vCrossSectionArea[1];
        m_vPredictedCrossSectionArea[0] = m_vCrossSectionArea[1] - dA;
        
        // Limitar área mínima
        if (m_vPredictedCrossSectionArea[0] < DRY_AREA) {
            m_vPredictedCrossSectionArea[0] = m_vCrossSectionArea[1];
        }
    }
    else if (nGetUpwardEstuarineCondition() == 2) {
        //! ✅ ELEVACIÓN IMPUESTA: Imponer área y calcular Q con Manning
        m_vPredictedCrossSectionArea[0] = m_dUpwardBoundaryValue;
        m_vCrossSectionHydraulicRadius[0] = linearInterpolation1d(
            m_vPredictedCrossSectionArea[0], 
            m_vEstuaryAreas[0], 
            m_vEstuaryHydraulicRadius[0]);
        
        double sign_S0 = (m_vCrossSectionBedSlope[0] >= 0) ? 1.0 : -1.0;
        m_vPredictedCrossSectionQ[0] = m_vPredictedCrossSectionArea[0] * 
                                      sqrt(fabs(m_vCrossSectionBedSlope[0]) + 1e-10) * 
                                      pow(m_vCrossSectionHydraulicRadius[0], 2.0/3.0) * 
                                      sign_S0 / (m_vManningN[0] + 1e-10);
    }
    else {
        //! ✅ CAUDAL IMPUESTO
        m_vPredictedCrossSectionQ[0] = m_dUpwardBoundaryValue;
        
        double dManningFactor = m_vPredictedCrossSectionQ[0] * m_vManningN[0] / 
                              (sqrt(fabs(m_vCrossSectionBedSlope[0]) + 1e-10));
        
        m_vPredictedCrossSectionArea[0] = linearInterpolation1d(
            dManningFactor, 
            m_vPrecalculatedSecondTerm[0], 
            m_vEstuaryAreas[0]);
    }

    //==============================================================================================================
    //! HIDROGRAMAS LATERALES
    //==============================================================================================================
    for (int i = 0; i < nGetHydrographsNumber(); i++) {
        int node = hydrographs[i].m_nNearestCrossSectionNo;
        if (node >= 0 && node < m_nCrossSectionsNumber) {
            m_vPredictedCrossSectionQ[node] += m_vLateralSourcesAtT[node];
        }
    }

    //==============================================================================================================
    //! CONDICIÓN DE FRONTERA AGUAS ABAJO (DOWNSTREAM)
    //==============================================================================================================
    int n = m_nCrossSectionsNumber - 1;
    
    if (nGetDownwardEstuarineCondition() == 0) {
        //! ✅ CONDICIÓN ABIERTA (FREE): Extrapolación de gradiente cero
        const double dQ = m_vCrossSectionQ[n-1] - m_vCrossSectionQ[n-2];
        const double dA = m_vCrossSectionArea[n-1] - m_vCrossSectionArea[n-2];
        
        m_vPredictedCrossSectionQ[n] = m_vCrossSectionQ[n-1] + dQ;
        m_vPredictedCrossSectionArea[n] = m_vCrossSectionArea[n-1] + dA;
        
        if (m_vPredictedCrossSectionArea[n] < DRY_AREA) {
            m_vPredictedCrossSectionArea[n] = m_vCrossSectionArea[n-1];
        }
    }
    else if (nGetDownwardEstuarineCondition() == 1) {
        //! ✅ CAUDAL IMPUESTO
        m_vPredictedCrossSectionQ[n] = m_dDownwardBoundaryValue;
        
        const double dManningFactor = m_vPredictedCrossSectionQ[n] * m_vManningN[n] / 
                                     (sqrt(fabs(m_vCrossSectionBedSlope[n]) + 1e-10));
        
        m_vPredictedCrossSectionArea[n] = linearInterpolation1d(
            dManningFactor, 
            m_vPrecalculatedSecondTerm[n],
            m_vEstuaryAreas[n]);
    }
    else {
        //! ✅ ELEVACIÓN/MAREA IMPUESTA: Imponer área en n, calcular Q con dh/dt de marea
        m_vPredictedCrossSectionArea[n] = m_dDownwardBoundaryValue;
        
        // Calcular dh/dt de la serie temporal de marea
        double dhdt = 0.0;
        if (m_vDownwardBoundaryConditionTime.size() > 1) {
            // Buscar intervalo actual en la serie de marea
            int idx = 0;
            for (size_t i = 0; i < m_vDownwardBoundaryConditionTime.size() - 1; i++) {
                if (m_dCurrentTime >= m_vDownwardBoundaryConditionTime[i] && 
                    m_dCurrentTime <= m_vDownwardBoundaryConditionTime[i+1]) {
                    idx = i;
                    break;
                }
            }
            
            // Calcular derivada temporal (dh/dt)
            double dt = m_vDownwardBoundaryConditionTime[idx+1] - m_vDownwardBoundaryConditionTime[idx];
            if (dt > 0) {
                double dh = m_vDownwardBoundaryConditionValue[idx+1] - m_vDownwardBoundaryConditionValue[idx];
                dhdt = dh / dt;
            }
        }
        
        // Q = A * dh/dt (conservación de masa en frontera)
        // Si dh/dt > 0: marea subiendo, Q entrante negativo
        // Si dh/dt < 0: marea bajando, Q saliente positivo
        double Q_tide = -m_vPredictedCrossSectionArea[n] * dhdt;
        
        // Combinar con extrapolación desde interior (peso 50/50 para suavizar)
        const double dQ = m_vCrossSectionQ[n-1] - m_vCrossSectionQ[n-2];
        double Q_extrap = m_vCrossSectionQ[n-1] + dQ;
        
        m_vPredictedCrossSectionQ[n] = 0.5 * Q_tide + 0.5 * Q_extrap;
    }


}

//======================================================================================================================
//! Update Corrector boundaries - Con suavizado
//======================================================================================================================
void CSimulation::updateCorrectorBoundaries() {
    //==============================================================================================================
    //! CONDICIÓN DE FRONTERA AGUAS ARRIBA (UPSTREAM) - Corrector
    //==============================================================================================================
    
    if (nGetUpwardEstuarineCondition() == 0) {
        //! ✅ CONDICIÓN ABIERTA (FREE): Extrapolación de gradiente cero
        const double dQ = m_vPredictedCrossSectionQ[2] - m_vPredictedCrossSectionQ[1];
        const double dA = m_vPredictedCrossSectionArea[2] - m_vPredictedCrossSectionArea[1];
        
        m_vCorrectedCrossSectionQ[0] = m_vPredictedCrossSectionQ[1] - dQ;
        m_vCorrectedCrossSectionArea[0] = m_vPredictedCrossSectionArea[1] - dA;
        
        if (m_vCorrectedCrossSectionArea[0] < DRY_AREA) {
            m_vCorrectedCrossSectionArea[0] = m_vPredictedCrossSectionArea[1];
        }
    }
    else if (nGetUpwardEstuarineCondition() == 1) {
        //! ✅ CONDICIÓN REFLECTANTE (WALL): Q = 0
        m_vCorrectedCrossSectionQ[0] = 0.0;
        
        // Extrapolar área desde interior (predictor)
        const double dA = m_vPredictedCrossSectionArea[2] - m_vPredictedCrossSectionArea[1];
        m_vCorrectedCrossSectionArea[0] = m_vPredictedCrossSectionArea[1] - dA;
        
        if (m_vCorrectedCrossSectionArea[0] < DRY_AREA) {
            m_vCorrectedCrossSectionArea[0] = m_vPredictedCrossSectionArea[1];
        }
    }
    else if (nGetUpwardEstuarineCondition() == 2) {
        //! ✅ ELEVACIÓN IMPUESTA: Imponer área y calcular Q con Manning
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
        //! ✅ CAUDAL IMPUESTO CON SPONGE LAYER SELECTIVO (Tipo 3)
        m_vCorrectedCrossSectionQ[0] = m_dUpwardBoundaryValue;
        
        double dManningFactor = m_vCorrectedCrossSectionQ[0] * m_vManningN[0] / 
                              (sqrt(fabs(m_vCrossSectionBedSlope[0]) + 1e-10));
        
        m_vCorrectedCrossSectionArea[0] = linearInterpolation1d(
            dManningFactor, 
            m_vPrecalculatedSecondTerm[0], 
            m_vEstuaryAreas[0]);
        
        // 🎯 SPONGE LAYER solo si discontinuidad > 5%
        if (m_nCrossSectionsNumber >= 3) {
            double Q_diff = fabs(m_vCorrectedCrossSectionQ[0] - m_vCrossSectionQ[1]);
            double Q_avg = 0.5 * fabs(m_vCorrectedCrossSectionQ[0]) + 0.5 * fabs(m_vCrossSectionQ[1]) + 1e-6;
            
            if (Q_diff / Q_avg > 0.05) {
                double alpha = 0.8;
                double Q_computed = m_vCorrectedCrossSectionQ[1];
                double Q_smooth = 0.5*(m_vCorrectedCrossSectionQ[0] + m_vCrossSectionQ[2]);
                m_vCorrectedCrossSectionQ[1] = alpha*Q_computed + (1-alpha)*Q_smooth;
            }
        }
    }

    //==============================================================================================================
    //! HIDROGRAMAS LATERALES
    //==============================================================================================================
    for (int i = 0; i < nGetHydrographsNumber(); i++) {
        int node = hydrographs[i].m_nNearestCrossSectionNo;
        if (node >= 0 && node < m_nCrossSectionsNumber) {
            m_vCorrectedCrossSectionQ[node] += m_vLateralSourcesAtT[node];
        }
    }

    //==============================================================================================================
    //! CONDICIÓN DE FRONTERA AGUAS ABAJO (DOWNSTREAM) - Corrector
    //==============================================================================================================
    int n = m_nCrossSectionsNumber - 1;
    
    if (nGetDownwardEstuarineCondition() == 0) {
        //! ✅ CONDICIÓN ABIERTA (FREE): Extrapolación de gradiente cero
        const double dQ = m_vPredictedCrossSectionQ[n-1] - m_vPredictedCrossSectionQ[n-2];
        const double dA = m_vPredictedCrossSectionArea[n-1] - m_vPredictedCrossSectionArea[n-2];
        
        m_vCorrectedCrossSectionQ[n] = m_vPredictedCrossSectionQ[n-1] + dQ;
        m_vCorrectedCrossSectionArea[n] = m_vPredictedCrossSectionArea[n-1] + dA;
        
        if (m_vCorrectedCrossSectionArea[n] < DRY_AREA) {
            m_vCorrectedCrossSectionArea[n] = m_vPredictedCrossSectionArea[n-1];
        }
    }
    else if (nGetDownwardEstuarineCondition() == 1) {
        //! ✅ CAUDAL IMPUESTO
        m_vCorrectedCrossSectionQ[n] = m_dDownwardBoundaryValue;
        
        const double dManningFactor = m_vCorrectedCrossSectionQ[n] * m_vManningN[n] / 
                                     (sqrt(fabs(m_vCrossSectionBedSlope[n]) + 1e-10));
        
        m_vCorrectedCrossSectionArea[n] = linearInterpolation1d(
            dManningFactor, 
            m_vPrecalculatedSecondTerm[n],
            m_vEstuaryAreas[n]);
    }
    else {
        //! ✅ ELEVACIÓN/MAREA IMPUESTA: Imponer área en n, calcular Q con dh/dt de marea
        m_vCorrectedCrossSectionArea[n] = m_dDownwardBoundaryValue;
        
        // Calcular dh/dt de la serie temporal de marea
        double dhdt = 0.0;
        if (m_vDownwardBoundaryConditionTime.size() > 1) {
            // Buscar intervalo actual en la serie de marea
            int idx = 0;
            for (size_t i = 0; i < m_vDownwardBoundaryConditionTime.size() - 1; i++) {
                if (m_dCurrentTime >= m_vDownwardBoundaryConditionTime[i] && 
                    m_dCurrentTime <= m_vDownwardBoundaryConditionTime[i+1]) {
                    idx = i;
                    break;
                }
            }
            
            // Calcular derivada temporal (dh/dt)
            double dt = m_vDownwardBoundaryConditionTime[idx+1] - m_vDownwardBoundaryConditionTime[idx];
            if (dt > 0) {
                double dh = m_vDownwardBoundaryConditionValue[idx+1] - m_vDownwardBoundaryConditionValue[idx];
                dhdt = dh / dt;
            }
        }
        
        // Q = A * dh/dt (conservación de masa en frontera)
        // Si dh/dt > 0: marea subiendo, Q entrante negativo
        // Si dh/dt < 0: marea bajando, Q saliente positivo
        double Q_tide = -m_vCorrectedCrossSectionArea[n] * dhdt;
        
        // Combinar con extrapolación desde interior (peso 50/50 para suavizar)
        const double dQ = m_vPredictedCrossSectionQ[n-1] - m_vPredictedCrossSectionQ[n-2];
        double Q_extrap = m_vPredictedCrossSectionQ[n-1] + dQ;
        
        m_vCorrectedCrossSectionQ[n] = 0.5 * Q_tide + 0.5 * Q_extrap;
    }
}


//======================================================================================================================
//! Update Boundaries
//======================================================================================================================
void CSimulation::updateBoundaries() {
    //! Legacy boundary update function (called from older code paths)
    //! Modern boundary handling is in updatePredictorBoundaries/updateCorrectorBoundaries
    if (nGetUpwardEstuarineCondition() == 0) {
        //! Open boundary condition
        m_vCrossSectionQ[0] = m_vCrossSectionQ[1];
        //! Corrector does not calculate i = 0, it is imposed
        m_vCrossSectionArea[0] = m_vCrossSectionArea[1];
    }
    else if (nGetUpwardEstuarineCondition() == 1) {
        //! Reflected boundary condition
        m_vCrossSectionQ[0] = m_dUpwardBoundaryValue;
        //! Corrector does not calculate i = 0, it es imposed
        m_vCrossSectionArea[0] = m_vCrossSectionArea[1];
    }
    else {
        //! Water elevation boundary condition (m_dUpwardBoundaryValue already converted to area)
        m_vCrossSectionArea[0] = m_dUpwardBoundaryValue;
        //! Corrector does not calculate i = 0, it es imposed
        m_vCrossSectionQ[0] = m_vCrossSectionQ[1];
    }

    if (nGetDownwardEstuarineCondition() == 0) {
        //! Open flow
        m_vCrossSectionArea[m_nCrossSectionsNumber-1] = m_vCrossSectionArea[m_nCrossSectionsNumber-2];
        m_vCrossSectionQ[m_nCrossSectionsNumber-1] = m_vCrossSectionQ[m_nCrossSectionsNumber-2];
    }
    else if (nGetDownwardEstuarineCondition() == 1) {
        //! Given the water flow
        m_vCrossSectionQ[m_nCrossSectionsNumber-1] = m_dDownwardBoundaryValue;
        m_vCrossSectionArea[m_nCrossSectionsNumber-1] = m_vCrossSectionArea[m_nCrossSectionsNumber-2];
    }
    else {
        //! Given tidal elevation seaward and previous to seaward
        m_vCrossSectionArea[m_nCrossSectionsNumber-1] = m_dDownwardBoundaryValue;
        m_vCrossSectionArea[m_nCrossSectionsNumber-2] = m_dNextDownwardBoundaryValue;
        m_vCrossSectionQ[m_nCrossSectionsNumber-1] = m_vCrossSectionQ[m_nCrossSectionsNumber-2];
    }
}

//======================================================================================================================
//! Merge Predictor and Corrector
//======================================================================================================================
void CSimulation::mergePredictorCorrector() {
    //! Merge predictor and corrector steps using McCormack averaging
    //! Applies TVD flux limiters to maintain monotonicity and avoid spurious oscillations
    if (bGetDoMcComarckLimiterFlux())
    {
        //! Include TVD-McComarck? - Usar vectores pre-alocados (miembros de clase)
        const int nCrossSectionsNumber = m_nCrossSectionsNumber;
        
        // ⚡ Resetear vectores a cero (sin reallocations)
        std::fill(m_vTVD_a1_med.begin(), m_vTVD_a1_med.end(), 0.0);
        std::fill(m_vTVD_a2_med.begin(), m_vTVD_a2_med.end(), 0.0);
        std::fill(m_vTVD_alfa1_med.begin(), m_vTVD_alfa1_med.end(), 0.0);
        std::fill(m_vTVD_alfa2_med.begin(), m_vTVD_alfa2_med.end(), 0.0);
        
        // Referencias para mantener compatibilidad con código existente
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


        for (int i = 0; i < m_nCrossSectionsNumber-1; i++)
        {
            if (alfa1_med[i] == 0.0) {
                r1_med[i] = 1.0;
            }
            else {
                //! Computing ri
                if (a1_med[i] < 0)
                {
                    if (i != m_nCrossSectionsNumber -2)
                    {
                        r1_med[i] = alfa1_med[i+1]/alfa1_med[i];
                    }
                    else
                   
                    {
                        r1_med[i] = 1.0;
                    }

                }
                else if (a1_med[i] == 0.0)
                {
                    r1_med[i] = 1.0;
                }
                else if (a1_med[i] > 0)
                {
                    if (i != 0)
                    {
                        r1_med[i] = alfa1_med[i-1]/alfa1_med[i];
                    }
                    else
                    {
                        r1_med[i] = 1.0;  // ✅ CORREGIDO: era 0, debe ser 1.0 para no desactivar limitador
                    }
                }
            }

            if (alfa2_med[i] == 0.0) {
                r2_med[i] = 1.0;
            }
            else {
                if (a2_med[i] < 0)
                {
                    if (i != m_nCrossSectionsNumber-1)
                    {
                        r2_med[i] = alfa2_med[i+1]/alfa2_med[i];
                    }
                    else
                    {
                        r2_med[i] = 1.0;
                    }
                }
                else if (alfa2_med[i] == 0.0)
                {
                    r2_med[i] = 1.0;
                }
                else if (a2_med[i] > 0)
                {
                    if (i != 0)
                    {
                        r2_med[i] = alfa2_med[i-1]/alfa2_med[i];
                    }
                    else
                    {
                        r2_med[i] = 1.0;
                    }
                }
            }
            
            // ✅ PROTECCIÓN: Limitar valores extremos de r que causan inestabilidad
            if (!std::isfinite(r1_med[i]) || fabs(r1_med[i]) > 1e6) r1_med[i] = 1.0;
            if (!std::isfinite(r2_med[i]) || fabs(r2_med[i]) > 1e6) r2_med[i] = 1.0;

            //! Computing fi_i
            if (nGetEquationLimiterFlux() == 1)
            {
                //! MinMod - ⚡ OPTIMIZADO: max(0, min(1, r)) sin crear vectores
                fi1_med[i] = std::max(0.0, std::min(1.0, r1_med[i]));
                fi2_med[i] = std::max(0.0, std::min(1.0, r2_med[i]));
            }
            else if (nGetEquationLimiterFlux() == 2)
            {
                //! Roe's Superbee
                fi1_med[i] = dMaxVectorValue({0.0, dMinVectorValue({2*r1_med[i], 1.}), dMinVectorValue({r1_med[i], 2.})});
                fi2_med[i] = dMaxVectorValue({0.0, dMinVectorValue({2*r2_med[i], 1.}), dMinVectorValue({r2_med[i], 2.})});
            }
            else if (nGetEquationLimiterFlux() == 3)

            {
                //! Van Leer
                fi1_med[i] = (fabs(r1_med[i]) + r1_med[i])/(1+fabs(r1_med[i]));
                fi2_med[i] = (fabs(r2_med[i]) + r2_med[i])/(1+fabs(r2_med[i]));
            }
            else if (nGetEquationLimiterFlux() == 4)
            {
                fi1_med[i] = (r1_med[i]*r1_med[i] + r1_med[i])/(1+r1_med[i]*r1_med[i]);
                fi2_med[i] = (r2_med[i]*r2_med[i] + r2_med[i])/(1+r2_med[i]*r2_med[i]);
            }
        }
        for (int i = 0; i < m_nCrossSectionsNumber-1; i++)
        {
            //! Computing factor
            vFactor1[i] = alfa1_med[i]*psi1_med[i]*(1-m_dLambda*fabs(a1_med[i]))*(1-fi1_med[i]);
            vFactor2[i] = alfa2_med[i]*psi2_med[i]*(1-m_dLambda*fabs(a2_med[i]))*(1-fi2_med[i]);

            m_vCrossSectionD1Factor[i+1] = 0.5*(vFactor1[i] + vFactor2[i]);
            m_vCrossSectionD2Factor[i+1] = 0.5*(vFactor1[i]*a1_med[i] + vFactor2[i]*a2_med[i]);
        }
        m_vCrossSectionD1Factor[0] = m_vCrossSectionD1Factor[1];
        m_vCrossSectionD2Factor[0] = m_vCrossSectionD2Factor[1];
        m_vCrossSectionD1Factor[m_nCrossSectionsNumber] = m_vCrossSectionD1Factor[m_nCrossSectionsNumber-1];
        m_vCrossSectionD2Factor[m_nCrossSectionsNumber] = m_vCrossSectionD2Factor[m_nCrossSectionsNumber-1];

        //! ✅ CORRECTO: Solo promediar nodos interiores (i=1 hasta i=n-2)
        //! Los nodos de frontera (i=0 e i=n-1) ya fueron asignados por las BC
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

//======================================================================================================================
//! Smooth solution
//======================================================================================================================
void CSimulation::smoothSolution() {
    // 🎯 SUAVIZADO ADAPTATIVO SOLO EN ZONAS CON SALTOS (para visualización)
    // No afecta el cálculo físico - solo el output
    
    const int n = m_nCrossSectionsNumber;
    vector<double> vAreaSmooth(n);
    vector<double> vQSmooth(n);
    
    // Copiar valores originales
    vAreaSmooth = m_vCrossSectionArea;
    vQSmooth = m_vCrossSectionQ;
    
    // Detectar y suavizar solo zonas con gradientes fuertes
    const int boundary_zone = 3;  // Suavizar 3 nodos cerca de contornos
    
    // 1️⃣ ZONA AGUAS ARRIBA (primeros 3 nodos)
    for (int i = 1; i < std::min(boundary_zone, n-1); i++) {
        // Filtro suave [0.25, 0.50, 0.25] solo cerca del contorno
        vAreaSmooth[i] = 0.25*m_vCrossSectionArea[i-1] + 
                         0.50*m_vCrossSectionArea[i] + 
                         0.25*m_vCrossSectionArea[i+1];
        vQSmooth[i] = 0.25*m_vCrossSectionQ[i-1] + 
                      0.50*m_vCrossSectionQ[i] + 
                      0.25*m_vCrossSectionQ[i+1];
    }
    
    // 2️⃣ ZONA AGUAS ABAJO (últimos 3 nodos)
    for (int i = std::max(n-boundary_zone, 1); i < n-1; i++) {
        vAreaSmooth[i] = 0.25*m_vCrossSectionArea[i-1] + 
                         0.50*m_vCrossSectionArea[i] + 
                         0.25*m_vCrossSectionArea[i+1];
        vQSmooth[i] = 0.25*m_vCrossSectionQ[i-1] + 
                      0.50*m_vCrossSectionQ[i] + 
                      0.25*m_vCrossSectionQ[i+1];
    }
    
    // 3️⃣ ZONA INTERIOR: Detectar saltos locales > 10% y suavizar
    for (int i = boundary_zone; i < n-boundary_zone; i++) {
        // Calcular variación relativa local
        double area_avg = (m_vCrossSectionArea[i-1] + m_vCrossSectionArea[i] + m_vCrossSectionArea[i+1]) / 3.0;
        double area_var = fabs(m_vCrossSectionArea[i] - area_avg);
        
        if (area_avg > 1e-6 && area_var / area_avg > 0.10) {
            // Hay salto > 10% - suavizar MUY suave [0.10, 0.80, 0.10]
            vAreaSmooth[i] = 0.10*m_vCrossSectionArea[i-1] + 
                             0.80*m_vCrossSectionArea[i] + 
                             0.10*m_vCrossSectionArea[i+1];
        }
        
        // Igual para caudal
        double Q_avg = (fabs(m_vCrossSectionQ[i-1]) + fabs(m_vCrossSectionQ[i]) + fabs(m_vCrossSectionQ[i+1])) / 3.0;
        double Q_var = fabs(fabs(m_vCrossSectionQ[i]) - Q_avg);
        
        if (Q_avg > 1e-6 && Q_var / Q_avg > 0.10) {
            vQSmooth[i] = 0.10*m_vCrossSectionQ[i-1] + 
                          0.80*m_vCrossSectionQ[i] + 
                          0.10*m_vCrossSectionQ[i+1];
        }
    }
    
    // Aplicar valores suavizados (manteniendo contornos i=0 y i=n-1 intactos)
    for (int i = 1; i < n-1; i++) {
        m_vCrossSectionArea[i] = vAreaSmooth[i];
        m_vCrossSectionQ[i] = vQSmooth[i];
    }
    
    // 🔄 Recalcular parámetros hidráulicos para consistencia
    // calculateHydraulicParameters();
}

//======================================================================================================================
//! Get the vector data given the name
//======================================================================================================================
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
    else if (strVariableName == "beta") {
        return m_vCrossSectionBeta;
    }
    else if (strVariableName == "DhDx") {
        return m_vCrossSectionDhDx;
    }
    else if (strVariableName == "rho") {
        return m_vCrossSectionRho;
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
    else if (strVariableName == "xl") {
        return m_vCrossSectionLeftRBLocation;
    }
    else if (strVariableName == "xr") {
        return m_vCrossSectionRightRBLocation;
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

//===============================================================================================================================
//! Displays information regarding the progress of the simulation
//===============================================================================================================================
void CSimulation::AnnounceProgress() {

    // ⚡ OPTIMIZACIÓN: Solo actualizar cada 100 timesteps para reducir overhead
    static int call_count = 0;
    if (++call_count % 100 != 0 && !m_bSaveTime) return;

    // Stdout is connected to a tty, so not running as a background job
    static double sdElapsed = 0;
    static double sdToGo = 0;
    const time_t tNow = time(nullptr);

    // Calculate time elapsed and remaining
    sdElapsed = difftime(tNow, m_tSysStartLoopTime);
    sdToGo = (sdElapsed * m_dSimDuration / m_dCurrentTime) - sdElapsed;

    if (m_nLogFileDetail == 2) {
        cout << "\r    - Elapsed[Remaining] Time: " << std::fixed << setprecision(3) << setw(6) << sdElapsed <<"[" << std::fixed << setprecision(3) << setw(6) << sdToGo << "] s - Progress: " << std::fixed << setprecision(3) << setw(6) << 100 * m_dCurrentTime / m_dSimDuration  << '%' << std::flush;
        LogStream << "TIMESTEP: " << m_nTimeId << " - Current time: " <<  m_dCurrentTime << endl;

    }
    else {
        if (m_bSaveTime) {
            cout << "\r    - Elapsed[Remaining] Time: " << std::fixed << setprecision(3) << setw(6) << sdElapsed <<"[" << std::fixed << setprecision(3) << setw(6) << sdToGo << "] s - Progress: " << std::fixed << setprecision(3) << setw(6) << 100 * m_dCurrentTime / m_dSimDuration  << '%' << std::flush;
        }
    }
}

//===============================================================================================================================
//! Compute the salinity.
//===============================================================================================================================
void CSimulation::calculate_salinity()
{
    // Add salinity variable to dry bed function
    for (int i = 0; i < m_nCrossSectionsNumber; i++)
    {
        m_vCrossSectionSalinity[i] = m_vCrossSectionSalinityASt[i]/m_vCrossSectionArea[i] + m_vCrossSectionSalinity[i];
        // Bound the minimum and maximum values of salinity to 0 a 35 psu

        if (m_vCrossSectionSalinity[i] < 0) m_vCrossSectionSalinity[i] = 0;
        if (m_vCrossSectionSalinity[i] > 35) m_vCrossSectionSalinity[i] = 35;
    }

    // ================================================================================================================
    // CONDICIONES DE FRONTERA AGUAS ARRIBA (UPSTREAM) - SALINIDAD
    // ================================================================================================================
    if (nGetUpwardSalinityCondition() == 0) {
        //! ✅ CONDICIÓN LIBRE (FREE): Solo permite entrada de agua dulce
        //! Si flujo entra (Q > 0): salinidad = 0 (agua dulce)
        //! Si flujo sale (Q < 0): usar gradiente cero (mantener valor actual sin difusión hacia atrás)
        if (m_vCrossSectionQ[0] > 0) {
            // Flujo entrante desde aguas arriba → agua dulce
            m_vCrossSectionSalinity[0] = 0.0;
        }
        // Si flujo saliente, mantener el valor calculado por advección (ya está en m_vCrossSectionSalinity[0])
    }
    else if (nGetUpwardSalinityCondition() == 1) {
        //! Salinidad NULA impuesta
        m_vCrossSectionSalinity[0] = 0.0;
    }
    else if (nGetUpwardSalinityCondition() == 2) {
        //! Salinidad OCEÁNICA impuesta
        m_vCrossSectionSalinity[0] = 35.0;
    }

    // ================================================================================================================
    // CONDICIONES DE FRONTERA AGUAS ABAJO (DOWNSTREAM) - SALINIDAD
    // ================================================================================================================
    if (nGetDownwardSalinityCondition() == 0) {
        //! CONDICIÓN LIBRE: Gradiente cero en frontera aguas abajo
        // Ya calculado por advección, no hacer nada
    }
    else if (nGetDownwardSalinityCondition() == 1) {
        //! Salinidad NULA impuesta
        m_vCrossSectionSalinity[m_nCrossSectionsNumber-1] = 0.0;
    }
    else if (nGetDownwardSalinityCondition() == 2) {
        //! Salinidad OCEÁNICA impuesta
        m_vCrossSectionSalinity[m_nCrossSectionsNumber-1] = 35.0;
    }
}


//===============================================================================================================================
//! Compute the salinity gradient using the backward, centered and upward finite differences
//! Diez-Minguito et al. (2013). ⚡ OPTIMIZADO: Usar vectores pre-alocados
//===============================================================================================================================
void CSimulation::calculate_salinity_gradient()
{
    // Usar vectores pre-alocados (miembros de clase)
    auto& vKAS_dif_forward = m_vSalinity_KAS_forward;
    auto& vKAS_dif_backward = m_vSalinity_KAS_backward;
    auto& vAUS_dif = m_vSalinity_AUS_diff;

    for (int i = 1; i < m_nCrossSectionsNumber-1; i++)
    {
        vKAS_dif_forward[i] = m_vCrossSectionArea[i+1]*(m_vCrossSectionSalinity[i+1]-m_vCrossSectionSalinity[i]);
        vKAS_dif_backward[i+1] = m_vCrossSectionArea[i]*(m_vCrossSectionSalinity[i+1]-m_vCrossSectionSalinity[i]);
        vAUS_dif[i] = (m_vCrossSectionQ[i+1]*m_vCrossSectionSalinity[i+1] - m_vCrossSectionQ[i-1]*m_vCrossSectionSalinity[i-1])*m_dLambda*0.5;
    }

    // Add the end and initial values
    vKAS_dif_forward[0] = m_vCrossSectionArea[1]*(m_vCrossSectionSalinity[1]-m_vCrossSectionSalinity[0]);
    vKAS_dif_forward[m_nCrossSectionsNumber-1] = vKAS_dif_forward[m_nCrossSectionsNumber-2];

    // ✅ CORRECCIÓN: En frontera upstream con condición FREE, no permitir difusión hacia atrás
    if (nGetUpwardSalinityCondition() == 0) {
        // Condición FREE: solo difusión hacia adelante, no hacia atrás (evita acumulación artificial)
        vKAS_dif_backward[1] = m_vCrossSectionArea[0]*(m_vCrossSectionSalinity[1]-m_vCrossSectionSalinity[0]);
        vKAS_dif_backward[0] = 0.0;  // SIN difusión hacia atrás de la frontera
    } else {
        // Condiciones impuestas: difusión normal
        vKAS_dif_backward[1] = m_vCrossSectionArea[0]*(m_vCrossSectionSalinity[1]-m_vCrossSectionSalinity[0]);
        vKAS_dif_backward[0] = vKAS_dif_backward[1];
    }

    // ✅ FRONTERA AGUAS ARRIBA: Esquema upwind basado en dirección del flujo
    if (nGetUpwardSalinityCondition() == 0) {
        // Condición FREE: agua dulce cuando entra, permite salida cuando refluye
        if (m_vCrossSectionQ[0] >= 0) {
            // Flujo entrante (Q > 0): advectar agua dulce desde frontera
            vAUS_dif[0] = m_vCrossSectionQ[0] * 0.0 * m_dLambda; // S[0]=0 para agua dulce entrante
        } else {
            // Flujo saliente (Q < 0): advectar la salinidad del propio nodo hacia fuera
            // CORRECCIÓN: Usar S[0] no S[1] para permitir lavado con la marea
            vAUS_dif[0] = m_vCrossSectionQ[0] * m_vCrossSectionSalinity[0] * m_dLambda;
        }
    } else {
        // Condiciones impuestas (tipo 1 o 2): advección normal
        vAUS_dif[0] = (m_vCrossSectionQ[1]*m_vCrossSectionSalinity[1] - m_vCrossSectionQ[0]*m_vCrossSectionSalinity[0])*m_dLambda;
    }
    
    // FRONTERA AGUAS ABAJO: advección normal
    vAUS_dif[m_nCrossSectionsNumber-1] = (m_vCrossSectionQ[m_nCrossSectionsNumber-1]*m_vCrossSectionSalinity[m_nCrossSectionsNumber-1] - m_vCrossSectionQ[m_nCrossSectionsNumber-2]*m_vCrossSectionSalinity[m_nCrossSectionsNumber-2])*m_dLambda;

    //! Calculate de ASt term (temporal gradient)
    for (int i = 0; i < m_nCrossSectionsNumber; i++)
    {
        m_vCrossSectionSalinityASt[i] = m_dLongitudinalDispersion*m_dLambda*m_dLambda/m_dTimestep*(vKAS_dif_forward[i] - vKAS_dif_backward[i]) - vAUS_dif[i];
    }
}


//===============================================================================================================================
//! Compute the bedload and suspended sediment transport using the van Rijn equation(van Rijn, 1992)
//===============================================================================================================================
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


//===============================================================================================================================
//! Compute the density given the salinity and sediment concentration
//===============================================================================================================================
void CSimulation::calculate_density()
{   
    for (int i = 0; i < m_nCrossSectionsNumber; i++)
    {
        double rhow = 1000 * (1 + dGetBetaSalinityConstant() * m_vCrossSectionSalinity[i]);
        if (m_bDoSedimentTransport) {
            m_vCrossSectionRho[i] = rhow + (1 - rhow/1000/m_vCrossSectionRhos[i])*m_vCrossSectionQt[i]/(m_vCrossSectionArea[i]*m_vCrossSectionDX[i])*m_dTimestep;
        }
        else {
            m_vCrossSectionRho[i] = rhow;
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


//======================================================================================================================
//! Precompute and cache estuary geometric and hydraulic data for efficient runtime access
//! Extracts data from estuary[] objects into flat vectors to avoid repeated virtual calls
//======================================================================================================================
void CSimulation::precomputeEstuaryData() {
    
    // Pre-calcular datos escalares básicos
    m_vElevationSectionsCount.resize(m_nCrossSectionsNumber);
    m_vBedZ.resize(m_nCrossSectionsNumber);
    m_vManningN.resize(m_nCrossSectionsNumber);
    m_vPositionX.resize(m_nCrossSectionsNumber);
    
    // ✅ MEJOR: Inicializar con dimensiones conocidas
    m_vWidth.assign(m_nCrossSectionsNumber, vector<double>());        // Vector vacío inicial
    m_vBeta.assign(m_nCrossSectionsNumber, vector<double>());         
    m_vLeftY.assign(m_nCrossSectionsNumber, vector<double>());        
    m_vRightY.assign(m_nCrossSectionsNumber, vector<double>());       
    m_vEstuaryAreas.assign(m_nCrossSectionsNumber, vector<double>());
    m_vEstuaryHydraulicRadius.assign(m_nCrossSectionsNumber, vector<double>());
    m_vEstuaryWaterDepths.assign(m_nCrossSectionsNumber, vector<double>());
    m_vEstuaryI1.assign(m_nCrossSectionsNumber, vector<double>());
    m_vPrecalculatedSecondTerm.assign(m_nCrossSectionsNumber, vector<double>());

    
    for (int i = 0; i < m_nCrossSectionsNumber; i++) {
        // Datos escalares básicos
        m_vBedZ[i] = estuary[i].dGetZ();
        m_vManningN[i] = estuary[i].dGetManningNumber();
        m_vPositionX[i] = estuary[i].dGetX();
        m_vElevationSectionsCount[i] = estuary[i].nGetElevationSectionsNumber();
        
        // Vectores de datos hidráulicos
        m_vWidth[i] = estuary[i].vGetWidth();   
        m_vEstuaryAreas[i] = estuary[i].vGetArea();
        m_vEstuaryHydraulicRadius[i] = estuary[i].vGetHydraulicRadius();
        m_vEstuaryWaterDepths[i] = estuary[i].vGetWaterDepth();
        m_vBeta[i] = estuary[i].vGetBeta();
        m_vLeftY[i] = estuary[i].vGetLeftRBLocation();
        m_vRightY[i] = estuary[i].vGetRightRBLocation();
        m_vEstuaryI1[i] = estuary[i].vGetI1();

        const auto& areas = m_vEstuaryAreas[i];
        const auto& hydraulicRadius = m_vEstuaryHydraulicRadius[i];
        
        m_vPrecalculatedSecondTerm[i].resize(areas.size());
        for (size_t j = 0; j < areas.size(); j++) {
            m_vPrecalculatedSecondTerm[i][j] = areas[j] * pow(hydraulicRadius[j], 2.0/3.0);
        }
    }
    }

//======================================================================================================================
//! GENERAR NOMBRE DE ARCHIVO AUTOMÁTICO BASADO EN PARÁMETROS
//======================================================================================================================
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