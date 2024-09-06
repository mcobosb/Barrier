#include <iostream>
#include <main.h>
#include <utils.h>
#include <windows.h>
using std::cerr;
using std::cout;
using std::endl;
using std::ios;

#include <iomanip>
using std::setprecision;

#include <string>
using std::to_string;

#include <cmath>
using std::sqrt;
using std::fabs;
using std::pow;

// #include "main.h"
#include "simulation.h"
#include "screen_presenter.h"
#include "hydrograph.h"
#include "cross_section.h"
#include "utils.h"

//===============================================================================================================================
//! The CSimulation constructor
//===============================================================================================================================
CSimulation::CSimulation ()
{
    m_dSimDuration =
    m_dSimTimestep = 0.0;

    m_nInitialEstuarineCondition = 0;

    vector<string> m_vOutputVariables;

    vector<CCrossSection> estuary;
    vector<CHydrograph> hydrographs;
}

//===============================================================================================================================
//! The CSimulation destructor
//===============================================================================================================================
CSimulation::~CSimulation () = default;


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

//! Method for getting the computation of water density
bool CSimulation::bGetDoWaterDensity() const {
    return m_bDoWaterDensity;
}

//! Method for setting the computation of water density
void CSimulation::bSetDoWaterDensity(bool computeWaterDensity) {
    m_bDoWaterDensity = computeWaterDensity;
};

//! Method for getting the along channel constant elevation
double CSimulation::dGetBetaSalinityConstant() const {
    return m_dBetaSalinityConstant;
}

//! Method for setting the along channel constant elevation
void CSimulation::dSetBetaSalinityConstant(double betaSalinity) {
    m_dBetaSalinityConstant = betaSalinity;
};

//! Method for getting the along channel constant elevation
double CSimulation::dGetLongitudinalDispersionConstant() const {
    return m_dLongitudinalDispersion;
}

//! Method for setting the along channel constant elevation
void CSimulation::dSetLongitudinalDispersionConstant(double longitudinalDispersion) {
    m_dLongitudinalDispersion = longitudinalDispersion;
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
    return m_bDoMackComarckLimiterFlux;
}
//! Method for setting if McComarck limiter flux is applied
void CSimulation::bSetDoMcComarckLimiterFlux(bool doMcComarckLimiterFlux) {
    m_bDoMackComarckLimiterFlux = doMcComarckLimiterFlux;
}

//! Method for getting equation limiter flux
int CSimulation::nGetEquationLimiterFlux() const {
    return m_nEquationMacComarckLimiterFlux;
}
//! Method for setting equation limiter flux
void CSimulation::nSetEquationLimiterFlux(int equationLimiterFlux) {
    m_nEquationMacComarckLimiterFlux = equationLimiterFlux;
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
bool CSimulation::bDoSimulation(int nArg, char const* pcArgv[]){

    initializeVectors();
    calculateBedSlope();
    calculateAlongEstuaryInitialConditions();
    calculateIs();

    //! Create the NetCDF file with dimensions, coordinates and variables
    int nRtn = writer.nDefineNetCDFFile(this);

    //! Check if any part of the estuary is dry
    if (bGetDoDryBed())
        dryArea();

    m_nTimeId = 0;
    //! TODO 013: Hydrographs along the estuary have to be reduced by dX of the input location
    while (m_dCurrentTime < m_dSimDuration)
    {
        cout << "Simulation time: " << m_dCurrentTime << endl;
        calculateBoundaryConditions();
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
        m_nPredictor = 1;
        //! Calculate the friction slope
        //! TODO 022: use the value computed in the following function
        calculateFrictionSlope();

        calculateGSAterms();
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
        calculateGSAterms();
        calculateFlowTerms();
        calculateSourceTerms();

        //! Change the values of the source terms if any part of the estuary is dry
        if (bGetDoDryBed())
            dryTerms();

        calculateCorrector();

        //! Boundary conditions
        updateCorrectorBoundaries();

        //! Check if any part of the estuary is dry (corrector)
        if (bGetDoDryBed())
            dryArea();

        //! Save corrected variables
        // if (m_bSaveTime) {
        //
        // }

        //==============================================================================================================
        //! STEP 03: update the following timestep
        //==============================================================================================================
        mergePredictorCorrector();

        smoothSolution();

        //! Compute the salinity gradient
        // if (bGetDoWaterDensity()) {
        //     //! TODO 016: include salinity
        // }

        if (m_bSaveTime) {
            nRtn = writer.nSetOutputData(this);
        }

        //! Increase counter
        m_dCurrentTime = m_dCurrentTime + m_dTimestep;
        if (m_bSaveTime) {
            m_nTimeId++;
        }
    }

    return nRtn;
}


//======================================================================================================================
//! Initialize vectorial variables
//======================================================================================================================
void CSimulation::initializeVectors()
{
    m_nPredictor = -1;

    const double nCrossSectionsNumber = this->m_nCrossSectionsNumber;
    const vector<double> vZeros(static_cast<size_t>(nCrossSectionsNumber), 0.0);

    m_vCrossSectionArea =
    m_vPredictedCrossSectionArea =
    m_vCorrectedCrossSectionArea =
    m_vCrossSectionQ =
    m_vPredictedCrossSectionQ =
    m_vCorrectedCrossSectionQ =
    m_vCrossSectionHydraulicRadius =
    m_vCrossSectionI1 =
    m_vCrossSectionI2 =
    m_vCrossSectionWidth =
    m_vCrossSectionBeta =
    m_vCrossSectionU =
    m_vCrossSectionC =
    m_vCrossSectionSalinity =
    m_vCrossSectionQb =
    m_vCrossSectionQs =
    m_vCrossSectionQt =
    m_vCrossSectionRho =
    m_vCrossSectionLeftRBLocation =
    m_vCrossSectionRightRBLocation =
    m_vCrossSectionBedSlope =
    m_vCrossSectionFrictionSlope =
    m_vCrossSectionDX =
    m_vCrossSectionManningNumber =
    m_vCrossSectiongAS0 =
    m_vCrossSectiongASf =
    m_vCrossSectionF0 =
    m_vCrossSectionF1 =
    m_vCrossSectionGv0 =
    m_vCrossSectionGv1 = vZeros;
    // m_vCrossSectionF0 =
    // m_vCrossSectionF1 =  vZeros;

    if (this->nGetDownwardEstuarineCondition() == 0) {
        m_vCrossSectionQ =
        m_vCrossSectionArea = vZeros;
    }
    else if (this->nGetDownwardEstuarineCondition() == 1)
    {
        m_vCrossSectionElevation = vZeros;
    }
    else
    {
        m_vCrossSectionQ = vZeros;
    }

    int nTimestepsNumber = static_cast<int> (m_dSimDuration/m_dSimTimestep);
    m_vOutputTimesIds = vector<int>(nTimestepsNumber, 0);

    for (int i = 0; i < nTimestepsNumber ; i++) {
        m_vOutputTimes.push_back(static_cast<double>(i)*m_dSimTimestep);
    }
}

//======================================================================================================================
//! Calculate the initial bed slope
//======================================================================================================================
void CSimulation::calculateBedSlope() {
    double dX1 = 0.0;
    double dX2 = 0.0;

    //Calculate S0
    for (int i = 0; i < m_nCrossSectionsNumber ; i++) {
        if (i == 0) {
            // Compute dx as forward difference
            dX1 = this->estuary[1].dGetX() - this->estuary[0].dGetX();
            //! TODO 009: check the sign of this->estuary[i+1].dGetZ() - this->estuary[i].dGetZ()
            this->m_vCrossSectionBedSlope[i] = (this->estuary[0].dGetZ() - this->estuary[1].dGetZ())/dX1;
            //! Save dX into a vector
            m_vCrossSectionDX[i] = dX1;
        }
        else if (i == m_nCrossSectionsNumber - 1){
            // Compute dx as backward difference
            dX2 = this->estuary[i].dGetX() - this->estuary[i-1].dGetX();
            this->m_vCrossSectionBedSlope[i] = (this->estuary[i-1].dGetZ() - this->estuary[i].dGetZ())/dX2;
            //! Save dX into a vector (last node backward)
            m_vCrossSectionDX[i] = dX2;
        }
        else {
            dX1 = this->estuary[i+1].dGetX() - this->estuary[i].dGetX();
            dX2 = this->estuary[i].dGetX() - this->estuary[i-1].dGetX();
            this->m_vCrossSectionBedSlope[i] = (this->estuary[i-1].dGetZ() - this->estuary[i+1].dGetZ())/(dX1 + dX2);
            //! Save dX into a vector
            m_vCrossSectionDX[i] = dX1;
        }

    }
}

//======================================================================================================================
//! Calculate the initial conditions along the estuary
//======================================================================================================================
void CSimulation::calculateAlongEstuaryInitialConditions() {

    if (m_nInitialEstuarineCondition == 1) {
        for (int i = 0; i < m_nCrossSectionsNumber; i++) {
            const double dManningFactor = this->m_vCrossSectionQ[i]*this->estuary[i].dGetManningNumber()/(sqrt(this->m_vCrossSectionBedSlope[i]));

            vector<double> vCrossSectionAreaTmp = this->estuary[i].vGetArea();
            vector<double> vCrossSectionHydraulicRadiusTmp = this->estuary[i].vGetHydraulicRadius();
            vector<double> dSecondTerm;

            //! Second term to obtain the area from slope equation in open channels
            for (size_t j = 0; j < vCrossSectionAreaTmp.size(); j++) {
                dSecondTerm.push_back(vCrossSectionAreaTmp[j]*pow(vCrossSectionHydraulicRadiusTmp[j], 2.0/3.0));
            }
            this->m_vCrossSectionArea[i] = linearInterpolation1d(dManningFactor, dSecondTerm, vCrossSectionAreaTmp);
        }
    }
    else if (m_nInitialEstuarineCondition == 2) {
        for (int i = 0; i < m_nCrossSectionsNumber; i++) {
            // this->m_vCrossSectionArea.push_back(this->m_vCrossSectionElevation[i]);
            vector<double> vCrossSectionAreaTmp = this->estuary[i].vGetArea();
            vector<double> vCrossSectionElevationTmp = this->estuary[i].vGetElevation();
            vector<double> vCrossSectionHydraulicRadiusTmp = this->estuary[i].vGetHydraulicRadius();

            //! m_vCrossSectionElevation is the elevation from the water depth m_dZ of every cross-section
            m_vCrossSectionArea[i] = linearInterpolation1d(this->m_vCrossSectionElevation[i],vCrossSectionElevationTmp, vCrossSectionAreaTmp);
            m_vCrossSectionHydraulicRadius[i] = linearInterpolation1d(this->m_vCrossSectionArea[i], vCrossSectionAreaTmp, vCrossSectionHydraulicRadiusTmp);

            //! Compute Q given the area
            m_vCrossSectionQ[i] = m_vCrossSectionArea[i]*sqrt(m_vCrossSectionBedSlope[i])*pow(m_vCrossSectionHydraulicRadius[i], 2.0/3.0)/this->estuary[i].dGetManningNumber();
        }
    }
    else {
        //! Water level in calm
        for (int i = 0; i < m_nCrossSectionsNumber; i++) {
            vector<double> vCrossSectionAreaTmp = this->estuary[i].vGetArea();
            vector<double> vCrossSectionElevationTmp = this->estuary[i].vGetElevation();

            m_vCrossSectionElevation.push_back(-estuary[i].dGetZ());
            m_vCrossSectionArea[i] = linearInterpolation1d(m_vCrossSectionElevation[i],vCrossSectionElevationTmp, vCrossSectionAreaTmp);
        }
    }

    for (int i = 0; i < m_nCrossSectionsNumber; i++) {
        //! Insert the Manning number onto the simulation object
        m_vCrossSectionManningNumber[i] = this->estuary[i].dGetManningNumber();
    }
}


//======================================================================================================================
//! Interpolate function
//======================================================================================================================
double CSimulation::linearInterpolation1d(const double dValue, const vector<double> &vX, const vector<double> &vY) {
    //! Encuentra el intervalo adecuado [x[i], x[i+1]]
    for (size_t i = 0; i < vX.size() - 1; i++) {
        if (dValue >= vX[i] && dValue <= vX[i+1]) {
            // Do the linear interpolation
            double slope = (vY[i+1] - vY[i]) / (vX[i+1] - vX[i]);
            return vY[i] + slope * (dValue - vX[i]);
        }
    }
    //! TODO 017: Error if dValue is outside the range ...
    return 0.0;
}
//======================================================================================================================
//! Obtain the hydraulic parameters from the cross-sectional data
//======================================================================================================================
void CSimulation::calculateHydraulicParameters() {
    //! Number of estuarine cross-sections
    const int nCrossSections = m_nCrossSectionsNumber;
    vector<int> vElevationSection;

    for (int i = 0; i < nCrossSections; i++) {
        //! TODO 008: Create getter and setter for cross-section area, hydraulic radius, elevation and slope
        double dArea = this->m_vCrossSectionArea[i];

        if (estuary[i].dGetArea(0) > dArea) {
            //! Take the first node if dArea < the Area of the first elevation node
            vElevationSection.push_back(0);
            getFirstHydraulicParameters(i);
        }
        else if (estuary[i].dGetArea(nCrossSections -1) < dArea) {
            //! Take the last node if dArea > the Area of the last elevation node
            vElevationSection.push_back(nCrossSections - 1);
            getLastHydraulicParameters(i);
        }
        else
        {
            for (int j = 0; j < estuary[i].nGetElevationSectionsNumber()-1; j++) {
                //! Getting the elevations node which Area is below dArea and next node Area higher than dArea
                if ((estuary[i].dGetArea(j) < dArea) && (estuary[i].dGetArea(j+1) > dArea)) {
                    vElevationSection.push_back(j);
                    interpolateHydraulicParameters(dArea, i, j);
                    break;
                }
            }
        }
    }
}

//======================================================================================================================
//! Interpolate the Hydraulic parameters between the elevation node given.
//======================================================================================================================
void CSimulation::interpolateHydraulicParameters(double dArea, int nCrossSection, int nElevationNode) {
    double dInterpolationFactor = (dArea - this->estuary[nCrossSection].dGetArea(nElevationNode)) / (this->estuary[nCrossSection].dGetArea(nElevationNode+1) - this->estuary[nCrossSection].dGetArea(nElevationNode));
    m_vCrossSectionHydraulicRadius[nCrossSection] = dInterpolationFactor*(this->estuary[nCrossSection].dGetHydraulicRadius(nElevationNode+1) - this->estuary[nCrossSection].dGetHydraulicRadius(nElevationNode)) + this->estuary[nCrossSection].dGetHydraulicRadius(nElevationNode);
    m_vCrossSectionWidth[nCrossSection] =  dInterpolationFactor*(this->estuary[nCrossSection].dGetWidth(nElevationNode+1) - this->estuary[nCrossSection].dGetWidth(nElevationNode)) + this->estuary[nCrossSection].dGetWidth(nElevationNode);
    m_vCrossSectionElevation[nCrossSection] = dInterpolationFactor*(this->estuary[nCrossSection].dGetElevation(nElevationNode+1) - this->estuary[nCrossSection].dGetElevation(nElevationNode)) + this->estuary[nCrossSection].dGetElevation(nElevationNode);
    m_vCrossSectionBeta[nCrossSection] = dInterpolationFactor*(this->estuary[nCrossSection].dGetBeta(nElevationNode+1) - this->estuary[nCrossSection].dGetBeta(nElevationNode)) + this->estuary[nCrossSection].dGetBeta(nElevationNode);
    m_vCrossSectionI1[nCrossSection] = dInterpolationFactor*(this->estuary[nCrossSection].dGetI1(nElevationNode+1) - this->estuary[nCrossSection].dGetI1(nElevationNode)) + this->estuary[nCrossSection].dGetI1(nElevationNode);
    m_vCrossSectionI2[nCrossSection] = dInterpolationFactor*(this->estuary[nCrossSection].dGetI2(nElevationNode+1) - this->estuary[nCrossSection].dGetI2(nElevationNode)) + this->estuary[nCrossSection].dGetI2(nElevationNode);
    m_vCrossSectionLeftRBLocation[nCrossSection] = dInterpolationFactor*(this->estuary[nCrossSection].dGetLeftY(nElevationNode+1) - this->estuary[nCrossSection].dGetLeftY(nElevationNode)) + this->estuary[nCrossSection].dGetLeftY(nElevationNode);
    m_vCrossSectionRightRBLocation[nCrossSection] = dInterpolationFactor*(this->estuary[nCrossSection].dGetRightY(nElevationNode+1) - this->estuary[nCrossSection].dGetRightY(nElevationNode)) + this->estuary[nCrossSection].dGetRightY(nElevationNode);
}


//======================================================================================================================
//! Get first node hydraulic parameters
//======================================================================================================================
void CSimulation::getFirstHydraulicParameters(const int nCrossSection) {
    m_vCrossSectionHydraulicRadius[nCrossSection] = this->estuary[nCrossSection].dGetHydraulicRadius(0);
    m_vCrossSectionWidth[nCrossSection] =  this->estuary[nCrossSection].dGetWidth(0);
    m_vCrossSectionElevation[nCrossSection] = this->estuary[nCrossSection].dGetElevation(0);
    m_vCrossSectionBeta[nCrossSection] = this->estuary[nCrossSection].dGetBeta(0);
    m_vCrossSectionI1[nCrossSection] = this->estuary[nCrossSection].dGetI1(0);
    m_vCrossSectionI2[nCrossSection] = this->estuary[nCrossSection].dGetI2(0);
    m_vCrossSectionLeftRBLocation[nCrossSection] = this->estuary[nCrossSection].dGetLeftY(0);
    m_vCrossSectionRightRBLocation[nCrossSection] = this->estuary[nCrossSection].dGetRightY(0);
}


//======================================================================================================================
//! Get last node hydraulic parameters
//======================================================================================================================
void CSimulation::getLastHydraulicParameters(const int nCrossSection) {
    const int nLastNode = this->estuary[nCrossSection].nGetElevationSectionsNumber() - 1;
     m_vCrossSectionHydraulicRadius[nCrossSection] = this->estuary[nCrossSection].dGetHydraulicRadius(nLastNode);
     m_vCrossSectionWidth[nCrossSection] =  this->estuary[nCrossSection].dGetWidth(nLastNode);
     m_vCrossSectionElevation[nCrossSection] = this->estuary[nCrossSection].dGetElevation(nLastNode);
     m_vCrossSectionBeta[nCrossSection] = this->estuary[nCrossSection].dGetBeta(nLastNode);
     m_vCrossSectionI1[nCrossSection] = this->estuary[nCrossSection].dGetI1(nLastNode);
     m_vCrossSectionI2[nCrossSection] = this->estuary[nCrossSection].dGetI2(nLastNode);
     m_vCrossSectionLeftRBLocation[nCrossSection] = this->estuary[nCrossSection].dGetLeftY(nLastNode);
     m_vCrossSectionRightRBLocation[nCrossSection] = this->estuary[nCrossSection].dGetRightY(nLastNode);
}


//======================================================================================================================
//! Compute the Courant number
//======================================================================================================================
void CSimulation::calculateTimestep() {
    double dTimestepTmp = 0.0;
    double dMinTimestep = 10000000.0;
    double dWaterDensityFactor = 1.0;
    double dX = 0.0;

    for (int i=0; i< this->m_nCrossSectionsNumber-1; i++) {
        if (m_vCrossSectionArea[i] != DRY_AREA) {
            //! TODO 016: if all cross-sections are dry?
            dX = (this->estuary[i+1].dGetX() - this->estuary[i].dGetX());
            //! Mean water flow
            m_vCrossSectionU[i] = m_vCrossSectionQ[i]/m_vCrossSectionArea[i];

            //! Perturbation celerity
            m_vCrossSectionC[i] = sqrt(G*m_vCrossSectionArea[i]/m_vCrossSectionWidth[i]);

            //! TODO 003: Courant number due to perturbation on density field
            if (m_bDoWaterDensity) {
                dWaterDensityFactor = m_dCourantNumber*dX/(m_vCrossSectionWidth[i]*m_dBetaSalinityConstant);
            }

            //! Compute the timestep given the Courant number
            dTimestepTmp = m_dCourantNumber*dX/(fabs(m_vCrossSectionU[i] + m_vCrossSectionC[i]));


            if (dTimestepTmp < dMinTimestep) {
                dMinTimestep = dTimestepTmp;
            }

            if (m_bDoWaterDensity && (dWaterDensityFactor < dMinTimestep)) {
                dMinTimestep = dWaterDensityFactor;
            }
        }
    }
    m_dTimestep = dMinTimestep;

    //! Check if m_dTimestep is lower than 1 seconds
    if (m_dTimestep < 1) {
        m_dTimestep = 1;
    }
    //! Check if m_dTimestep is higher than save timestep
    if (m_dTimestep > m_dSimTimestep) {
        m_dTimestep = m_dSimTimestep;
    }

    //! Check if the timestep must be reduced to the following save step
    if ((m_nTimeId < static_cast<int>(m_vOutputTimes.size()-1)) && (m_dCurrentTime + m_dTimestep > m_vOutputTimes[m_nTimeId+1])) {

        // double dNextTimestep = m_vOutputTimes[m_nTimeId+1];
        // if (dNextTimestep - (m_dTimestep + m_dCurrentTime) < 0.0) {
        m_dTimestep = m_vOutputTimes[m_nTimeId+1] - m_dCurrentTime;
        m_bSaveTime = true;
        // }
    }
    else {
            m_bSaveTime = false;
    }
    m_dLambda = m_dTimestep / dMinVectorValue(m_vCrossSectionDX);
}

//===============================================================================================================================
//! Compute the source terms I1 (difference between pressure thrusts applied over the frontiers x1 and x2) and
//! I2(pressure forces due the channel width variation)
//===============================================================================================================================
void CSimulation::calculateIs()
{


    //Calculate I1
    for (int i = 0; i < m_nCrossSectionsNumber; i++) {
        double dValue = 0.0;
        this->estuary[i].dAppend2Vector("I1", dValue);
        for (int j = 1; j < estuary[i].nGetElevationSectionsNumber(); j++) {
            dValue = (this->estuary[i].dGetSigma(j) + this->estuary[i].dGetSigma(j-1)) / 2.0*this->estuary[i].dGetElevation(j)*(this->estuary[i].dGetElevation(j) - this->estuary[i].dGetElevation(j-1));
            this->estuary[i].dAppend2Vector("I1", dValue);
        }
    }


    double dI1dx = 0.0;
    double dX1 = 0.0;
    double dX2 = 0.0;

    //Calculate I2
    for (int i = 0; i < m_nCrossSectionsNumber ; i++) {

        if (i == 0) {
            // Compute dx as forward difference
            dX1 = this->estuary[i+1].dGetX() - this->estuary[i].dGetX();
        }
        else if (i == m_nCrossSectionsNumber - 1){
            // Compute dx as backward difference
            dX2 = this->estuary[i].dGetX() - this->estuary[i-1].dGetX();
        }
        else {
            dX1 = this->estuary[i+1].dGetX() - this->estuary[i].dGetX();
            dX2 = this->estuary[i].dGetX() - this->estuary[i-1].dGetX();
        }

        for (int j = 0; j < estuary[i].nGetElevationSectionsNumber(); j++) {
            double dhdx = 0.0;
            if (i == 0) {
                // Upward finite-difference
                dhdx = (this->estuary[i+1].dGetElevation(j) - this->estuary[i].dGetElevation(j))/dX1;
                dI1dx = (this->estuary[i+1].dGetI1(j) - this->estuary[i].dGetI1(j))/dX1;
            }
            else if (i == m_nCrossSectionsNumber - 1) {
                // Downward finite-difference
                dhdx = (this->estuary[i].dGetElevation(j) - this->estuary[i-1].dGetElevation(j))/dX2;
                dI1dx = (this->estuary[i].dGetI1(j) - this->estuary[i-1].dGetI1(j))/dX2;
            }
            else {
                // Central finite-difference
                dhdx = (this->estuary[i+1].dGetElevation(j) - this->estuary[i-1].dGetElevation(j))/(dX1 +dX2);
                dI1dx = (this->estuary[i+1].dGetI1(j) - this->estuary[i-1].dGetI1(j))/(dX1 + dX2);

            }
            double dValue = dI1dx - this->estuary[i].dGetArea(j)*dhdx;
            this->estuary[i].dAppend2Vector("I2", dValue);
        }
    }
}


//===============================================================================================================================
//! Calculate the boundary conditions at time t
//===============================================================================================================================
void CSimulation::calculateBoundaryConditions() {
    //! Calculate upward boundary condition at time t
    if (this->nGetUpwardEstuarineCondition() == 1) {
        //! Water flow given upstream, punctual water flow means distributed between sections
        m_dUpwardBoundaryValue = linearInterpolation1d(m_dTimestep, m_vUpwardBoundaryConditionTime, m_vUpwardBoundaryConditionValue);
        m_dNextUpwardBoundaryValue = m_dUpwardBoundaryValue;
    }
    else if (this->nGetUpwardEstuarineCondition() == 2) {
        //! Elevation given upstream
        double dUpwardBoundaryValue = linearInterpolation1d(m_dTimestep, m_vUpwardBoundaryConditionTime, m_vUpwardBoundaryConditionValue);
        m_dUpwardBoundaryValue = linearInterpolation1d(-estuary[0].dGetZ()  + dUpwardBoundaryValue, estuary[0].vGetElevation(), estuary[0].vGetArea());
        m_dNextUpwardBoundaryValue = linearInterpolation1d(-estuary[1].dGetZ()  + dUpwardBoundaryValue, estuary[1].vGetElevation(), estuary[1].vGetArea());
    }
    else {};

    //! Calculate downward boundary condition at time t
    if (this->nGetDownwardEstuarineCondition() == 1) {
        //! Water flow given upstream
        m_dDownwardBoundaryValue = linearInterpolation1d(m_dTimestep, m_vDownwardBoundaryConditionTime, m_vDownwardBoundaryConditionValue);
        m_dNextDownwardBoundaryValue = m_dDownwardBoundaryValue;
    }
    else if (this->nGetDownwardEstuarineCondition() == 2) {
        double dDownwardBoundaryValue = linearInterpolation1d(m_dTimestep, m_vDownwardBoundaryConditionTime, m_vDownwardBoundaryConditionValue);
        m_dDownwardBoundaryValue = linearInterpolation1d(-estuary[m_nCrossSectionsNumber-1].dGetZ() + dDownwardBoundaryValue, estuary[m_nCrossSectionsNumber-1].vGetElevation(), estuary[m_nCrossSectionsNumber-1].vGetArea());
        m_dNextDownwardBoundaryValue = linearInterpolation1d(-estuary[m_nCrossSectionsNumber-2].dGetZ() + dDownwardBoundaryValue, estuary[m_nCrossSectionsNumber-2].vGetElevation(), estuary[m_nCrossSectionsNumber-2].vGetArea());
    }
    else {};
}

//===============================================================================================================================
//! Change Area and Water flow if soil is dry
//===============================================================================================================================
void CSimulation::dryArea()
{
    if (m_nPredictor == 1) {
        for (int i = 0; i < this->m_nCrossSectionsNumber; i++)
        {
            //! Define the dry area as a water elevation of 1 cm
            if (m_vPredictedCrossSectionArea[i] <= DRY_AREA) {
                m_vPredictedCrossSectionArea[i] = DRY_AREA;
                m_vPredictedCrossSectionQ[i] = DRY_Q;
            }
        }
    }
    else if (m_nPredictor == 0) {
        for (int i = 0; i < this->m_nCrossSectionsNumber; i++)
        {
            //! Define the dry area as a water elevation of 1 cm
            if (m_vCorrectedCrossSectionArea[i] <= DRY_AREA) {
                m_vCorrectedCrossSectionArea[i] = DRY_AREA;
                m_vCorrectedCrossSectionQ[i] = DRY_Q;
            }
        }
    }
    else {
        for (int i = 0; i < this->m_nCrossSectionsNumber; i++)
        {
            //! Define the dry area as a water elevation of 1 cm
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
        for (int i = 0; i < this->m_nCrossSectionsNumber; i++)
        {
            //! Define the dry area as a water elevation of 1 cm
            if (m_vPredictedCrossSectionArea[i] <= DRY_AREA) {
                m_vCrossSectionF0[i] = 0.0;
                m_vCrossSectionF1[i] = 0.0;
                m_vCrossSectionGv1[i] = 0.0;
            }
        }
    }
    else {
        for (int i = 0; i < this->m_nCrossSectionsNumber; i++)
        {
            //! Define the dry area as a water elevation of 1 cm
            if (m_vCorrectedCrossSectionArea[i] <= DRY_AREA) {
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
    double dValue = 0.0;
    for (int i = 0; i < this->m_nCrossSectionsNumber; i++)
    {
        dValue = CDX_MURILLO*sqrt(2*pow(m_vCrossSectionArea[i], 2.0/3.0)/(G*m_vCrossSectionDX[i]));

        if (dValue < estuary[i].dGetManningNumber()) {
            const double murilloFactor = dValue/m_vCrossSectionManningNumber[i];
            m_vCrossSectionManningNumber[i] = dValue;
            m_vCrossSectionI1[i] = m_vCrossSectionI1[i]*murilloFactor;
        }
    }
}


//======================================================================================================================
//! Calculate the initial bed slope
//======================================================================================================================
void CSimulation::calculateFrictionSlope() {

    //Calculate Sf
    for (int i = 0; i < m_nCrossSectionsNumber ; i++) {
        m_vCrossSectionFrictionSlope[i] = m_vCrossSectionQ[i]*fabs(m_vCrossSectionQ[i])*pow(m_vCrossSectionManningNumber[i], 2.0)/(pow(m_vCrossSectionArea[i],2.0)*pow(m_vCrossSectionHydraulicRadius[i], 4.0/3.0));
    }
}

//======================================================================================================================
//! Calculate gAS terms
//======================================================================================================================
void CSimulation::calculateGSAterms() {

    //! Check if do surface term balance
    if (bGetDoSurfaceTermBalance()) {
        double dMeanHydraulicRadius = 0.0;
        double dMeanQ = 0.0;
        double dMeanElevation = 0.0;
        double dMeanArea = 0.0;

        //! Average between sections
        for (int i = 0; i < m_nCrossSectionsNumber; i++) {
            if (i < m_nCrossSectionsNumber - 1) {
                dMeanArea = (m_vCrossSectionArea[i] + m_vCrossSectionArea[i+1])/2.0;
                dMeanElevation = (m_vCrossSectionElevation[i] + m_vCrossSectionElevation[i+1])/2.0;
                dMeanQ = (m_vCrossSectionQ[i] + m_vCrossSectionQ[i+1])/2.0;
                dMeanHydraulicRadius = (m_vCrossSectionHydraulicRadius[i] + m_vCrossSectionHydraulicRadius[i+1])/2.0;
            }
            else {
                dMeanArea = (m_vCrossSectionArea[i] + m_vCrossSectionArea[i-1])/2.0;
                dMeanElevation = (m_vCrossSectionElevation[i] + m_vCrossSectionElevation[i-1])/2.0;
                dMeanQ = (m_vCrossSectionQ[i] + m_vCrossSectionQ[i-1])/2.0;
                dMeanHydraulicRadius = (m_vCrossSectionHydraulicRadius[i] + m_vCrossSectionHydraulicRadius[i-1])/2.0;
            }
            m_vCrossSectiongAS0[i] = G*dMeanArea*dMeanElevation;
            m_vCrossSectiongASf[i] = G*pow(m_vCrossSectionManningNumber[i], 2.0)*pow(dMeanQ, 2.0)/(dMeanArea*pow(dMeanHydraulicRadius, 4.0/3.0));
        }
    }
    else {
        for (int i = 0; i < m_nCrossSectionsNumber; i++) {
            m_vCrossSectiongAS0[i] = G*m_vCrossSectionArea[i]*m_vCrossSectionArea[i];
            if (m_nPredictor == 1) {
                m_vCrossSectiongASf[i] = G*pow(m_vCrossSectionManningNumber[i], 2.0)*pow(m_vCrossSectionQ[i], 2.0)/(m_vCrossSectionArea[i]*pow(m_vCrossSectionHydraulicRadius[i], 4.0/3.0));
            }
            else {
                m_vCrossSectiongASf[i] = G*pow(m_vCrossSectionManningNumber[i], 2.0)*pow(m_vPredictedCrossSectionQ[i], 2.0)/(m_vPredictedCrossSectionArea[i]*pow(m_vCrossSectionHydraulicRadius[i], 4.0/3.0));
            }

        }
    }
}

//======================================================================================================================
//! Calculate F terms
//======================================================================================================================
void CSimulation::calculateFlowTerms() {

    for (int i = 0; i < m_nCrossSectionsNumber; i++) {
        m_vCrossSectionF0[i] = m_vCrossSectionQ[i];
        m_vCrossSectionF1[i] = (pow(m_vCrossSectionQ[i], 2.0)/m_vCrossSectionArea[i] + G*m_vCrossSectionI1[i])*m_vCrossSectionBeta[i];
    }
}

//======================================================================================================================
//! Calculate Gv terms
//======================================================================================================================
void CSimulation::calculateSourceTerms() {

    for (int i = 0; i < m_nCrossSectionsNumber; i++) {
        m_vCrossSectionGv0[i] = m_vCrossSectionQ[i]/m_vCrossSectionDX[i];
        m_vCrossSectionGv1[i] = G*m_vCrossSectionI2[i] + m_vCrossSectiongAS0[i] - m_vCrossSectiongASf[i];
    }
}

//======================================================================================================================
//! Calculate the predictor
//======================================================================================================================
void CSimulation::calculatePredictor() {
    //! TODO 014: Add sediment density
    for (int i = 0; i < m_nCrossSectionsNumber - 1; i++) {
        m_vPredictedCrossSectionArea[i] = m_vCrossSectionArea[i] - m_dLambda*(m_vCrossSectionF0[i+1] - m_vCrossSectionF0[i]) + m_dTimestep*m_vCrossSectionGv0[i];
        m_vPredictedCrossSectionQ[i] = m_vCrossSectionQ[i] - m_dLambda*(m_vCrossSectionF1[i+1] - m_vCrossSectionF1[i]) + m_dTimestep*m_vCrossSectionGv1[i];
    }
}


//======================================================================================================================
//! Calculate the corrector
//======================================================================================================================
void CSimulation::calculateCorrector() {
    //! TODO 014: Add sediment density
    for (int i = 1; i < m_nCrossSectionsNumber; i++) {
        m_vCorrectedCrossSectionArea[i] = m_vCrossSectionArea[i] - m_dLambda*(m_vCrossSectionF0[i] - m_vCrossSectionF0[i-1]) + m_dTimestep*m_vCrossSectionGv0[i];
        m_vCorrectedCrossSectionQ[i] = m_vCrossSectionQ[i] - m_dLambda*(m_vCrossSectionF1[i] - m_vCrossSectionF1[i-1]) + m_dTimestep*m_vCrossSectionGv1[i];
    }
}

//======================================================================================================================
//! Update Predictor boundaries
//======================================================================================================================
void CSimulation::updatePredictorBoundaries() {
    if (nGetUpwardEstuarineCondition() == 0) {
        //! Open boundary condition
        m_vPredictedCrossSectionQ[0] = m_vPredictedCrossSectionQ[1];
    }
    else if (nGetUpwardEstuarineCondition() == 1) {
        //! Reflected boundary condition
        m_vPredictedCrossSectionQ[0] = m_dUpwardBoundaryValue;
        // m_vPredictedCrossSectionQ[1] = m_dNextUpwardBoundaryValue;
    }
    else {
        //! Water elevation boundary condition
        m_vPredictedCrossSectionArea[0] = m_dUpwardBoundaryValue;
        m_vPredictedCrossSectionArea[1] = m_dNextUpwardBoundaryValue;
    }

    if (nGetDownwardEstuarineCondition() == 0) {
        //! Open flow
        m_vPredictedCrossSectionArea[m_nCrossSectionsNumber-1] = m_vPredictedCrossSectionArea[m_nCrossSectionsNumber-2];
        m_vPredictedCrossSectionQ[m_nCrossSectionsNumber-1] = m_vPredictedCrossSectionQ[m_nCrossSectionsNumber-2];
    }
    else if (nGetDownwardEstuarineCondition() == 1) {
        //! Given the water flow
        m_vPredictedCrossSectionQ[m_nCrossSectionsNumber-1] = m_dDownwardBoundaryValue;
        // m_vPredictedCrossSectionQ[m_nCrossSectionsNumber-2] = m_dNextDownwardBoundaryValue;
    }
    else {
        //! Given tidal elevation seaward and previous to seaward
        m_vPredictedCrossSectionArea[m_nCrossSectionsNumber-1] = m_dDownwardBoundaryValue;
        m_vPredictedCrossSectionArea[m_nCrossSectionsNumber-2] = m_dNextDownwardBoundaryValue;
    }
}

//======================================================================================================================
//! Update Corrector boundaries
//======================================================================================================================
void CSimulation::updateCorrectorBoundaries() {
    m_vCorrectedCrossSectionQ[0] = m_vCorrectedCrossSectionQ[1];
    if (nGetUpwardEstuarineCondition() == 0) {
        //! Open boundary condition
        m_vCorrectedCrossSectionQ[0] = m_vCorrectedCrossSectionQ[1];
    }
    else if (nGetUpwardEstuarineCondition() == 1) {
        //! Reflected boundary condition
        m_vCorrectedCrossSectionQ[0] = m_vPredictedCrossSectionQ[0];
        m_vCorrectedCrossSectionQ[1] = m_vPredictedCrossSectionQ[1];
    }
    else {
        //! Water elevation boundary condition
        m_vCorrectedCrossSectionArea[0] = m_vPredictedCrossSectionArea[0];
        m_vCorrectedCrossSectionArea[1] = m_vPredictedCrossSectionArea[1];
    }

    if (nGetDownwardEstuarineCondition() == 0) {
        //! Open flow
        m_vCorrectedCrossSectionArea[m_nCrossSectionsNumber-1] = m_vCorrectedCrossSectionArea[m_nCrossSectionsNumber-2];
        m_vCorrectedCrossSectionQ[m_nCrossSectionsNumber-1] = m_vCorrectedCrossSectionQ[m_nCrossSectionsNumber-2];
    }
    else if (nGetDownwardEstuarineCondition() == 1) {
        //! Given the water flow
        m_vCorrectedCrossSectionQ[m_nCrossSectionsNumber-1] = m_dDownwardBoundaryValue;
        m_vCorrectedCrossSectionQ[m_nCrossSectionsNumber-2] = m_dNextDownwardBoundaryValue;
    }
    else {
        //! Given tidal elevation seaward and previous to seaward
        m_vCorrectedCrossSectionArea[m_nCrossSectionsNumber-1] = m_dDownwardBoundaryValue;
        m_vCorrectedCrossSectionArea[m_nCrossSectionsNumber-2] = m_dNextDownwardBoundaryValue;
        // m_vCorrectedCrossSectionQ[m_nCrossSectionsNumber-1] = m_vPredictedCrossSectionQ[m_nCrossSectionsNumber-1];
    }
}

//======================================================================================================================
//! Merge Predictor and Corrector
//======================================================================================================================
void CSimulation::mergePredictorCorrector() {
    if (bGetDoMcComarckLimiterFlux()) {
        //! TODO 015: Include TVD-McComarck
    }
    else {
        for (int i = 0; i < m_nCrossSectionsNumber; i++) {
            m_vCrossSectionArea[i] = 0.5*(m_vPredictedCrossSectionArea[i] + m_vCorrectedCrossSectionArea[i]);
            m_vCrossSectionQ[i] = 0.5*(m_vPredictedCrossSectionQ[i] + m_vCorrectedCrossSectionQ[i]);
        }
    }
}

//======================================================================================================================
//! Smooth solution
//======================================================================================================================
void CSimulation::smoothSolution() {
    for (int i = 0; i < m_nCrossSectionsNumber-1; i++) {
        m_vCrossSectionArea[i] = 0.5*(m_vCrossSectionArea[i+1] + m_vCrossSectionArea[i]);
        m_vCrossSectionQ[i] = 0.5*(m_vCrossSectionQ[i+1] + m_vCrossSectionQ[i]);
    }
}

//======================================================================================================================
//! Get the vector data given the name
//======================================================================================================================
vector<double> CSimulation::vGetVariable(string strVariableName) {

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
    else if (strVariableName == "eta") {
        return m_vCrossSectionElevation;
    }
    else if (strVariableName == "beta") {
        return m_vCrossSectionBeta;
    }
    else if (strVariableName == "I1") {
        return m_vCrossSectionI1;
    }
    else if (strVariableName == "I2") {
        return m_vCrossSectionI2;
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
//! Notifies the user that the simulation has ended, asks for keypress if necessary, and if compiled under GNU can send an email
//===============================================================================================================================
void CSimulation::DoSimulationEnd(int const nRtn)
{};
