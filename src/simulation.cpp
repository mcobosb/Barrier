#include <iostream>
#include <main.h>
#include <utils.h>
#include <windows.h>
using std::cerr;
using std::cout;
using std::endl;
using std::ios;

#include <io.h>

#include <iomanip>
using std::setprecision;
using std::setw;

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

//===============================================================================================================================
//! The CSimulation constructor
//===============================================================================================================================
CSimulation::CSimulation ()
{
    m_dSimDuration =
    m_dSimTimestep = 0.0;

    m_nInitialEstuarineCondition =
    m_nLogFileDetail =
    m_nTimeLogId =
    m_nTimeId =
    m_nPredictor = 0;

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
    writer.nDefineNetCDFFile(this);

    //! Check if any part of the estuary is dry
    if (bGetDoDryBed())
        dryArea();

    // nRtn = writer.nSetOutputData(this);

    m_nTimeId = 0;
    //! Save initial timestep
    writer.nSetOutputData(this);
    m_nTimeId++;

    int m_nStep = 1;
    //! TODO 013: Hydrographs along the estuary have to be reduced by dX of the input location
    while (m_dCurrentTime < m_dSimDuration)
    {
        LogStream << "Step No. " << m_nStep << " - Simulation time: " << m_dCurrentTime << endl;
        // AnnounceProgress();

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

        //==============================================================================================================
        //! STEP 03: update the following timestep
        //==============================================================================================================
        mergePredictorCorrector();

        smoothSolution();

        m_nPredictor = 0;
        //! Check if any part of the estuary is dry (corrector)
        if (bGetDoDryBed())
            dryArea();

        //! Compute the salinity gradient
        // if (bGetDoWaterDensity())
        // }

        if (m_bSaveTime || (m_nLogFileDetail == 2)) {
            writer.nSetOutputData(this);
        }

        //! Increase counter
        m_dCurrentTime = m_dCurrentTime + m_dTimestep;
        if (m_bSaveTime) {
            m_nTimeId++;
        }
        m_nStep++;
    }

    writer.nCloseNetCDFFile(this);
    LogStream << "Step No. " << m_nStep << " - Simulation time: " << m_dCurrentTime << endl;
    return false;
}


//======================================================================================================================
//! Initialize vectorial variables
//======================================================================================================================
void CSimulation::initializeVectors()
{
    m_nPredictor = -1;
    m_nTimeLogId = 0;

    const double nCrossSectionsNumber = this->m_nCrossSectionsNumber;
    const vector<double> vZeros(static_cast<size_t>(nCrossSectionsNumber), 0.0);
    const vector<double> vOnes(static_cast<size_t>(nCrossSectionsNumber), 1.0);

    m_vCrossSectionArea =
    m_vPredictedCrossSectionArea =
    m_vCorrectedCrossSectionArea =
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
    m_vCrossSection_gAS0 =
    m_vCrossSection_gASf =
    m_vCrossSectionF0 =
    m_vCrossSectionF1 =
    m_vCrossSectionGv0 =
    m_vCrossSectionGv1 =
    m_vCrossSectionMurilloFactor = vZeros;

    m_vCrossSectionD1Factor =
    m_vCrossSectionD2Factor = vOnes;

    int nTimestepsNumber = static_cast<int> (m_dSimDuration/m_dSimTimestep) + 1;
    m_vOutputTimesIds = vector<int>(nTimestepsNumber, 0);
    for (int i = 0; i < nTimestepsNumber; i++)
    {
        m_vOutputTimesIds[i] = i;
    }

    for (int i = 0; i < nTimestepsNumber ; i++) {
        m_vOutputTimes.push_back(static_cast<double>(i)*m_dSimTimestep);
    }

    m_dCurrentTime = 0.0;
}

//======================================================================================================================
//! Calculate the initial bed slope
//======================================================================================================================
void CSimulation::calculateBedSlope() {

    //Calculate S0
    for (int i = 0; i < m_nCrossSectionsNumber ; i++) {
        if (i == 0) {
            // Compute dx as forward difference
            double dX1 = this->estuary[1].dGetX() - this->estuary[0].dGetX();
            //! TODO 009: check the sign of this->estuary[i+1].dGetZ() - this->estuary[i].dGetZ()
            this->m_vCrossSectionBedSlope[i] = (this->estuary[0].dGetZ() - this->estuary[1].dGetZ())/dX1;
            //! Save dX into a vector
            m_vCrossSectionDX[i] = dX1;
        }
        else if (i == m_nCrossSectionsNumber - 1){
            // Compute dx as backward difference
            double dX2 = this->estuary[i].dGetX() - this->estuary[i-1].dGetX();
            this->m_vCrossSectionBedSlope[i] = (this->estuary[i-1].dGetZ() - this->estuary[i].dGetZ())/dX2;
            //! Save dX into a vector (last node backward)
            m_vCrossSectionDX[i] = dX2;
        }
        else {
            double dX1 = this->estuary[i+1].dGetX() - this->estuary[i].dGetX();
            double dX2 = this->estuary[i].dGetX() - this->estuary[i-1].dGetX();
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
        //! Along estuary water flow given
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
            this->m_vCrossSectionElevation.push_back(m_vCrossSectionArea[i]/m_vCrossSectionWidth[i]);
        }
    }
    else if (m_nInitialEstuarineCondition == 2) {
        //! Along estuary elevation given
        for (int i = 0; i < m_nCrossSectionsNumber; i++) {
            // this->m_vCrossSectionArea.push_back(this->m_vCrossSectionElevation[i]);
            vector<double> vCrossSectionAreaTmp = this->estuary[i].vGetArea();
            vector<double> vCrossSectionElevationTmp = this->estuary[i].vGetElevation();
            vector<double> vCrossSectionHydraulicRadiusTmp = this->estuary[i].vGetHydraulicRadius();

            //! m_vCrossSectionElevation is the elevation from the water depth m_dZ of every cross-section
            m_vCrossSectionArea[i] = linearInterpolation1d(this->m_vCrossSectionElevation[i] - this->estuary[i].dGetZ(),vCrossSectionElevationTmp, vCrossSectionAreaTmp);
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
    //! Found the interval bettween [x[i], x[i+1]]
    for (size_t i = 0; i < vX.size() - 1; i++) {
        if (dValue >= vX[i] && dValue <= vX[i+1]) {
            // Do the linear interpolation
            const double slope = (vY[i+1] - vY[i]) / (vX[i+1] - vX[i]);
            return vY[i] + slope * (dValue - vX[i]);
        }
    }
    cerr << "LinearInterpolation1d error: value " << dValue << " outside the range (" << vX[0] << ", " << vX[vX.size() - 1] << ")" << endl;
    return true;
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
        const int nElevationSectionsNumber = estuary[i].nGetElevationSectionsNumber();

        if (estuary[i].dGetArea(0) > dArea) {
            //! Take the first node if dArea < the Area of the first elevation node
            vElevationSection.push_back(0);
            getFirstHydraulicParameters(i);
        }
        else if (estuary[i].dGetArea(nElevationSectionsNumber -1) < dArea) {
            //! Take the last node if dArea > the Area of the last elevation node
            vElevationSection.push_back(nCrossSections - 1);
            getLastHydraulicParameters(i);
        }
        else
        {
            for (int j = 0; j < nElevationSectionsNumber; j++) {
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
    double dMinTimestep = 10000000.0;
    double dWaterDensityFactor = 1.0;

    for (int i=0; i< this->m_nCrossSectionsNumber-1; i++) {
        if (m_vCrossSectionArea[i] != DRY_AREA) {
            double dX = (this->estuary[i+1].dGetX() - this->estuary[i].dGetX());
            //! Mean water flow
            m_vCrossSectionU[i] = m_vCrossSectionQ[i]/m_vCrossSectionArea[i];

            //! Perturbation celerity
            m_vCrossSectionC[i] = sqrt(G*m_vCrossSectionArea[i]/m_vCrossSectionWidth[i]);

            //! TODO 003: Courant number due to perturbation on density field
            if (m_bDoWaterDensity) {
                dWaterDensityFactor = m_dCourantNumber*dX/(m_vCrossSectionWidth[i]*m_dBetaSalinityConstant);
            }

            //! Compute the timestep given the Courant number


            if (const double dTimestepTmp = m_dCourantNumber*dX/(fabs(m_vCrossSectionU[i] + m_vCrossSectionC[i])); dTimestepTmp < dMinTimestep) {
                dMinTimestep = dTimestepTmp;
            }

            if (m_bDoWaterDensity && (dWaterDensityFactor < dMinTimestep)) {
                dMinTimestep = dWaterDensityFactor;
            }
        }
        else {
            //! Mean water flow
            m_vCrossSectionU[i] = DRY_Q/DRY_AREA;

            //! Perturbation celerity
            m_vCrossSectionC[i] = sqrt(G*DRY_AREA/m_vCrossSectionWidth[i]);
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
    if ((m_nTimeId < static_cast<int>(m_vOutputTimes.size())) && (m_dCurrentTime + m_dTimestep > m_vOutputTimes[m_nTimeId])) {

        // double dNextTimestep = m_vOutputTimes[m_nTimeId+1];
        // if (dNextTimestep - (m_dTimestep + m_dCurrentTime) < 0.0) {
        m_dTimestep = m_vOutputTimes[m_nTimeId] - m_dCurrentTime;
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
            double dh_dx;
            double dI1dx;

            if (i == 0) {
                // Upward finite-difference
                dh_dx = (this->estuary[i+1].dGetElevation(j) - this->estuary[i].dGetElevation(j))/dX1;
                dI1dx = (this->estuary[i+1].dGetI1(j) - this->estuary[i].dGetI1(j))/dX1;
            }
            else if (i == m_nCrossSectionsNumber - 1) {
                // Downward finite-difference
                dh_dx = (this->estuary[i].dGetElevation(j) - this->estuary[i-1].dGetElevation(j))/dX2;
                dI1dx = (this->estuary[i].dGetI1(j) - this->estuary[i-1].dGetI1(j))/dX2;
            }
            else {
                // Central finite-difference
                dh_dx = (this->estuary[i+1].dGetElevation(j) - this->estuary[i-1].dGetElevation(j))/(dX1 +dX2);
                dI1dx = (this->estuary[i+1].dGetI1(j) - this->estuary[i-1].dGetI1(j))/(dX1 + dX2);

            }
            const double dValue = dI1dx - this->estuary[i].dGetArea(j)*dh_dx;
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
        m_dUpwardBoundaryValue = linearInterpolation1d(m_dCurrentTime, m_vUpwardBoundaryConditionTime, m_vUpwardBoundaryConditionValue);
        m_dNextUpwardBoundaryValue = m_dUpwardBoundaryValue;
    }
    else if (this->nGetUpwardEstuarineCondition() == 2) {
        //! Elevation given upstream
        double dUpwardBoundaryValue = linearInterpolation1d(m_dCurrentTime, m_vUpwardBoundaryConditionTime, m_vUpwardBoundaryConditionValue);
        m_dUpwardBoundaryValue = linearInterpolation1d(-estuary[0].dGetZ()  + dUpwardBoundaryValue, estuary[0].vGetElevation(), estuary[0].vGetArea());
        m_dNextUpwardBoundaryValue = linearInterpolation1d(-estuary[1].dGetZ()  + dUpwardBoundaryValue, estuary[1].vGetElevation(), estuary[1].vGetArea());
    }
    else {};

    //! Calculate downward boundary condition at time t
    if (this->nGetDownwardEstuarineCondition() == 1) {
        //! Water flow given upstream
        m_dDownwardBoundaryValue = linearInterpolation1d(m_dCurrentTime, m_vDownwardBoundaryConditionTime, m_vDownwardBoundaryConditionValue);
        m_dNextDownwardBoundaryValue = m_dDownwardBoundaryValue;
    }
    else if (this->nGetDownwardEstuarineCondition() == 2) {
        //! Elevation given upstream
        double dDownwardBoundaryValue = linearInterpolation1d(m_dCurrentTime, m_vDownwardBoundaryConditionTime, m_vDownwardBoundaryConditionValue);
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
    else if (m_nPredictor == 2) {
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
            if (m_vCrossSectionArea[i] <= DRY_AREA) {
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
        for (int i = 0; i < this->m_nCrossSectionsNumber; i++)
        {
            const double dValue = CDX_MURILLO * sqrt(
                                      2 * pow(m_vCrossSectionHydraulicRadius[i], 2.0 / 3.0) / (G * m_vCrossSectionDX[i]));

            if (dValue < estuary[i].dGetManningNumber()) {
                m_vCrossSectionMurilloFactor[i] = dValue/m_vCrossSectionManningNumber[i];
                m_vCrossSectionManningNumber[i] = dValue;
                m_vCrossSectionI1[i] = m_vCrossSectionI1[i]*m_vCrossSectionMurilloFactor[i] ;
            }
            else {
                m_vCrossSectionMurilloFactor[i] = 1.0;
            }
        }
    }
    else {
        for (int i = 0; i < this->m_nCrossSectionsNumber; i++) {
            m_vCrossSectionI1[i] = m_vCrossSectionI1[i]*m_vCrossSectionMurilloFactor[i];
        }
    }

}

//======================================================================================================================
//! Calculate gAS terms
//======================================================================================================================
void CSimulation::calculate_GS_A_terms() {

    //! Check if do surface term balance
    if (bGetDoSurfaceTermBalance()) {
        double dMeanHydraulicRadius;
        double dMeanQ;
        double dMeanArea;

        //! Average between sections
        if (m_nPredictor == 1) {
            for (int i = 0; i < m_nCrossSectionsNumber; i++) {
                if (i < m_nCrossSectionsNumber - 1) {
                    dMeanArea = (m_vCrossSectionArea[i] + m_vCrossSectionArea[i+1])/2.0;
                    dMeanQ = (m_vCrossSectionQ[i] + m_vCrossSectionQ[i+1])/2.0;
                    dMeanHydraulicRadius = (m_vCrossSectionHydraulicRadius[i] + m_vCrossSectionHydraulicRadius[i+1])/2.0;
                }
                else {
                    dMeanArea = m_vCrossSectionArea[i];
                    dMeanQ = m_vCrossSectionQ[i];
                    dMeanHydraulicRadius = m_vCrossSectionHydraulicRadius[i];
                }
                m_vCrossSection_gAS0[i] = G*dMeanArea*m_vCrossSectionBedSlope[i];
                m_vCrossSectionFrictionSlope[i] = dMeanQ*fabs(dMeanQ)/(dMeanArea*pow(dMeanHydraulicRadius, 4.0/3.0));
                m_vCrossSection_gASf[i] = G*pow(m_vCrossSectionManningNumber[i], 2.0)*m_vCrossSectionFrictionSlope[i];
            }
        }
        else
        {
            for (int i = 0; i < m_nCrossSectionsNumber; i++) {
                if (i > 0) {
                    dMeanArea = (m_vPredictedCrossSectionArea[i] + m_vPredictedCrossSectionArea[i-1])/2.0;
                    dMeanQ = (m_vPredictedCrossSectionQ[i] + m_vPredictedCrossSectionQ[i-1])/2.0;
                    dMeanHydraulicRadius = (m_vCrossSectionHydraulicRadius[i] + m_vCrossSectionHydraulicRadius[i-1])/2.0;
                }
                else {
                    dMeanArea = m_vPredictedCrossSectionArea[i];
                    dMeanQ = m_vPredictedCrossSectionQ[i];
                    dMeanHydraulicRadius = m_vCrossSectionHydraulicRadius[i];
                }
                m_vCrossSection_gAS0[i] = G*dMeanArea*m_vCrossSectionBedSlope[i];
                m_vCrossSectionFrictionSlope[i] = dMeanQ*fabs(dMeanQ)/(dMeanArea*pow(dMeanHydraulicRadius, 4.0/3.0));
                m_vCrossSection_gASf[i] = G*pow(m_vCrossSectionManningNumber[i], 2.0)*pow(dMeanQ, 2.0)/(dMeanArea*pow(dMeanHydraulicRadius, 4.0/3.0));
            }
        }

    }
    else {
        for (int i = 0; i < m_nCrossSectionsNumber; i++) {
            m_vCrossSection_gAS0[i] = G*m_vCrossSectionArea[i]*m_vCrossSectionArea[i];
            if (m_nPredictor == 1) {
                m_vCrossSectionFrictionSlope[i] = m_vCrossSectionQ[i]*fabs(m_vCrossSectionQ[i])/(m_vCrossSectionArea[i]*pow(m_vCrossSectionHydraulicRadius[i], 4.0/3.0));
                m_vCrossSection_gASf[i] = G*pow(m_vCrossSectionManningNumber[i], 2.0)*m_vCrossSectionFrictionSlope[i];
            }
            else {
                m_vCrossSectionFrictionSlope[i] = m_vPredictedCrossSectionQ[i]*fabs(m_vPredictedCrossSectionQ[i])/(m_vPredictedCrossSectionArea[i]*pow(m_vCrossSectionHydraulicRadius[i], 4.0/3.0));
                m_vCrossSection_gASf[i] = G*pow(m_vCrossSectionManningNumber[i], 2.0)*pow(m_vPredictedCrossSectionQ[i], 2.0)/(m_vPredictedCrossSectionArea[i]*pow(m_vCrossSectionHydraulicRadius[i], 4.0/3.0));
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
            m_vCrossSectionF1[i] = (pow(m_vCrossSectionQ[i], 2.0)/m_vCrossSectionArea[i] + G*m_vCrossSectionI1[i])*m_vCrossSectionBeta[i];
        }
    }
    else {
        for (int i = 0; i < m_nCrossSectionsNumber; i++) {
            m_vCrossSectionF0[i] = m_vPredictedCrossSectionQ[i];
            m_vCrossSectionF1[i] = (pow(m_vPredictedCrossSectionQ[i], 2.0)/m_vPredictedCrossSectionArea[i] + G*m_vCrossSectionI1[i])*m_vCrossSectionBeta[i];
        }

    }
}

//======================================================================================================================
//! Calculate Gv terms
//======================================================================================================================
void CSimulation::calculateSourceTerms() {

    for (int i = 0; i < m_nCrossSectionsNumber; i++) {
        m_vCrossSectionGv0[i] = 0.0; //! TODO 025: verify
        m_vCrossSectionGv1[i] = G*m_vCrossSectionI2[i] + m_vCrossSection_gAS0[i] - m_vCrossSection_gASf[i];
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
        //! TODO 028: verify if  m_vCrossSectionArea[i] and m_vCrossSectionQ[i] might be Predicted
        m_vCorrectedCrossSectionArea[i] = m_vPredictedCrossSectionArea[i] - m_dLambda*(m_vCrossSectionF0[i] - m_vCrossSectionF0[i-1]) + m_dTimestep*m_vCrossSectionGv0[i];
        m_vCorrectedCrossSectionQ[i] = m_vPredictedCrossSectionQ[i] - m_dLambda*(m_vCrossSectionF1[i] - m_vCrossSectionF1[i-1]) + m_dTimestep*m_vCrossSectionGv1[i];
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
        m_vPredictedCrossSectionQ[m_nCrossSectionsNumber-1] = m_vPredictedCrossSectionQ[m_nCrossSectionsNumber-2];
    }
}

//======================================================================================================================
//! Update Corrector boundaries
//======================================================================================================================
void CSimulation::updateCorrectorBoundaries() {
    // m_vCorrectedCrossSectionQ[0] = m_vCorrectedCrossSectionQ[1];
    if (nGetUpwardEstuarineCondition() == 0) {
        //! Open boundary condition
        m_vCorrectedCrossSectionQ[0] = m_vCorrectedCrossSectionQ[1];
        //! Corrector does not calculate i = 0, it is imposed
        m_vCorrectedCrossSectionArea[0] = m_vPredictedCrossSectionArea[0];
    }
    else if (nGetUpwardEstuarineCondition() == 1) {
        //! Reflected boundary condition
        m_vCorrectedCrossSectionQ[0] = m_dUpwardBoundaryValue;
        //! Corrector does not calculate i = 0, it is imposed
        m_vCorrectedCrossSectionArea[0] = m_vPredictedCrossSectionArea[0];
    }
    else {
        //! Water elevation boundary condition
        m_vCorrectedCrossSectionArea[0] = m_vPredictedCrossSectionArea[0];
        m_vCorrectedCrossSectionArea[1] = m_vPredictedCrossSectionArea[1];
        //! Corrector does not calculate i = 0, it is imposed
        m_vCorrectedCrossSectionQ[0] = m_vPredictedCrossSectionQ[0];
    }

    if (nGetDownwardEstuarineCondition() == 0) {
        //! Open flow
        m_vCorrectedCrossSectionArea[m_nCrossSectionsNumber-1] = m_vCorrectedCrossSectionArea[m_nCrossSectionsNumber-2];
        m_vCorrectedCrossSectionQ[m_nCrossSectionsNumber-1] = m_vCorrectedCrossSectionQ[m_nCrossSectionsNumber-2];
    }
    else if (nGetDownwardEstuarineCondition() == 1) {
        //! Given the water flow
        m_vCorrectedCrossSectionQ[m_nCrossSectionsNumber-1] = m_dDownwardBoundaryValue;
        // m_vCorrectedCrossSectionQ[m_nCrossSectionsNumber-2] = m_dNextDownwardBoundaryValue;
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
    if (bGetDoSurfaceGradientMethod()) {
        for (int i = 0; i < m_nCrossSectionsNumber; i++) {
            m_vCrossSectionElevation[i] = m_vCrossSectionElevation[i] - estuary[i].dGetZ();
        }
    }
    if (bGetDoMcComarckLimiterFlux()) {
        //! Include TVD-McComarck?
        const double nCrossSectionsNumber = this->m_nCrossSectionsNumber;
        vector<double> a1_med (static_cast<size_t>(nCrossSectionsNumber), 0.0);
        vector<double> a2_med (static_cast<size_t>(nCrossSectionsNumber), 0.0);

        vector<double> alfa1_med (static_cast<size_t>(nCrossSectionsNumber), 0.0);
        vector<double> alfa2_med (static_cast<size_t>(nCrossSectionsNumber), 0.0);

        vector<double> psi1_med (static_cast<size_t>(nCrossSectionsNumber), 0.0);
        vector<double> psi2_med (static_cast<size_t>(nCrossSectionsNumber), 0.0);

        vector<double> r1_med (static_cast<size_t>(nCrossSectionsNumber), 0.0);
        vector<double> r2_med (static_cast<size_t>(nCrossSectionsNumber), 0.0);

        vector<double> fi1_med (static_cast<size_t>(nCrossSectionsNumber), 0.0);
        vector<double> fi2_med (static_cast<size_t>(nCrossSectionsNumber), 0.0);

        vector<double> vFactor1 (static_cast<size_t>(nCrossSectionsNumber), 0.0);
        vector<double> vFactor2 (static_cast<size_t>(nCrossSectionsNumber), 0.0);


        //! García-Navarro Psi formula
        double delta1 = m_dDeltaValue;
        double delta2 = m_dDeltaValue;
        for (int i = 1; i < m_nCrossSectionsNumber; i++)
        {
            const double u_med = (m_vCrossSectionQ[i] / sqrt(m_vCrossSectionArea[i]) + m_vCrossSectionQ[i - 1] /
                                  sqrt(m_vCrossSectionArea[i - 1])) / (
                                     sqrt(m_vCrossSectionArea[i]) + sqrt(m_vCrossSectionArea[i - 1]));
            const double A_med = (m_vCrossSectionArea[i] + m_vCrossSectionArea[i - 1]) / 2.0;
            const double c_med = (m_vCrossSectionC[i] + m_vCrossSectionC[i - 1]) / 2.0;
            a1_med[i-1] = u_med + c_med;
            a2_med[i-1] = u_med - c_med;

            if (!bGetDoSurfaceGradientMethod()) {
                alfa1_med[i-1] = ((m_vCrossSectionQ[i] - m_vCrossSectionQ[i-1]) - a2_med[i-1]*(A_med - m_vCrossSectionArea[i-1]))/(2.0*c_med);
                alfa2_med[i-1] = -((m_vCrossSectionQ[i] - m_vCrossSectionQ[i-1]) - a1_med[i-1]*(A_med - m_vCrossSectionArea[i-1]))/(2.0*c_med);
            }
            else {
                alfa1_med[i-1] = ((m_vCrossSectionQ[i]/m_vCrossSectionWidth[i] - m_vCrossSectionQ[i-1]/m_vCrossSectionWidth[i-1])-a2_med[i-1]*(m_vCrossSectionElevation[i] - m_vCrossSectionElevation[i-1]))/(2.0*c_med);
                alfa2_med[i-1] = -((m_vCrossSectionQ[i]/m_vCrossSectionWidth[i] - m_vCrossSectionQ[i-1]/m_vCrossSectionWidth[i-1])-a1_med[i-1]*(m_vCrossSectionElevation[i] - m_vCrossSectionElevation[i-1]))/(2.0*c_med);
            }

            if (nGetPsiFormula() != 1) {
                //! Tseng Psi formula
                vector<double> deltaValues = {0.0, a1_med[i] - (m_vCrossSectionU[i-1] + m_vCrossSectionC[i-1]), m_vCrossSectionU[i] + m_vCrossSectionC[i] - a1_med[i]};
                delta1 = dMaxVectorValue(deltaValues);

                deltaValues ={0.0, a1_med[i] - (m_vCrossSectionU[i-1] + m_vCrossSectionC[i-1]), m_vCrossSectionU[i] + m_vCrossSectionC[i] - a1_med[i]};
                delta2 = dMaxVectorValue(deltaValues);
            }

            //! Computing psi_i
            if (abs(a1_med[i-1])>= delta1)
            {
                psi1_med[i-1] = abs(a1_med[i-1]);
            }
            else
            {
                psi1_med[i-1] = delta1;
            }

            if (abs(a2_med[i-1])>= delta2)
            {
                psi2_med[i-1] = abs(a2_med[i-1]);
            }
            else
            {
                psi2_med[i-1] = delta2;
            }
        }


        for (int i = 1; i  < m_nCrossSectionsNumber; i++)
        {
            if (alfa1_med[i] == 0.0) {
                r1_med[i-1] = 1.0;
            }
            else {
                //! Computing ri
                if (a1_med[i-1] < 0)
                {
                    r1_med[i-1] = alfa1_med[i]/alfa1_med[i-1];
                }
                else if (a1_med[i-1] == 0.0)
                {
                    r1_med[i-1] = 1.0;
                }
                else if ((i != 1) && (a1_med[i-1] > 0))
                {
                    r1_med[i-1] = alfa1_med[i-1]/alfa1_med[i];
                }
            }

            if (alfa2_med[i] == 0.0) {
                r2_med[i-1] = 1.0;
            }
            else {
                if (a2_med[i-1] < 0)
                {
                    r2_med[i-1] = alfa2_med[i]/alfa2_med[i-1];
                }
                else if (alfa2_med[i-1] == 0.0)
                {
                    r2_med[i-1] = 1.0;
                }
                else if ((i != 1) && (a2_med[i-1] > 0))
                {
                    r2_med[i-1] = alfa2_med[i-1]/alfa2_med[i];
                }
            }

            //! Computing fi_i
            if (nGetEquationLimiterFlux() == 1)
            {
                //! MinMod
                fi1_med[i-1] = dMaxVectorValue({0.0, dMinVectorValue({1., r1_med[i-1]})});
                fi2_med[i-1] = dMaxVectorValue({0.0, dMinVectorValue({1., r2_med[i-1]})});
            }
            else if (nGetEquationLimiterFlux() == 2)
            {
                //! Roe's Superbee
                fi1_med[i-1] = dMaxVectorValue({0.0, dMinVectorValue({2*r1_med[i-1], 1.}), dMinVectorValue({r1_med[i-1], 2.})});
                fi2_med[i-1] = dMaxVectorValue({0.0, dMinVectorValue({2*r2_med[i-1], 1.}), dMinVectorValue({r2_med[i-1], 2.})});
            }
            else if (nGetEquationLimiterFlux() == 3)
            {
                //! Van Leer
                fi1_med[i-1] = (abs(r1_med[i-1]) + r1_med[i-1])/(1+abs(r1_med[i-1]));
                fi2_med[i-1] = (abs(r2_med[i-1]) + r2_med[i-1])/(1+abs(r2_med[i-1]));
            }
            else if (nGetEquationLimiterFlux() == 4)
            {
                fi1_med[i-1] = (r1_med[i-1]*r1_med[i-1] + r1_med[i-1])/(1+r1_med[i-1]*r1_med[i-1]);
                fi2_med[i-1] = (r2_med[i-1]*r2_med[i-1] + r2_med[i-1])/(1+r2_med[i-1]*r2_med[i-1]);
            }

            //! Computing factor
            vFactor1[i-1] = alfa1_med[i-1]*psi1_med[i-1]*(1-m_dLambda*abs(a1_med[i-1]))*(1-fi1_med[i-1]);
            vFactor2[i-1] = alfa2_med[i-1]*psi2_med[i-1]*(1-m_dLambda*abs(a2_med[i-1]))*(1-fi2_med[i-1]);

            m_vCrossSectionD1Factor[i-1] = 0.5*(vFactor1[i-1] + vFactor2[i-1]);
            m_vCrossSectionD2Factor[i-1] = 0.5*(vFactor1[i-1]*a1_med[i-1] + vFactor2[i-1]*a2_med[i-1]);
        }
        // m_vCrossSectionD1Factor[0] = m_vCrossSectionD1Factor[1];
        // m_vCrossSectionD2Factor[0] = m_vCrossSectionD2Factor[1];
    }
    for (int i = 0; i < m_nCrossSectionsNumber-1; i++) {
        m_vCrossSectionArea[i] = 0.5*(m_vPredictedCrossSectionArea[i] + m_vCorrectedCrossSectionArea[i]) +m_dLambda*(m_vCrossSectionD1Factor[i+1] - m_vCrossSectionD1Factor[i]);
        m_vCrossSectionQ[i] = 0.5*(m_vPredictedCrossSectionQ[i] + m_vCorrectedCrossSectionQ[i]) +m_dLambda*(m_vCrossSectionD2Factor[i+1] - m_vCrossSectionD2Factor[i]);
    }
    m_vCrossSectionArea[m_nCrossSectionsNumber-1] = 0.5*(m_vPredictedCrossSectionArea[m_nCrossSectionsNumber-1] + m_vCorrectedCrossSectionArea[m_nCrossSectionsNumber-1]) +m_dLambda*(m_vCrossSectionD1Factor[m_nCrossSectionsNumber-1] - m_vCrossSectionD1Factor[m_nCrossSectionsNumber-2]);
    m_vCrossSectionQ[m_nCrossSectionsNumber-1] = 0.5*(m_vPredictedCrossSectionQ[m_nCrossSectionsNumber-1] + m_vCorrectedCrossSectionQ[m_nCrossSectionsNumber-1]) +m_dLambda*(m_vCrossSectionD2Factor[m_nCrossSectionsNumber-1] - m_vCrossSectionD2Factor[m_nCrossSectionsNumber-2]);

}

//======================================================================================================================
//! Smooth solution
//======================================================================================================================
void CSimulation::smoothSolution() {
    for (int i = 1; i < m_nCrossSectionsNumber-1; i++) {
        m_vCrossSectionArea[i] = 0.25*(m_vCrossSectionArea[i+1] + 2*m_vCrossSectionArea[i] + m_vCrossSectionArea[i-1]);
        m_vCrossSectionQ[i] = 0.25*(m_vCrossSectionQ[i+1] +2*m_vCrossSectionQ[i]+ m_vCrossSectionQ[i-1]);
    }
    // m_vCrossSectionArea[0] = 0.5*(m_vCrossSectionArea[1] + m_vCrossSectionArea[0]);
    // m_vCrossSectionQ[0] = 0.5*(m_vCrossSectionQ[1] + m_vCrossSectionQ[0]);
    //
    // m_vCrossSectionArea[m_nCrossSectionsNumber-1] = 0.5*(m_vCrossSectionArea[m_nCrossSectionsNumber-1] + m_vCrossSectionArea[m_nCrossSectionsNumber-2]);
    // m_vCrossSectionQ[m_nCrossSectionsNumber-1] = 0.5*(m_vCrossSectionQ[m_nCrossSectionsNumber-1] + m_vCrossSectionQ[m_nCrossSectionsNumber-2]);
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
//! Displays information regarding the progress of the simulation
//===============================================================================================================================
void CSimulation::AnnounceProgress() const {

    // Stdout is connected to a tty, so not running as a background job
    static double sdElapsed = 0;
    static double sdToGo = 0;
    const time_t tNow = time(nullptr);

    // Calculate time elapsed and remaining
    sdElapsed = difftime(tNow, m_tSysStartTime);
    sdToGo = (sdElapsed * m_dSimDuration / m_dCurrentTime) - sdElapsed;

    // Tell the user about progress (note need to make several separate calls to cout here, or MS VC++ compiler appears to get confused)


    if (m_bSaveTime) {
        cout << SIMULATING;
        cout << std::fixed << setprecision(3) << setw(9) << 100 * m_dCurrentTime / m_dSimDuration;
        cout << setw(9) << "To NetCDF:" + to_string(m_nTimeId);
    }

    // else {
    //     cout << setw(9) << SPACE;
    // }
}

//===============================================================================================================================
//! Notifies the user that the simulation has ended, asks for keypress if necessary, and if compiled under GNU can send an email
//===============================================================================================================================
void CSimulation::DoSimulationEnd(int const nRtn)
{};
