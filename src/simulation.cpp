#include <iostream>
#include <main.h>
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

//===============================================================================================================================
//! The CSimulation constructor
//===============================================================================================================================
CSimulation::CSimulation ()
{
    m_dSimDuration =
    m_dSimTimestep = 0.0;

    m_strInitialEstuarineCondition = "";

    vector<string> m_vOutputVariables;

    vector<CCrossSection> estuary;
    vector<CHydrograph> hydrographs;
}

//===============================================================================================================================
//! The CSimulation destructor
//===============================================================================================================================
CSimulation::~CSimulation () = default;


//! Method for getting the simulation duration
double CSimulation::dGetSimulationDuration() {
    return m_dSimDuration;
}

//! Method for setting the simulation duration
void CSimulation::dSetSimulationDuration(double simDuration) {
    m_dSimDuration = simDuration;
};

//! Method for getting the simulation timestep
double CSimulation::dGetSimulationTimestep() {
    return m_dSimTimestep;
}

//! Method for setting the simulation timestep
void CSimulation::dSetSimulationTimestep(double simTimestep) {
    m_dSimTimestep = simTimestep;
};

//! Method for getting the initial estuarine condition
string CSimulation::strGetInitialEstuarineCondition() {
    return m_strInitialEstuarineCondition;
}

//! Method for setting the initial estuarine condition
void CSimulation::strSetInitialEstuarineCondition(string initialCondition) {
    m_strInitialEstuarineCondition = initialCondition;
};

//! Method for getting the computation of water density
bool CSimulation::bGetDoWaterDensity() {
    return m_bDoWaterDensity;
}

//! Method for setting the computation of water density
void CSimulation::bSetDoWaterDensity(bool computeWaterDensity) {
    m_bDoWaterDensity = computeWaterDensity;
};

//! Method for getting the along channel constant elevation
double CSimulation::dGetBetaSalinityConstant() {
    return m_dBetaSalinityConstant;
}

//! Method for setting the along channel constant elevation
void CSimulation::dSetBetaSalinityConstant(double betaSalinity) {
    m_dBetaSalinityConstant = betaSalinity;
};

//! Method for getting the along channel constant elevation
double CSimulation::dGetLongitudinalDispersionConstant() {
    return m_dLongitudinalDispersion;
}

//! Method for setting the along channel constant elevation
void CSimulation::dSetLongitudinalDispersionConstant(double longitudinalDispersion) {
    m_dLongitudinalDispersion = longitudinalDispersion;
};

//! Method for getting the upward estuarine condition
int CSimulation::nGetUpwardEstuarineCondition() {
    return m_nUpwardEstuarineCondition;
}

//! Method for setting the upward estuarine condition
void CSimulation::nSetUpwardEstuarineCondition(int upwardEstuarineCondition) {
    m_nUpwardEstuarineCondition = upwardEstuarineCondition;
};

//! Method for getting the downward estuarine condition
string CSimulation::strGetDownwardEstuarineCondition() {
    return m_strDownwardEstuarineCondition;
}

//! Method for setting the downward estuarine condition
void CSimulation::strSetDownwardEstuarineCondition(string downwardEstuarineCondition) {
    m_strDownwardEstuarineCondition = downwardEstuarineCondition;
};

//! Method for getting the courant number
double CSimulation::dGetCourantNumber() {
    return m_dCourantNumber;
}

//! Method for setting the courant number
void CSimulation::dSetCourantNumber(double courantNumber) {
    m_dCourantNumber = courantNumber;
}

//! Method for getting if McComarck limiter flux is applied
bool CSimulation::bGetDoMcComarckLimiterFlux() {
    return m_bDoMackComarckLimiterFlux;
}
//! Method for setting if McComarck limiter flux is applied
void CSimulation::bSetDoMcComarckLimiterFlux(bool doMcComarckLimiterFlux) {
    m_bDoMackComarckLimiterFlux = doMcComarckLimiterFlux;
}

//! Method for getting equation limiter flux
int CSimulation::nGetEquationLimiterFlux() {
    return m_nEquationMacComarckLimiterFlux;
}
//! Method for setting equation limiter flux
void CSimulation::nSetEquationLimiterFlux(int equationLimiterFlux) {
    m_nEquationMacComarckLimiterFlux = equationLimiterFlux;
}

//! Method for getting Psi formula
int CSimulation::nGetPsiFormula() {
    return m_nPsiFormula;
}
//! Method for setting Psi formula
void CSimulation::nSetPsiFormula(int psiFormula) {
    m_nPsiFormula = psiFormula;
}

//! Method for getting Delta Value
double CSimulation::dGetDeltaValue() {
    return m_dDeltaValue;
}
//! Method for setting Delta Value
void CSimulation::dSetDeltaValue(double deltaValue) {
    m_dDeltaValue = deltaValue;
}

//! Method for getting if surface gradient method is applied
bool CSimulation::bGetDoSurfaceGradientMethod() {
    return m_bDoSurfaceGradientMethod;
}
//! Method for setting if surface gradient method is applied
void CSimulation::bSetDoSurfaceGradientMethod(bool doSurfaceGradientMethod) {
    m_bDoSurfaceGradientMethod = doSurfaceGradientMethod;
}

//! Method for getting if source term balance is applied
bool CSimulation::bGetDoSurfaceTermBalance() {
    return m_bDoSourceTermBalance;
}
//! Method for setting if source term balance is applied
void CSimulation::bSetDoSurfaceTermBalance(bool doSourceTermBalance) {
    m_bDoSourceTermBalance = doSourceTermBalance;
}

//! Method for getting if beta coefficient is applied
bool CSimulation::bGetDoBetaCoefficient() {
    return m_bDoBetaCoefficient;
}
//! Method for setting if beta coefficient is applied
void CSimulation::bSetDoBetaCoefficient(bool doBetaCoefficient) {
    m_bDoBetaCoefficient = doBetaCoefficient;
}

//! Method for getting if dry bed is applied
bool CSimulation::bGetDoDryBed() {
    return m_bDoDryBed;
}
//! Method for setting if dry bed is applied
void CSimulation::bSetDoDryBed(bool doDryBed) {
    m_bDoDryBed = doDryBed;
}

//! Method for getting if Murillo condition is applied
bool CSimulation::bGetDoMurilloCondition() {
    return m_bDoMurilloCondition;
}
//! Method for setting if Murillo condition is applied
void CSimulation::bSetDoMurilloCondition(bool doMurilloCondition) {
    m_bDoMurilloCondition = doMurilloCondition;
}

//===============================================================================================================================
//! Appends output variables
//===============================================================================================================================
void CSimulation::strAddOutputVariable(string strItem){
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
int CSimulation::nGetHydrographsNo(){
    return m_nHydrographsNumber;
}
//! Method for setting the number of hydrographs
void CSimulation::nSetHydrographsNo(int nValue) {
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


    double dTimeTmp = 0.0;
    double dTimestepTmp =  0.0;
    double dArea = 0.0;

    calculateBedSlope();
    calculateAlongEstuaryInitialConditions();

    while (dTimeTmp < m_dSimDuration)
    {

        //! Calculate the hydraulic parameters for prediction
        calculateHydraulicParameters(true);

        //! TODO 003: Check if the entire estuary is dry
        calculateCourantNumber();

        calculateIs();



        //! Increase counter
        dTimeTmp = dTimeTmp + dTimestepTmp;
    }

    return false;
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
            dX1 = this->estuary[i+1].dGetX() - this->estuary[i].dGetX();
            //! TODO 009: check the sign of this->estuary[i+1].dGetZ() - this->estuary[i].dGetZ()
            this->m_vCrossSectionBedSlope.push_back((this->estuary[i+1].dGetZ() - this->estuary[i].dGetZ())/dX1);
        }
        else if (i == m_nCrossSectionsNumber - 1){
            // Compute dx as backward difference
            dX2 = this->estuary[i].dGetX() - this->estuary[i-1].dGetX();
            this->m_vCrossSectionBedSlope.push_back((this->estuary[i].dGetZ() - this->estuary[i-1].dGetZ())/dX2);
        }
        else {
            dX1 = this->estuary[i+1].dGetX() - this->estuary[i].dGetX();
            dX2 = this->estuary[i].dGetX() - this->estuary[i-1].dGetX();
            this->m_vCrossSectionBedSlope.push_back((this->estuary[i+1].dGetZ() - this->estuary[i-1].dGetZ())/(dX1 + dX2));
        }
    }
}

//======================================================================================================================
//! Calculate the initial conditions along the estuary
//======================================================================================================================
void CSimulation::calculateAlongEstuaryInitialConditions() {

    if (m_strInitialEstuarineCondition == "Q") {
        for (int i = 0; i < m_nCrossSectionsNumber; i++) {
            double dManningFactor = this->m_vCrossSectionQ[i]*this->estuary[i].dGetManningNumber()/(sqrt(this->m_vCrossSectionBedSlope[i]));

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
    else {
        for (int i = 0; i < m_nCrossSectionsNumber; i++) {
            this->m_vCrossSectionArea[i] = this->m_vCrossSectionElevation[i];
            vector<double> vCrossSectionAreaTmp = this->estuary[i].vGetArea();
            vector<double> vCrossSectionElevationTmp = this->estuary[i].vGetElevation();
            vector<double> vCrossSectionHydraulicRadiusTmp = this->estuary[i].vGetHydraulicRadius();

            //! m_vCrossSectionElevation is the elevation from the water depth m_dZ of every cross-section
            this->m_vCrossSectionArea[i] = linearInterpolation1d(this->m_vCrossSectionElevation[i],vCrossSectionElevationTmp, vCrossSectionAreaTmp);
            this->m_vCrossSectionHydraulicRadius[i] = linearInterpolation1d(this->m_vCrossSectionArea[i], vCrossSectionAreaTmp, vCrossSectionHydraulicRadiusTmp);

            //! Compute Q given the area
            this->m_vCrossSectionQ[i] = m_vCrossSectionArea[i]*sqrt(this->m_vCrossSectionBedSlope[i])*pow(vCrossSectionHydraulicRadiusTmp[i], 2.0/3.0)/this->estuary[i].dGetManningNumber();
        }
    }
}


//======================================================================================================================
//! Interpolate function
//======================================================================================================================
double CSimulation::linearInterpolation1d(double dValue, vector<double> vX, vector<double> vY) {
    //! Encuentra el intervalo adecuado [x[i], x[i+1]]
    for (size_t i = 0; i < vX.size() - 1; ++i) {
        if (dValue >= vX[i] && dValue <= vX[i+1]) {
            // Do the linear interpolation
            double slope = (vY[i+1] - vY[i]) / (vX[i+1] - vX[i]);
            return vY[i] + slope * (dValue - vX[i]);
        }
        else {
            //! TODO 005: Return an interpolation error
            return 0.0;
        }
    }
}
//======================================================================================================================
//! Obtain the hydraulic parameters from the cross-sectional data
//======================================================================================================================
void CSimulation::calculateHydraulicParameters(const bool bPredictedParameters) {
    //! Number of estuarine cross-sections
    const int nCrossSections = m_nCrossSectionsNumber;
    vector<int> vElevationSection;

    for (int i = 0; i < nCrossSections; i++) {
        //! TODO 008: Create getter and setter for cross-section area, hydraulic radius, elevation and slope
        double dArea = this->m_vCrossSectionArea[i];

        if (estuary[i].dGetArea(0) > dArea) {
            //! Take the first node if dArea < the Area of the first elevation node
            vElevationSection.push_back(0);
            getFirstHydraulicParameters(i, bPredictedParameters);
            break;
        }
        else if (estuary[i].dGetArea(nCrossSections -1) < dArea) {
            //! Take the last node if dArea > the Area of the last elevation node
            vElevationSection.push_back(nCrossSections - 1);
            getLastHydraulicParameters(i, bPredictedParameters);
            break;
        }

        for (int j = 0; j < estuary[i].nGetElevationSectionsNumber()-1; j++) {
            //! Getting the elevations node which Area is below dArea and next node Area higher than dArea
            if ((estuary[i].dGetArea(j) < dArea) && (estuary[i].dGetElevation(j+1) > dArea)) {
                vElevationSection.push_back(j);
                interpolateHydraulicParameters(dArea, i, j, bPredictedParameters);
                break;
            }
        }
    }
}

//======================================================================================================================
//! Interpolate the Hydraulic parameters between the elevation node given.
//======================================================================================================================
void CSimulation::interpolateHydraulicParameters(double dArea, int nCrossSection, int nElevationNode, bool bPredictedParameters) {
    double dInterpolationFactor = (dArea - this->estuary[nCrossSection].dGetArea(nElevationNode)) / (this->estuary[nCrossSection].dGetArea(nElevationNode+1) - this->estuary[nCrossSection].dGetArea(nElevationNode));
    if (bPredictedParameters) {
        m_vPredictedCrossSectionHydraulicRadius[nCrossSection] = dInterpolationFactor*(this->estuary[nCrossSection].dGetHydraulicRadius(nElevationNode+1) - this->estuary[nCrossSection].dGetHydraulicRadius(nElevationNode)) + this->estuary[nCrossSection].dGetHydraulicRadius(nElevationNode);
        m_vPredictedCrossSectionWidth[nCrossSection] =  dInterpolationFactor*(this->estuary[nCrossSection].dGetWidth(nElevationNode+1) - this->estuary[nCrossSection].dGetWidth(nElevationNode)) + this->estuary[nCrossSection].dGetWidth(nElevationNode);
        m_vPredictedCrossSectionElevation[nCrossSection] = dInterpolationFactor*(this->estuary[nCrossSection].dGetElevation(nElevationNode+1) - this->estuary[nCrossSection].dGetElevation(nElevationNode)) + this->estuary[nCrossSection].dGetElevation(nElevationNode);
        m_vPredictedCrossSectionBeta[nCrossSection] = dInterpolationFactor*(this->estuary[nCrossSection].dGetBeta(nElevationNode+1) - this->estuary[nCrossSection].dGetElevation(nElevationNode)) + this->estuary[nCrossSection].dGetElevation(nElevationNode);
        m_vPredictedCrossSectionI1[nCrossSection] = dInterpolationFactor*(this->estuary[nCrossSection].dGetI1(nElevationNode+1) - this->estuary[nCrossSection].dGetElevation(nElevationNode)) + this->estuary[nCrossSection].dGetElevation(nElevationNode);
        m_vPredictedCrossSectionI2[nCrossSection] = dInterpolationFactor*(this->estuary[nCrossSection].dGetI2(nElevationNode+1) - this->estuary[nCrossSection].dGetElevation(nElevationNode)) + this->estuary[nCrossSection].dGetElevation(nElevationNode);
        m_vPredictedCrossSectionLeftRBLocation[nCrossSection] = dInterpolationFactor*(this->estuary[nCrossSection].dGetLeftY(nElevationNode+1) - this->estuary[nCrossSection].dGetElevation(nElevationNode)) + this->estuary[nCrossSection].dGetElevation(nElevationNode);
        m_vPredictedCrossSectionRightRBLocation[nCrossSection] = dInterpolationFactor*(this->estuary[nCrossSection].dGetRightY(nElevationNode+1) - this->estuary[nCrossSection].dGetElevation(nElevationNode)) + this->estuary[nCrossSection].dGetElevation(nElevationNode);


    }
    else {
        m_vCorrectedCrossSectionHydraulicRadius[nCrossSection] = dInterpolationFactor*(this->estuary[nCrossSection].dGetHydraulicRadius(nElevationNode+1) - this->estuary[nCrossSection].dGetHydraulicRadius(nElevationNode)) + this->estuary[nCrossSection].dGetHydraulicRadius(nElevationNode);
        m_vCorrectedCrossSectionWidth[nCrossSection] =  dInterpolationFactor*(this->estuary[nCrossSection].dGetWidth(nElevationNode+1) - this->estuary[nCrossSection].dGetWidth(nElevationNode)) + this->estuary[nCrossSection].dGetWidth(nElevationNode);
        m_vCorrectedCrossSectionElevation[nCrossSection] = dInterpolationFactor*(this->estuary[nCrossSection].dGetElevation(nElevationNode+1) - this->estuary[nCrossSection].dGetElevation(nElevationNode)) + this->estuary[nCrossSection].dGetElevation(nElevationNode);
        m_vCorrectedCrossSectionBeta[nCrossSection] = dInterpolationFactor*(this->estuary[nCrossSection].dGetBeta(nElevationNode+1) - this->estuary[nCrossSection].dGetElevation(nElevationNode)) + this->estuary[nCrossSection].dGetElevation(nElevationNode);
        m_vCorrectedCrossSectionI1[nCrossSection] = dInterpolationFactor*(this->estuary[nCrossSection].dGetI1(nElevationNode+1) - this->estuary[nCrossSection].dGetElevation(nElevationNode)) + this->estuary[nCrossSection].dGetElevation(nElevationNode);
        m_vCorrectedCrossSectionI2[nCrossSection] = dInterpolationFactor*(this->estuary[nCrossSection].dGetI2(nElevationNode+1) - this->estuary[nCrossSection].dGetElevation(nElevationNode)) + this->estuary[nCrossSection].dGetElevation(nElevationNode);
        m_vCorrectedCrossSectionLeftRBLocation[nCrossSection] = dInterpolationFactor*(this->estuary[nCrossSection].dGetLeftY(nElevationNode+1) - this->estuary[nCrossSection].dGetElevation(nElevationNode)) + this->estuary[nCrossSection].dGetElevation(nElevationNode);
        m_vCorrectedCrossSectionRightRBLocation[nCrossSection] = dInterpolationFactor*(this->estuary[nCrossSection].dGetRightY(nElevationNode+1) - this->estuary[nCrossSection].dGetElevation(nElevationNode)) + this->estuary[nCrossSection].dGetElevation(nElevationNode);
    }
}


//======================================================================================================================
//! Get first node hydraulic parameters
//======================================================================================================================
void CSimulation::getFirstHydraulicParameters(const int nCrossSection, const bool bPredictedParameters) {
    if (bPredictedParameters) {
        m_vPredictedCrossSectionHydraulicRadius[nCrossSection] = this->estuary[nCrossSection].dGetHydraulicRadius(0);
        m_vPredictedCrossSectionWidth[nCrossSection] =  this->estuary[nCrossSection].dGetWidth(0);
        m_vPredictedCrossSectionElevation[nCrossSection] = this->estuary[nCrossSection].dGetElevation(0);
        m_vPredictedCrossSectionBeta[nCrossSection] = this->estuary[nCrossSection].dGetBeta(0);
        m_vPredictedCrossSectionI1[nCrossSection] = this->estuary[nCrossSection].dGetI1(0);
        m_vPredictedCrossSectionI2[nCrossSection] = this->estuary[nCrossSection].dGetI2(0);
        m_vPredictedCrossSectionLeftRBLocation[nCrossSection] = this->estuary[nCrossSection].dGetLeftY(0);
        m_vPredictedCrossSectionRightRBLocation[nCrossSection] = this->estuary[nCrossSection].dGetRightY(0);
    }
    else {
        m_vCorrectedCrossSectionHydraulicRadius[nCrossSection] = this->estuary[nCrossSection].dGetHydraulicRadius(0);
        m_vCorrectedCrossSectionWidth[nCrossSection] =  this->estuary[nCrossSection].dGetWidth(0);
        m_vCorrectedCrossSectionElevation[nCrossSection] = this->estuary[nCrossSection].dGetElevation(0);
        m_vCorrectedCrossSectionBeta[nCrossSection] = this->estuary[nCrossSection].dGetBeta(0);
        m_vCorrectedCrossSectionI1[nCrossSection] = this->estuary[nCrossSection].dGetI1(0);
        m_vCorrectedCrossSectionI2[nCrossSection] = this->estuary[nCrossSection].dGetI2(0);
        m_vCorrectedCrossSectionLeftRBLocation[nCrossSection] = this->estuary[nCrossSection].dGetLeftY(0);
        m_vCorrectedCrossSectionRightRBLocation[nCrossSection] = this->estuary[nCrossSection].dGetRightY(0);
    }
}


//======================================================================================================================
//! Get last node hydraulic parameters
//======================================================================================================================
void CSimulation::getLastHydraulicParameters(const int nCrossSection, const bool bPredictedParameters) {
    const int nLastNode = this->estuary[nCrossSection].nGetElevationSectionsNumber() - 1;
    if (bPredictedParameters) {
         m_vPredictedCrossSectionHydraulicRadius[nCrossSection] = this->estuary[nCrossSection].dGetHydraulicRadius(nLastNode);
         m_vPredictedCrossSectionWidth[nCrossSection] =  this->estuary[nCrossSection].dGetWidth(nLastNode);
         m_vPredictedCrossSectionElevation[nCrossSection] = this->estuary[nCrossSection].dGetElevation(nLastNode);
         m_vPredictedCrossSectionBeta[nCrossSection] = this->estuary[nCrossSection].dGetBeta(nLastNode);
         m_vPredictedCrossSectionI1[nCrossSection] = this->estuary[nCrossSection].dGetI1(nLastNode);
         m_vPredictedCrossSectionI2[nCrossSection] = this->estuary[nCrossSection].dGetI2(nLastNode);
         m_vPredictedCrossSectionLeftRBLocation[nCrossSection] = this->estuary[nCrossSection].dGetLeftY(nLastNode);
         m_vPredictedCrossSectionRightRBLocation[nCrossSection] = this->estuary[nCrossSection].dGetRightY(nLastNode);
    }
    else {
         m_vCorrectedCrossSectionHydraulicRadius[nCrossSection] = this->estuary[nCrossSection].dGetHydraulicRadius(nLastNode);
         m_vCorrectedCrossSectionWidth[nCrossSection] =  this->estuary[nCrossSection].dGetWidth(nLastNode);
         m_vCorrectedCrossSectionElevation[nCrossSection] = this->estuary[nCrossSection].dGetElevation(nLastNode);
         m_vCorrectedCrossSectionBeta[nCrossSection] = this->estuary[nCrossSection].dGetBeta(nLastNode);
         m_vCorrectedCrossSectionI1[nCrossSection] = this->estuary[nCrossSection].dGetI1(nLastNode);
         m_vCorrectedCrossSectionI2[nCrossSection] = this->estuary[nCrossSection].dGetI2(nLastNode);
         m_vCorrectedCrossSectionLeftRBLocation[nCrossSection] = this->estuary[nCrossSection].dGetLeftY(nLastNode);
         m_vCorrectedCrossSectionRightRBLocation[nCrossSection] = this->estuary[nCrossSection].dGetRightY(nLastNode);
    }
}


//======================================================================================================================
//! Compute the Courant number
//======================================================================================================================
void CSimulation::calculateCourantNumber() {
    double dTimestepTmp = 0.0;
    double dMinTimestep = 10000000.0;
    double dWaterDensityFactor = 1.0;
    double dX = 0.0;

    for (int i=0; i< this->m_nCrossSectionsNumber; i++) {
        if (m_vCrossSectionArea[i] != DRY_BED) {
            dX = static_cast<double>(this->estuary[i].dGetX());
            //! Mean water flow
            m_vCrossSectionU[i] = m_vCrossSectionQ[i]/m_vCrossSectionArea[i];

            //! Perturbation celerity
            m_vCrossSectionC[i] = sqrt(G*m_vCrossSectionArea[i]/m_vPredictedCrossSectionWidth[i]);

            //! TODO 003: Courant number due to perturbation on density field
            if (m_bDoWaterDensity) {
                dWaterDensityFactor = m_dCourantNumber*dX/(m_vPredictedCrossSectionWidth[i]*m_dBetaSalinityConstant);
            }

            //! Compute the timestep given the Courant number
            dTimestepTmp = m_dCourantNumber*dX/(fabs(m_vCrossSectionU[i] + m_vCrossSectionC[i]));


            if (dTimestepTmp < dMinTimestep) {
                dMinTimestep = dTimestepTmp;
            }

            if (dWaterDensityFactor < dMinTimestep) {
                dMinTimestep = dWaterDensityFactor;
            }
        }
    }
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

    double dhdx = 0.0;
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
//! Notifies the user that the simulation has ended, asks for keypress if necessary, and if compiled under GNU can send an email
//===============================================================================================================================
void CSimulation::DoSimulationEnd(int const nRtn)
{};
