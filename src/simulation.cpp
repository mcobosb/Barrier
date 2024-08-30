#include <iostream>
using std::cerr;
using std::cout;
using std::endl;
using std::ios;

#include <iomanip>
using std::setprecision;

#include <string>
using std::to_string;

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
    m_dSimTimestep =
    m_dInitialConstantWaterFlow =
    m_dInitialConstantElevation = 0.0;

    m_nInitialEstuarineCondition = 0;

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
int CSimulation::nGetInitialEstuarineCondition() {
    return m_nInitialEstuarineCondition;
}

//! Method for setting the initial estuarine condition
void CSimulation::nSetInitialEstuarineCondition(int initialCondition) {
    m_nInitialEstuarineCondition = initialCondition;
};

//! Method for getting the along channel constant water flow
double CSimulation::dGetInitialConstantWaterFlow() {
    return m_dInitialConstantWaterFlow;
}

//! Method for setting the along channel constant water flow
void CSimulation::dSetInitialConstantWaterFlow(double constWaterFlow) {
    m_dInitialConstantWaterFlow = constWaterFlow;
};

//! Method for getting the along channel constant elevation
double CSimulation::dGetInitialConstantElevation() {
    return m_dInitialConstantElevation;
}

//! Method for setting the along channel constant elevation
void CSimulation::dSetInitialConstantElevation(double constElevation) {
    m_dInitialConstantElevation = constElevation;
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
int CSimulation::nGetDownwardEstuarineCondition() {
    return m_nDownwardEstuarineCondition;
}

//! Method for setting the downward estuarine condition
void CSimulation::nSetDownwardEstuarineCondition(int downwardEstuarineCondition) {
    m_nDownwardEstuarineCondition = downwardEstuarineCondition;
};

//! Method for getting the downward fix water flow
double CSimulation::dGetDownwardWaterFlow() {
    return m_dDownwardWaterFlow;
}

//! Method for setting the downward fix water flow
void CSimulation::dSetDownwardWaterFlow(double downwardWaterFlow) {
    m_dDownwardWaterFlow = downwardWaterFlow;
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
//! Appends a CHydrograph objet to the CSimulation
//===============================================================================================================================
void CSimulation::AddHydrograph(){
    // Create the
    CHydrograph hydrograph;
    hydrographs.push_back(hydrograph);
}


//! Method for getting the number of hydrographs
int CSimulation::nGetHydrographsNo(){
    return m_nHydrographsNo;
}
//! Method for setting the number of hydrographs
void CSimulation::nSetHydrographsNo(int nValue) {
    m_nHydrographsNo = nValue;
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

    calculateIs();

    return false;
}


//===============================================================================================================================
//! Compute the source terms I1 (difference between pressure thrusts applied over the frontiers x1 and x2) and
// I2(pressure forces due the channel width variation)
//===============================================================================================================================
void CSimulation::calculateIs()
{


    //Calculate I1
    int nCrossSections = estuary.size();
    for (int i = 0; i < nCrossSections; i++) {
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
    for (int i = 0; i < nCrossSections ; i++) {

        if (i == 0) {
            // Compute dx as forward difference
            dX1 = this->estuary[i+1].dGetX() - this->estuary[i].dGetX();
        }
        else if (i == nCrossSections - 1){
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
            else if (i == nCrossSections - 1) {
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


    // //Calculate I12
    // // Central finite-difference
    // dhdx, dI1dx = np.zeros([db.sizes["x"], db.sizes["z"]]), np.zeros(
    //     [db.sizes["x"], db.sizes["z"]]
    // )
    // dhdx[1:-1, :] = (db["eta"][2:, :] - db["eta"][:-2, :]) / (2 * dConfig["dx"])
    // dI1dx[1:-1, :] = (db["I1"][2:, :] - db["I1"][:-2, :]) / (2 * dConfig["dx"])
    // // Upward finite-difference
    // dhdx[0, :] = (db["eta"][1, :] - db["eta"][0, :]) / dConfig["dx"]
    // dI1dx[0, :] = (db["I1"][1, :] - db["I1"][0, :]) / dConfig["dx"]
    //
    // // Downward finite-difference
    // dhdx[-1, :] = (db["eta"][-1, :] - db["eta"][-2, :]) / dConfig["dx"]
    // dI1dx[-1, :] = (db["I1"][-1, :] - db["I1"][-2, :]) / dConfig["dx"]
    //
    // db["I2"][:] = dI1dx - db["A"][:].values * dhdx

}


//===============================================================================================================================
//! Notifies the user that the simulation has ended, asks for keypress if necessary, and if compiled under GNU can send an email
//===============================================================================================================================
void CSimulation::DoSimulationEnd(int const nRtn)
{};
