/*!
*
 * \class CSimulation
 * \brief This class runs Saint Venant simulations
 * \details TODO 001 This is a more detailed description of the CSimulation class
 * \author Manuel Cobos Budia

 * \date 2024
 * \copyright GNU General Public License
 *
 * \file simulation.h
 * \brief Contains CSimulation definitions
 *
 */

#ifndef SIMULATION_H
#define SIMULATION_H
/*===============================================================================================================================

===============================================================================================================================*/
#include <ctime>
using std::localtime;
using std::time;
using std::time_t;

#include <string>
using std::string;

#include <vector>
using std::vector;

#include <fstream>
using std::ofstream;

// #include "data_reader.h"
#include "hydrograph.h"
// #include "estuary.h"
#include "cross_section.h"

class CSimulation
{
    // friend class CHydrograph;
    friend class CHydrograph;

    // // friend class CEstuary;
    // friend class CEstuary;

    // friend class CCrossSection;
    friend class CCrossSection;

    public:
    ofstream LogStream;

    //! Duration of simulation, in hours
    double m_dSimDuration;

    //! The length of an iteration (a timestep) in hours
    double m_dSimTimestep;

    //! The initial estuarine condition, IEC [0: filename, 1: constant water flow and 2:constant elevation]
    int m_nInitialEstuarineCondition;

    //! The constant initial along channel water flow [if IEC = 1]
    double m_dInitialConstantWaterFlow;

    //! The constant initial along channel elevation [if IEC = 2]
    double m_dInitialConstantElevation;

    //! Compute water density?
    bool m_bDoWaterDensity;

    //! The beta salinity constant [if compute water density]
    double m_dBetaSalinityConstant;

    //! The longitudinal dispersion, kh [if compute water density]
    double m_dLongitudinalDispersion;

    //! The upward estuarine condition
    int m_nUpwardEstuarineCondition;

    //! The downward estuarine condition
    int m_nDownwardEstuarineCondition;

    //! The downward fix water flow
    double m_dDownwardWaterFlow;

    //! The Courant Number
    double m_dCourantNumber;

    //! Do MackComarck Limiter Flux?
    bool m_bDoMackComarckLimiterFlux;

    //! The equation for MackComarck Limiter Flux
    int m_nEquationMacComarckLimiterFlux;

    //! Psi Formula
    int m_nPsiFormula;

    //! Delta Value
    double m_dDeltaValue;

    //! Do Surface Gradient Method?
    bool m_bDoSurfaceGradientMethod;

    //! Do Source Term Balance?
    bool m_bDoSourceTermBalance;

    //! Do beta coefficient?
    bool m_bDoBetaCoefficient;

    //! Do Dry Bed?
    bool m_bDoDryBed;

    //! Do Murillo condition?
    bool m_bDoMurilloCondition;

    //! Number of hydrographs
    int m_nHydrographsNo;

    //! Names of output variables
    vector<string> m_vOutputVariables;

    //! A vector with cross-sections objects along the estuary
    vector<CCrossSection> estuary;

    void AddCrossSection();


    CSimulation();
    ~CSimulation();



    //! Method for getting the simulation duration
    [[nodiscard]] double dGetSimulationDuration();
    //! Method for setting the simulation duration
    void dSetSimulationDuration(double simDuration);

    //! Method for getting the simulation timestep
    [[nodiscard]] double dGetSimulationTimestep();
    //! Method for setting the simulation timestep
    void dSetSimulationTimestep(double simTimestep);

    //! Method for getting the initial estuarine condition
    [[nodiscard]] int nGetInitialEstuarineCondition();
    //! Method for setting the initial estuarine condition
    void nSetInitialEstuarineCondition(int initialCondition);

    //! Method for getting the along channel constant water flow
    [[nodiscard]] double dGetInitialConstantWaterFlow();
    //! Method for setting the along channel constant water flow
    void dSetInitialConstantWaterFlow(double constWaterFlow);

    //! Method for getting the along channel constant elevation
    [[nodiscard]] double dGetInitialConstantElevation();
    //! Method for setting the along channel constant elevation
    void dSetInitialConstantElevation(double constElevation);

    //! Method for getting the compute water density
    [[nodiscard]] bool bGetDoWaterDensity();
    //! Method for setting the compute water density
    void bSetDoWaterDensity(bool doWaterDensity);

    //! Method for getting the beta salinity constant
    [[nodiscard]] double dGetBetaSalinityConstant();
    //! Method for setting the beta salinity constant
    void dSetBetaSalinityConstant(double salinityConstant);

    //! Method for getting the along channel constant elevation
    [[nodiscard]] double dGetLongitudinalDispersionConstant();
    //! Method for setting the along channel constant elevation
    void dSetLongitudinalDispersionConstant(double longitudinalDispersion);

    //! Method for getting the upward estuarine condition
    [[nodiscard]] int nGetUpwardEstuarineCondition();
    //! Method for setting the upward estuarine condition
    void nSetUpwardEstuarineCondition(int upwardEstuarineCondition);

    //! Method for getting the downward estuarine condition
    [[nodiscard]] int nGetDownwardEstuarineCondition();
    //! Method for setting the downward estuarine condition
    void nSetDownwardEstuarineCondition(int downwardEstuarineCondition);

    //! Method for getting the downward water flow
    [[nodiscard]] double dGetDownwardWaterFlow();
    //! Method for setting the downward estuarine condition
    void dSetDownwardWaterFlow(double downwardWaterFlow);

    //! Method for getting the courant number
    [[nodiscard]] double dGetCourantNumber();
    //! Method for setting the courant number
    void dSetCourantNumber(double courantNumber);

    //! Method for getting if McComarck limiter flux is applied
    [[nodiscard]] bool bGetDoMcComarckLimiterFlux();
    //! Method for setting if McComarck limiter flux is applied
    void bSetDoMcComarckLimiterFlux(bool doMcComarckLimiterFlux);

    //! Method for getting equation limiter flux
    [[nodiscard]] int nGetEquationLimiterFlux();
    //! Method for setting equation limiter flux
    void nSetEquationLimiterFlux(int equationLimiterFlux);

    //! Method for getting Psi formula
    [[nodiscard]] int nGetPsiFormula();
    //! Method for setting Psi formula
    void nSetPsiFormula(int psiFormula);

    //! Method for getting Delta Value
    [[nodiscard]] double dGetDeltaValue();
    //! Method for setting Delta Value
    void dSetDeltaValue(double deltaValue);

    //! Method for getting if surface gradient method is applied
    [[nodiscard]] bool bGetDoSurfaceGradientMethod();
    //! Method for setting if surface gradient method is applied
    void bSetDoSurfaceGradientMethod(bool doSurfaceGradientMethod);

    //! Method for getting if source term balance is applied
    [[nodiscard]] bool bGetDoSurfaceTermBalance();
    //! Method for setting if source term balance is applied
    void bSetDoSurfaceTermBalance(bool doSourceTermBalance);

    //! Method for getting if beta coefficient is applied
    [[nodiscard]] bool bGetDoBetaCoefficient();
    //! Method for setting if beta coefficient is applied
    void bSetDoBetaCoefficient(bool doBetaCoefficient);

    //! Method for getting if dry bed is applied
    [[nodiscard]] bool bGetDoDryBed();
    //! Method for setting if dry bed is applied
    void bSetDoDryBed(bool doDryBed);

    //! Method for getting if Murillo condition is applied
    [[nodiscard]] bool bGetDoMurilloCondition();
    //! Method for setting if Murillo condition is applied
    void bSetDoMurilloCondition(bool doMurilloCondition);

    //! Add output variable
    void strAddOutputVariable(string strItem);

    //! A vector with hydrograph objects
    vector<CHydrograph> hydrographs;

    //! Method for getting the number of hydrographs
    [[nodiscard]] int nGetHydrographsNo();
    //! Method for setting the number of hydrographs
    void nSetHydrographsNo(int nValue);

    void AddHydrograph();



    //! Runs the simulation
    bool bDoSimulation(int, char const* []);
    void calculateIs();

    //! Carries out end-of-simulation tidying (error messages etc.)
    void DoSimulationEnd(int);
};
#endif // SIMULATION_H