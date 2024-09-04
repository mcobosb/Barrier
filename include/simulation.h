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

    //! Timestep to be saved in hours
    double m_dSimTimestep;

    //! Computational timestep obtained from Courant number
    double m_dTimestep;

    //! Lambda Value
    double m_dLambda;

    //! The initial estuarine condition, IEC [0 = in calm, 1 = water flow or 2 = elevation]
    int m_nInitialEstuarineCondition;

    //! Compute water density?
    bool m_bDoWaterDensity;

    //! The beta salinity constant [if compute water density]
    double m_dBetaSalinityConstant;

    //! The longitudinal dispersion, kh [if compute water density]
    double m_dLongitudinalDispersion;

    //! The upward estuarine condition
    int m_nUpwardEstuarineCondition;

    //! The upward estuarine condition filename
    string m_strUpwardBoundaryConditionFilename;

    //! The upward estuarine condition time vector
    vector<double> m_vUpwardBoundaryConditionTime;

    //! The upward estuarine condition value vector
    vector<double> m_vUpwardBoundaryConditionValue;

    //! The upward estuarine condition value at time t
    double m_dUpwardBoundaryValue;

    //! The downward estuarine condition
    int m_nDownwardEstuarineCondition;

    //! The downward estuarine condition filename
    string m_strDownwardBoundaryConditionFilename;

    //! The downward estuarine condition time vector
    vector<double> m_vDownwardBoundaryConditionTime;

    //! The downward estuarine condition value vector
    vector<double> m_vDownwardBoundaryConditionValue;

    //! The downward estuarine condition value at time t
    double m_dDownwardBoundaryValue;

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

    //! Number of cross-sections
    int m_nCrossSectionsNumber;

    //! Number of hydrographs
    int m_nHydrographsNumber;

    //! Vector with output variable names
    vector<string> m_vOutputVariables;

    //! Cross-section Bed slope
    vector<double> m_vCrossSectionBedSlope;

    //! Cross-section Bed slope
    vector<double> m_vCrossSectionFrictionSlope;

    //! Cross-section Manning number
    vector<double> m_vCrossSectionManningNumber;

    //! Cross-section areas
    vector<double> m_vCrossSectionArea;
    //! Predicted cross-section areas
    vector<double> m_vPredictedCrossSectionArea;
    //! Corrected cross-section areas
    vector<double> m_vCorrectedCrossSectionArea;

    //! Cross-section water flows
    vector<double> m_vCrossSectionQ;
    //! Predicted cross-section water flows
    vector<double> m_vPredictedCrossSectionQ;
    //! Corrected cross-section water flows
    vector<double> m_vCorrectedCrossSectionQ;

    //! Cross-section hydraulic radius
    vector<double> m_vCrossSectionHydraulicRadius;
    // //! Predicted cross-section hydraulic radius
    // vector<double> m_vPredictedCrossSectionHydraulicRadius;
    // //! Corrected cross-section hydraulic radius
    // vector<double> m_vCorrectedCrossSectionHydraulicRadius;

    //! Cross-section hydraulic dX
    vector<double> m_vCrossSectionDX;

    //! Cross-section hydraulic widths
    vector<double> m_vCrossSectionWidth;
    // //! Corrected cross-section widths
    // vector<double> m_vCorrectedCrossSectionWidth;

    //! Cross-section elevation
    vector<double> m_vCrossSectionElevation;
    // //! Predicted cross-section elevations
    // vector<double> m_vPredictedCrossSectionElevation;
    // //! Corrected cross-section elevations
    // vector<double> m_vCorrectedCrossSectionElevation;

    //! Cross-section betas
    vector<double> m_vCrossSectionBeta;
    // //! Corrected cross-section betas
    // vector<double> m_vCorrectedCrossSectionBeta;

    //! Cross-section I1s
    vector<double> m_vCrossSectionI1;
    // //! Corrected cross-section I1s
    // vector<double> m_vCorrectedCrossSectionI1;

    //! Cross-section I2s
    vector<double> m_vCrossSectionI2;
    // //! Corrected cross-section I2s
    // vector<double> m_vCorrectedCrossSectionI2;

    //! Cross-section water densities
    vector<double> m_vCrossSectionRho;
    // //! Corrected cross-section water densities
    // vector<double> m_vCorrectedCrossSectionRho;

    //! Cross-section left river bank locations
    vector<double> m_vCrossSectionLeftRBLocation;
    // //! Corrected cross-section left river bank locations
    // vector<double> m_vCorrectedCrossSectionLeftRBLocation;

    //! Cross-section right river bank locations
    vector<double> m_vCrossSectionRightRBLocation;
    // //! Corrected cross-section right river bank locations
    // vector<double> m_vCorrectedCrossSectionRightRBLocation;

    //! Cross-section mean water velocity
    vector<double> m_vCrossSectionU;

    //! Cross-section perturbation water velocities
    vector<double> m_vCrossSectionC;

    //! Cross-section salinities
    vector<double> m_vCrossSectionSalinity;

    //! Cross-section bottom sediment transport
    vector<double> m_vCrossSectionQb;

    //! Cross-section suspended sediment transport
    vector<double> m_vCrossSectionQs;

    //! Cross-section total sediment transport
    vector<double> m_vCrossSectionQt;

    //! gAS0 terms
    vector<double> m_vCrossSectiongAS0;

    //! gASf terms
    vector<double> m_vCrossSectiongASf;

    //! F0 terms
    vector<double> m_vCrossSectionF0;

    //! F1 terms
    vector<double> m_vCrossSectionF1;

    //! Gv0 terms
    vector<double> m_vCrossSectionGv0;

    //! Gv1 terms
    vector<double> m_vCrossSectionGv1;

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
    void initializeVectors();
    void calculateBedSlope();
    void calculateAlongEstuaryInitialConditions();
    double linearInterpolation1d(double dValue, vector<double> vX, vector<double> vY);
    void calculateHydraulicParameters();
    void interpolateHydraulicParameters(double dArea, int nCrossSection, int nElevationNode);
    void getFirstHydraulicParameters(int nCrossSection);
    void getLastHydraulicParameters(int nCrossSection);
    void calculateTimestep();
    void calculateIs();
    void calculateBoundaryConditions();
    void dryArea();
    void dryTerms();
    void doMurilloCondition();
    void calculateFrictionSlope();
    void calculateGSAterms();
    void calculateFlowTerms();
    void calculateSourceTerms();
    void calculatePredictor();
    void calculateCorrector();
    void updatePredictorBoundaries();
    void updateCorrectorBoundaries();
    void mergePredictorCorrector();
    void smoothSolution();

    //! Carries out end-of-simulation tidying (error messages etc.)
    void DoSimulationEnd(int);
};
#endif // SIMULATION_H