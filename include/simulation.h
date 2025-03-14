/*!
*
 * \class CSimulation
 * \brief This class runs Saint-Venant simulations
 * \details Description of CSimulation class
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

#include "hydrograph.h"
#include "cross_section.h"
#include "data_writer.h"
#include "data_reader.h"
#include "screen_presenter.h"

class CSimulation
{
    // friend class CHydrograph;
    friend class CHydrograph;

    // friend class CDataWriter;
    friend class CDataWriter;

    // friend class CCrossSection;
    friend class CCrossSection;

    // friend class CDataReader;
    friend class CDataReader;

    // friend class CScreenPresenter;
    friend class CScreenPresenter;


    public:
    //! System start-simulation time
    time_t m_tSysStartTime{};

    //! Sytem loop-start-simulation time
    time_t m_tSysStartLoopTime{};

    ofstream LogStream;

    //! Detail of the log file
    int m_nLogFileDetail{};

    //! Output Log time id
    int m_nTimeLogId{};

    //! Duration of simulation, in seconds
    double m_dSimDuration;

    //! Timestep to be saved in hours
    double m_dSimTimestep;

    //! Current time in seconds
    double m_dCurrentTime{};

    //! Factor for updating the time
    double m_dTimeFactor;

    //! Vector of output times
    vector<double> m_vOutputTimes;

    //! Vector with the ids of output times
    vector<int> m_vOutputTimesIds;

    //! Output timestep id
    int m_nTimeId;

    //! Name of main output file
    string m_strOutFile;

    //! Name of main log file
    string m_strLogFile;

    //! Is this timestep saved?
    bool m_bSaveTime;

    //! Computational timestep obtained from Courant number
    double m_dTimestep;

    //! Lambda Value
    double m_dLambda;

    //! The initial estuarine condition, IEC [0 = in calm, 1 = water flow or 2 = elevation]
    int m_nInitialEstuarineCondition;

    //! Compute sediment transport?
    bool m_bDoSedimentTransport;

    //! The equation for Sediment Transport
    int m_nEquationSedimentTransport;

    //! Compute water salinity?
    bool m_bDoWaterSalinity{};

    //! Name of the initial salinity condition file [if compute water salinity]
    string m_strInitialSalinityConditionFilename;

    //! Upward salinity condition
    int m_nUpwardSalinityCondition;

    //! Downward salinity condition
    int m_nDownwardSalinityCondition;

    //! The beta salinity constant [if compute water salinity]
    double m_dBetaSalinityConstant{};

    //! The longitudinal dispersion, kh [if compute water salinity]
    double m_dLongitudinalDispersion{};

    //! The upward estuarine condition
    int m_nUpwardEstuarineCondition{};

    //! The upward estuarine condition filename
    string m_strUpwardBoundaryConditionFilename;

    //! The upward estuarine condition time vector
    vector<double> m_vUpwardBoundaryConditionTime;

    //! The upward estuarine condition value vector
    vector<double> m_vUpwardBoundaryConditionValue;

    //! The upward estuarine condition value at time t
    double m_dUpwardBoundaryValue{};

    //! The upward estuarine condition value at second node and time t
    double m_dNextUpwardBoundaryValue{};

    //! The downward estuarine condition
    int m_nDownwardEstuarineCondition{};

    //! The downward estuarine condition filename
    string m_strDownwardBoundaryConditionFilename;

    //! The downward estuarine condition time vector
    vector<double> m_vDownwardBoundaryConditionTime;

    //! The downward estuarine condition value vector
    vector<double> m_vDownwardBoundaryConditionValue;

    //! The downward estuarine condition value at time t
    double m_dDownwardBoundaryValue{};

    //! The downward estuarine condition value at previous the last node and time t
    double m_dNextDownwardBoundaryValue{};

    //! The Courant Number
    double m_dCourantNumber{};

    //! Do McCormack Limiter Flux?
    bool m_bDoMcCormackLimiterFlux{};

    //! The equation for McCormack Limiter Flux
    int m_nEquationMcCormackLimiterFlux{};

    //! Psi Formula
    int m_nPsiFormula{};

    //! Delta Value
    double m_dDeltaValue{};

    //! Do Surface Gradient Method?
    bool m_bDoSurfaceGradientMethod{};

    //! Do Source Term Balance?
    bool m_bDoSourceTermBalance{};

    //! Do beta coefficient?
    bool m_bDoBetaCoefficient{};

    //! Do Dry Bed?
    bool m_bDoDryBed{};

    //! Do Murillo condition?
    bool m_bDoMurilloCondition{};

    //! Compute water density?
    bool m_bDoWaterDensity{};

    //! Number of cross-sections
    int m_nCrossSectionsNumber{};

    //! Number of hydrographs
    int m_nHydrographsNumber{};

    //! Predictor phase? (0: initial calculation; 1: yes, 2: correction phase
    int m_nPredictor{};

    //! Vector with output variable names
    vector<string> m_vOutputVariables;

    //! Cross-section Location X
    vector<double> m_vCrossSectionX;

    //! Cross-section Bed slope
    vector<double> m_vCrossSectionBedSlope;

    //! Cross-section Bed slope Direction +/- ve, 1/-1
    vector<int> m_vCrossSectionBedSlopeDirection;

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

    //! Cross-section dX
    vector<double> m_vCrossSectionDX;

    //! Cross-section widths
    vector<double> m_vCrossSectionWidth;

    //! Cross-section water depth
    vector<double> m_vCrossSectionWaterDepth;

    //! Cross-section water elevation (over the mean water level)
    vector<double> m_vCrossSectionWaterElevation;

    //! Cross-section betas
    vector<double> m_vCrossSectionBeta;

    //! Cross-section DhDx
    vector<double> m_vCrossSectionDhDx;

    //! Cross-section water densities
    vector<double> m_vCrossSectionRho;

    //! Cross-section left river bank locations
    vector<double> m_vCrossSectionLeftRBLocation;

    //! Cross-section right river bank locations
    vector<double> m_vCrossSectionRightRBLocation;

    //! Cross-section mean water velocity
    vector<double> m_vCrossSectionU;

    //! Cross-section perturbation water velocities
    vector<double> m_vCrossSectionC;

    //! Cross-section salinity
    vector<double> m_vCrossSectionSalinity;

    //! Cross-section salinity temporal gradient (ASt term)
    vector<double> m_vCrossSectionSalinityASt;

    //! Cross-section bottom sediment transport
    vector<double> m_vCrossSectionQb;

    //! Cross-section suspended sediment transport
    vector<double> m_vCrossSectionQs;

    //! Cross-section total sediment transport
    vector<double> m_vCrossSectionQt;

    //! gAS0 terms
    vector<double> m_vCrossSection_gAS0;

    //! gASf terms
    vector<double> m_vCrossSection_gASf;

    //! F0 terms
    vector<double> m_vCrossSectionF0;

    //! F1 terms
    vector<double> m_vCrossSectionF1;

    //! Gv0 terms
    vector<double> m_vCrossSectionGv0;

    //! Gv1 terms
    vector<double> m_vCrossSectionGv1;

    //! D1 terms
    vector<double> m_vCrossSectionD1Factor;

    //! D2 terms
    vector<double> m_vCrossSectionD2Factor;

    //! Murillo Factor vector
    vector<double> m_vCrossSectionMurilloFactor;

    //! D50
    vector<double> m_vCrossSectionD50;

    //! D90
    vector<double> m_vCrossSectionD90;

    //! DiamX
    vector<double> m_vCrossSectionDiamX;

    //! Daveraged
    vector<double> m_vCrossSectionDaveraged;

    //! Sediment sigma
    vector<double> m_vCrossSectionSedimentSigma;

    //! Parameter of relative density
    vector<double> m_vCrossSectionRhos;

    //! Thickness of sediment layer
    vector<double> m_vCrossSectionThickness;

    //! Lateral sources from hydro at current time
    vector<double> m_vLateralSourcesAtT;

    //! Code for error handling
    int m_nStringError;

    //! Boolean for error
    bool m_bReturnError;

    //! Text attachment for error handling
    string m_strErrorAttachment;

    //! A vector with cross-sections objects along the estuary
    vector<CCrossSection> estuary;

    CScreenPresenter presenter;
    CDataReader reader;
    CDataWriter writer;

    void AddCrossSection();


    CSimulation();
    ~CSimulation();



    //! Method for getting the simulation duration
    [[nodiscard]] double dGetSimulationDuration() const;
    //! Method for setting the simulation duration
    void dSetSimulationDuration(double simDuration);

    //! Method for getting the simulation timestep
    [[nodiscard]] double dGetSimulationTimestep() const;
    //! Method for setting the simulation timestep
    void dSetSimulationTimestep(double simTimestep);

    //! Method for getting the initial estuarine condition
    [[nodiscard]] int nGetInitialEstuarineCondition() const;
    //! Method for setting the initial estuarine condition
    void nSetInitialEstuarineCondition(int initialCondition);

    //! Method for getting the compute sediment transport
    [[nodiscard]] bool bGetDoSedimentTransport() const;
    //! Method for setting the compute sediment transport
    void bSetDoSedimentTransport(bool doSedimentTransport);

    //! Method for getting equation of sediment transport
    [[nodiscard]] int nGetEquationSedimentTransport() const;
    //! Method for setting equation of sediment transport
    void nSetEquationSedimentTransport(int equationSedimentTransport);

    //! Method for getting equation of sediment transport
    [[nodiscard]] double dGetSedimentDensity() const;
    //! Method for setting equation of sediment transport
    void dSetSedimentDensity(double sedimentDensity);

    //! Method for getting the compute water  salinity
    [[nodiscard]] bool bGetDoWaterSalinity() const;
    //! Method for setting the compute water salinity
    void bSetDoWaterSalinity(bool doWaterSalinity);

    //! Method for getting the compute water salinity
    [[nodiscard]] int nGetUpwardSalinityCondition() const;
    //! Method for setting the compute water salinity
    void nSetUpwardSalinityCondition(int nUpwardCondition);

    //! Method for getting the compute water salinity
    [[nodiscard]] int nGetDownwardSalinityCondition() const;
    //! Method for setting the compute water salinity
    void nSetDownwardSalinityCondition(int nDownwardCondition);

    //! Method for getting the beta salinity constant
    [[nodiscard]] double dGetBetaSalinityConstant() const;
    //! Method for setting the beta salinity constant
    void dSetBetaSalinityConstant(double salinityConstant);

    //! Method for getting the along channel constant elevation
    [[nodiscard]] double dGetLongitudinalDispersionConstant() const;
    //! Method for setting the along channel constant elevation
    void dSetLongitudinalDispersionConstant(double longitudinalDispersion);

    //! Method for getting the upward estuarine condition
    [[nodiscard]] int nGetUpwardEstuarineCondition() const;
    //! Method for setting the upward estuarine condition
    void nSetUpwardEstuarineCondition(int upwardEstuarineCondition);

    //! Method for getting the downward estuarine condition
    [[nodiscard]] int nGetDownwardEstuarineCondition() const;
    //! Method for setting the downward estuarine condition
    void nSetDownwardEstuarineCondition(int downwardEstuarineCondition);

    //! Method for getting the courant number
    [[nodiscard]] double dGetCourantNumber() const;
    //! Method for setting the courant number
    void dSetCourantNumber(double courantNumber);

    //! Method for getting if McComarck limiter flux is applied
    [[nodiscard]] bool bGetDoMcComarckLimiterFlux() const;
    //! Method for setting if McComarck limiter flux is applied
    void bSetDoMcComarckLimiterFlux(bool doMcComarckLimiterFlux);

    //! Method for getting equation limiter flux
    [[nodiscard]] int nGetEquationLimiterFlux() const;
    //! Method for setting equation limiter flux
    void nSetEquationLimiterFlux(int equationLimiterFlux);

    //! Method for getting Psi formula
    [[nodiscard]] int nGetPsiFormula() const;
    //! Method for setting Psi formula
    void nSetPsiFormula(int psiFormula);

    //! Method for getting Delta Value
    [[nodiscard]] double dGetDeltaValue() const;
    //! Method for setting Delta Value
    void dSetDeltaValue(double deltaValue);

    //! Method for getting if surface gradient method is applied
    [[nodiscard]] bool bGetDoSurfaceGradientMethod() const;
    //! Method for setting if surface gradient method is applied
    void bSetDoSurfaceGradientMethod(bool doSurfaceGradientMethod);

    //! Method for getting if source term balance is applied
    [[nodiscard]] bool bGetDoSurfaceTermBalance() const;
    //! Method for setting if source term balance is applied
    void bSetDoSurfaceTermBalance(bool doSourceTermBalance);

    //! Method for getting if beta coefficient is applied
    [[nodiscard]] bool bGetDoBetaCoefficient() const;
    //! Method for setting if beta coefficient is applied
    void bSetDoBetaCoefficient(bool doBetaCoefficient);

    //! Method for getting if dry bed is applied
    [[nodiscard]] bool bGetDoDryBed() const;
    //! Method for setting if dry bed is applied
    void bSetDoDryBed(bool doDryBed);

    //! Method for getting if Murillo condition is applied
    [[nodiscard]] bool bGetDoMurilloCondition() const;
    //! Method for setting if Murillo condition is applied
    void bSetDoMurilloCondition(bool doMurilloCondition);

    //! Method for getting the compute water density
    [[nodiscard]] bool bGetDoWaterDensity() const;
    //! Method for setting the compute water density
    void bSetDoWaterDensity(bool doWaterDensity);

    //! Add output variable
    void strAddOutputVariable(const string& strItem);

    //! Read hydrograph input?
    bool m_bHydroFile;

    //! A vector with hydrograph objects
    vector<CHydrograph> hydrographs;

    //! Method for getting the number of hydrographs
    [[nodiscard]] int nGetHydrographsNumber() const;
    //! Method for setting the number of hydrographs
    void nSetHydrographsNumber(int nValue);

    //! Get the vector of a variable
    vector<double> vGetVariable(const string& strVariableName) const;

    void AddHydrograph();



    //! Runs the simulation
    void bDoSimulation(int, char const* []);
    void initializeVectors();
    void calculateBedSlope();
    void calculateAlongEstuaryInitialConditions();
    static double linearInterpolation1d(double dValue, const vector<double> &vX, const vector<double> &vY);
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
    void calculate_sediment_transport();
    void calculate_density();
    void calculate_GS_A_terms();
    void calculateFlowTerms();
    void calculateSourceTerms();
    void calculatePredictor();
    void calculateCorrector();
    void updatePredictorBoundaries();
    void updateCorrectorBoundaries();
    void updateBoundaries();
    void mergePredictorCorrector();
    void smoothSolution();
    void calculate_salinity_gradient();
    void calculate_salinity();

    void AnnounceProgress();

    //! Carries out end-of-simulation tidying (error messages etc.)
    void bDoSimulationEnd();
};
#endif // SIMULATION_H