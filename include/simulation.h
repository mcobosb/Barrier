/*!
*
* \class CSimulation
* \brief This class runs Saint-Venant simulations
* \details Description of CSimulation class
* \author Manuel Cobos Budia

 * \date 2026
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

#include <iomanip>  // Para setfill, setw
using std::setfill;
using std::setw;

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
public:
    //! If true, Manning is calculated as a function of water level
    bool m_bManningDependsOnLevel = false;
    //! Getter for Manning level-dependence
    [[nodiscard]] bool bGetManningDependsOnLevel() const { return m_bManningDependsOnLevel; }
    //! Setter for Manning level-dependence
    void bSetManningDependsOnLevel(bool bValue) { m_bManningDependsOnLevel = bValue; }

    // === Temperatura y balance de energía ===
    //! ¿Calcular temperatura del agua?
    bool m_bDoWaterTemperature{};

    //! Nombre del archivo de condición inicial de temperatura (si aplica)
    std::string m_strInitialTemperatureConditionFilename;

    //! Término temporal de advección-difusión de temperatura (ASt)
    std::vector<double> m_vCrossSectionTemperatureASt;

    //! Tipo de condición de frontera aguas arriba para temperatura (0=libre, 1=impuesta, 2=serie temporal)
    int m_nUpwardTemperatureCondition{};

    //! Tipo de condición de frontera aguas abajo para temperatura
    int m_nDownwardTemperatureCondition{};

    //! Archivo de condición de frontera aguas arriba (serie temporal)
    std::string m_strUpwardTemperatureBoundaryConditionFilename;

    //! Vector de tiempos de la condición de frontera aguas arriba
    std::vector<double> m_vUpwardTemperatureBoundaryConditionTime;

    //! Vector de valores de la condición de frontera aguas arriba
    std::vector<double> m_vUpwardTemperatureBoundaryConditionValue;

    //! Valor de temperatura aguas arriba (si es constante)
    double m_dUpwardTemperatureBoundaryValue{};

    //! Archivo de condición de frontera aguas abajo (serie temporal)
    std::string m_strDownwardTemperatureBoundaryConditionFilename;

    //! Vector de tiempos de la condición de frontera aguas abajo
    std::vector<double> m_vDownwardTemperatureBoundaryConditionTime;

    //! Vector de valores de la condición de frontera aguas abajo
    std::vector<double> m_vDownwardTemperatureBoundaryConditionValue;

    //! Valor de temperatura aguas abajo (si es constante)
    double m_dDownwardTemperatureBoundaryValue{};

    //! Coeficiente de difusión térmica (m²/s)
    double m_dThermalDispersion{};

    // === Forzamientos de balance de energía superficial ===


    //! Archivo único de forzamiento de balance de energía superficial (Tair, humedad relativa, viento)
    std::string m_strHeatFluxFile;

    //! Coeficiente de transferencia de calor sensible (CS)
    double m_dHeatFlux_CS{};

    //! Coeficiente de transferencia de calor latente (CL)
    double m_dHeatFlux_CL{};

    // === Series temporales de forzamiento ===

    //! Tiempos y valores de forzamiento de heat flux (Tair, humedad relativa, viento)
    std::vector<double> m_vHeatFluxTime;
    std::vector<double> m_vHeatFluxAirTemp;
    std::vector<double> m_vHeatFluxRelHumidity;
    std::vector<double> m_vHeatFluxWind;
    // === Opciones para continuar simulación ===
    bool m_bContinueSimulation = false;
    std::string m_strContinueNetcdfPath;
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

    //! Save at every computational timestep (for debugging)?
    bool m_bSaveAllTimesteps;

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
    //! The beta temperature constant [if compute water temperature]
    double m_dBetaTemperatureConstant{};

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

    //! The upward estuarine salinity condition filename
    string m_strUpwardSalinityBoundaryConditionFilename;

    //! The upward estuarine salinity condition time vector
    vector<double> m_vUpwardSalinityBoundaryConditionTime;

    //! The upward estuarine salinity condition value vector
    vector<double> m_vUpwardSalinityBoundaryConditionValue;

    //! The upward estuarine salinity condition value at time t
    double m_dUpwardSalinityBoundaryValue{};

    //! The upward estuarine salinity condition value at previous the last node and time t
    double m_dNextUpwardSalinityBoundaryValue{};

    //! The downward estuarine salinity condition filename
    string m_strDownwardSalinityBoundaryConditionFilename;

    //! The downward estuarine salinity condition time vector
    vector<double> m_vDownwardSalinityBoundaryConditionTime;

    //! The downward estuarine salinity condition value vector
    vector<double> m_vDownwardSalinityBoundaryConditionValue;

    //! The downward estuarine salinity condition value at time t
    double m_dDownwardSalinityBoundaryValue{};

    //! The downward estuarine salinity condition value at previous the last node and time t
    double m_dNextDownwardSalinityBoundaryValue{};
    
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
    
    //! Smooth bathymetry before simulation?
    bool m_bDoSmoothBathymetry{};
    int m_nBathymetrySmoothingPasses{1};
    double m_dBathymetrySmoothingAlpha{0.25};

        //! Nivel máximo de marea astronómica calculado
    double m_dMaxAstronomicalTide = 0.0;

    int m_nThresholddBdeta = 0;
    //! Vector de eta donde se supera el threshold de gradiente de B(eta) por sección
    vector<double> m_vEtaWidthGradientThreshold;

    //! Smooth solution (regularization) during simulation?
    bool m_bDoSmoothSolution{};
    int m_nSolutionSmoothingPasses{1};
    double m_dSolutionSmoothingAlpha{0.25};

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

    //! Cross-section Bed slope for predictor (forward difference) - balance de términos fuente
    vector<double> m_vCrossSectionBedSlopePredictor;

    //! Cross-section Bed slope for corrector (backward difference) - balance de términos fuente
    vector<double> m_vCrossSectionBedSlopeCorrector;

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

    //! Cross-section water depth (predicted, for corrector/baroclinic coupling)
    vector<double> m_vPredictedCrossSectionWaterDepth;

    //! Cross-section water elevation (over the mean water level)
    vector<double> m_vCrossSectionWaterElevation;

    //! Cross-section I1 pressure integral
    vector<double> m_vCrossSectionI1;
    //! Predicted cross-section I1 pressure integral
    vector<double> m_vPredictedCrossSectionI1;

    //! Cross-section I2 momentum integral
    vector<double> m_vCrossSectionI2;

    //! Cross-section dI1/dx (pressure gradient for non-prismatic channels)
    vector<double> m_vCrossSectionDI1Dx;

    //! Cross-section DhDx
    vector<double> m_vCrossSectionDhDx;

    //! Cross-section water densities (estado actual)
    vector<double> m_vCrossSectionDensity;
    //! Cross-section water densities (predicho, para corrector/baroclínico)
    vector<double> m_vPredictedCrossSectionDensity;
    //! Gradiente de densidad dRho/dx (baroclínico)
    vector<double> m_vCrossSectionDRhoDx;

    //! Cross-section left river bank locations
    vector<double> m_vCrossSectionLeftRBLocation;

    //! Cross-section right river bank locations
    vector<double> m_vCrossSectionRightRBLocation;
    
    //! Cross-section left river bank UTM X coordinates
    vector<double> m_vCrossSectionLeftRBLocation_UTM_X;
    
    //! Cross-section left river bank UTM Y coordinates
    vector<double> m_vCrossSectionLeftRBLocation_UTM_Y;
    
    //! Cross-section right river bank UTM X coordinates
    vector<double> m_vCrossSectionRightRBLocation_UTM_X;
    
    //! Cross-section right river bank UTM Y coordinates
    vector<double> m_vCrossSectionRightRBLocation_UTM_Y;

    //! Cross-section mean water velocity
    vector<double> m_vCrossSectionU;

    //! Cross-section perturbation water velocities
    vector<double> m_vCrossSectionC;

    //! Vector de salinidad por sección transversal (g/kg o PSU)
    vector<double> m_vCrossSectionSalinity;

    //! Vector de temperatura por sección transversal (°C)
    vector<double> m_vCrossSectionTemperature;

    //! Predicted Cross-section salinity
    vector<double> m_vPredictedCrossSectionS;

    // Corrected Cross-section salinity
    vector<double> m_vCorrectedCrossSectionS;
    
    // Predicted Cross-section temperature
    vector<double> m_vPredictedCrossSectionT;

    // Corrected Cross-section temperature
    vector<double> m_vCorrectedCrossSectionT;
    
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

    // ⚡ TVD limiter work arrays (pre-allocated to avoid per-timestep allocations)
    vector<double> m_vTVD_a1_med;
    vector<double> m_vTVD_a2_med;
    vector<double> m_vTVD_alfa1_med;
    vector<double> m_vTVD_alfa2_med;
    vector<double> m_vTVD_psi1_med;
    vector<double> m_vTVD_psi2_med;
    vector<double> m_vTVD_r1_med;
    vector<double> m_vTVD_r2_med;
    vector<double> m_vTVD_fi1_med;
    vector<double> m_vTVD_fi2_med;
    vector<double> m_vTVD_Factor1;
    vector<double> m_vTVD_Factor2;

    // ⚡ Salinity gradient work arrays (pre-allocated)
    vector<double> m_vSalinity_KAS_forward;
    vector<double> m_vSalinity_KAS_backward;
    vector<double> m_vSalinity_AUS_diff;

    // ⚡ Precalculated constants (computed once at initialization)
    vector<double> m_vManningNumberSquared;      // Manning² (usado en fricción ~4000 veces/día)
    vector<double> m_vInvDX;                     // 1/ΔX (usado en gradientes ~8000 veces/día)
    vector<double> m_vDxSum;                     // ΔX[i] + ΔX[i+1] (diferencias centradas)
    vector<double> m_vInvDxSum;                  // 1/(ΔX[i] + ΔX[i+1])
    vector<double> m_vGtimesDX;                  // g*ΔX (término constante)

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

    //! Method for getting the compute water temperature
    [[nodiscard]] bool bGetDoWaterTemperature() const;
    //! Method for setting the compute water temperature
    void bSetDoWaterTemperature(bool doWaterTemperature);

    //! Method for getting the compute water salinity
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
    //! Method for getting the beta temperature constant
    [[nodiscard]] double dGetBetaTemperatureConstant() const { return m_dBetaTemperatureConstant; }
    //! Method for setting the beta temperature constant
    void dSetBetaTemperatureConstant(double tempConstant) { m_dBetaTemperatureConstant = tempConstant; }

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
    
    //! Method for getting smooth bathymetry flag
    [[nodiscard]] bool bGetDoSmoothBathymetry() const;
    //! Method for setting smooth bathymetry flag
    void bSetDoSmoothBathymetry(bool doSmoothBathymetry);
    
    //! Method for getting smooth solution flag
    [[nodiscard]] bool bGetDoSmoothSolution() const;
    //! Method for setting smooth solution flag
    void bSetDoSmoothSolution(bool doSmoothSolution);

    //! Method for getting save all timesteps flag
    [[nodiscard]] bool bGetSaveAllTimesteps() const;
    //! Method for setting save all timesteps flag
    void bSetSaveAllTimesteps(bool saveAllTimesteps);

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
    string generateOutputFileName() const;
    static double linearInterpolation1d(double dValue, const vector<double> &vX, const vector<double> &vY);
    void calculateHydraulicParameters();
    void calculateRiverBankUTMCoordinates();
    void interpolateHydraulicParameters(double dArea, int nCrossSection, int nElevationNode);
    void getFirstHydraulicParameters(int nCrossSection);
    void getLastHydraulicParameters(int nCrossSection);
    void calculateTimestep();
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
    void mergeTracerPredictorCorrector();
    void smoothSolution();
    void smoothBathymetry();
    // Predictor-corrector para salinidad y temperatura
    void calculate_salinity(); // Wrapper legacy (llama a predictor/corrector según m_nPredictor)
    void calculate_salinity_predictor();
    void calculate_salinity_corrector();
    void calculate_salinity_gradient();

    void calculateRadiativeFluxes();
    void calculate_temperature(); // Wrapper legacy (llama a predictor/corrector según m_nPredictor)
    void calculate_temperature_predictor();
    void calculate_temperature_corrector();

    static double n_eta(double eta, double maxTide, double etaMaxGrad);
    static double getMaxAstronomicalTide(const std::string& tidesFile);

    void AnnounceProgress();

    //! Carries out end-of-simulation tidying (error messages etc.)
    void bDoSimulationEnd();

    private:
    // ✅ AÑADIR: Propiedades de fecha (que vendrán de CDataReader)
    int m_nSimStartSec;
    int m_nSimStartMin;
    int m_nSimStartHour;
    int m_nSimStartDay;
    int m_nSimStartMonth;
    int m_nSimStartYear;

    public:
    // ✅ AÑADIR: Getters para fecha de inicio
    int nGetSimStartSec() const { return m_nSimStartSec; }
    int nGetSimStartMin() const { return m_nSimStartMin; }
    int nGetSimStartHour() const { return m_nSimStartHour; }
    int nGetSimStartDay() const { return m_nSimStartDay; }
    int nGetSimStartMonth() const { return m_nSimStartMonth; }
    int nGetSimStartYear() const { return m_nSimStartYear; }

    // ✅ AÑADIR: Setters para fecha de inicio
    void nSetSimStartSec(int sec) { m_nSimStartSec = sec; }
    void nSetSimStartMin(int min) { m_nSimStartMin = min; }
    void nSetSimStartHour(int hour) { m_nSimStartHour = hour; }
    void nSetSimStartDay(int day) { m_nSimStartDay = day; }
    void nSetSimStartMonth(int month) { m_nSimStartMonth = month; }
    void nSetSimStartYear(int year) { m_nSimStartYear = year; }

    //! Set complete simulation start date and time at once
    void setSimulationStartDateTime(int year, int month, int day, int hour, int min, int sec) {
        m_nSimStartYear = year;
        m_nSimStartMonth = month;
        m_nSimStartDay = day;
        m_nSimStartHour = hour;
        m_nSimStartMin = min;
        m_nSimStartSec = sec;
    }

    //! Get simulation start date and time as formatted string (ISO 8601)
    string getSimulationStartDateTimeString() const {
        ostringstream oss;
        oss << setfill('0') << setw(4) << m_nSimStartYear << "-"
            << setw(2) << m_nSimStartMonth << "-"
            << setw(2) << m_nSimStartDay << " "
            << setw(2) << m_nSimStartHour << ":"
            << setw(2) << m_nSimStartMin << ":"
            << setw(2) << m_nSimStartSec;
        return oss.str();
    }

    private:
    // ✅ CORREGIR: Estas deben ser vector<vector<double>>
    vector<double> m_vBedZ;
    vector<double> m_vManningN;
    vector<double> m_vPositionX;
    vector<double> m_vBeta;  // ✅ ELIMINADO de vector<double> a vector<vector<double>>
    
    // ✅ CAMBIAR de vector<double> a vector<vector<double>>
    vector<vector<double>> m_vWidth;      // ✅ CORREGIDO
    vector<vector<double>> m_vLeftY;      // ✅ CORREGIDO
    vector<vector<double>> m_vRightY;     // ✅ CORREGIDO
    
    vector<vector<double>> m_vEstuaryAreas;
    vector<vector<double>> m_vEstuaryHydraulicRadius;
    vector<vector<double>> m_vEstuaryWaterDepths;
    vector<vector<double>> m_vEstuaryI1;  // Pressure integral I1 for each elevation
    vector<vector<double>> m_vPrecalculatedSecondTerm;
    vector<int> m_vElevationSectionsCount;

    public:
    // ✅ Solo método de optimización simple
    void precomputeEstuaryData();
};
#endif // SIMULATION_H