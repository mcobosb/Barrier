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

#include <cmath>    // Para std::cos, std::pow en funciones inline

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

class CSimulation {
public:
    // === PHYSICAL CONSTANTS ===
    // Mathematical constants
    static constexpr double PI = 3.14159265358979323846;
    static constexpr double DEG_TO_RAD = PI / 180.0;
    static constexpr double RAD_TO_DEG = 180.0 / PI;
    
    // Water properties
    static constexpr double WATER_DENSITY = 1000.0;           // Water density (kg/m³)
    static constexpr double WATER_SPECIFIC_HEAT = 4186.0;     // Specific heat capacity (J/(kg·°C))
    static constexpr double WATER_EMISSIVITY = 0.97;          // Water surface emissivity
    
    // Air properties (at standard conditions)
    static constexpr double AIR_SPECIFIC_HEAT = 1005.0;       // Air specific heat at constant pressure (J/(kg·K))
    static constexpr double AIR_GAS_CONSTANT = 287.05;        // Specific gas constant for dry air (J/(kg·K))
    static constexpr double ATM_PRESSURE_DEFAULT = 101325.0;  // Standard atmospheric pressure (Pa)
    
    // Radiation constants
    static constexpr double STEFAN_BOLTZMANN = 5.67e-8;       // Stefan-Boltzmann constant (W/(m²·K⁴))
    static constexpr double SOLAR_CONSTANT = 1367.0;          // Solar constant at top of atmosphere (W/m²)
    static constexpr double ATM_TRANSMISSIVITY = 0.75;        // Clear-sky atmospheric transmissivity
    
    // === INLINE HELPER FUNCTIONS (Efficient for small computations) ===
    
    //! Calculate latent heat of vaporization (J/kg) from water temperature (°C)
    //! Empirical formula: L_v = (2501 - 2.37*T) * 1000 J/kg
    static inline double calc_Lv(double T_water) {
        return (2501.0 - 2.37 * T_water) * 1000.0;
    }

    //! Calculate air density (kg/m³) from temperature (°C) and pressure (Pa)
    //! Uses ideal gas law: ρ = P / (R * T)
    static inline double calc_rho_air(double T_air, double Pressure_Pa = ATM_PRESSURE_DEFAULT) {
        return Pressure_Pa / (AIR_GAS_CONSTANT * (T_air + 273.15));
    }

    //! Calculate water albedo using Briegleb et al. (1986) model
    //! Separates direct and diffuse components with cloud-dependent weighting
    //! zenith_rad: solar zenith angle (0 = overhead, π/2 = horizon)
    //! cloud_cover: cloud fraction (0.0 = clear, 1.0 = overcast)
    static inline double calc_albedo_briegleb(double zenith_rad, double cloud_cover = 0.0) {
        double cos_theta = std::cos(zenith_rad);
        
        // Night or sun below horizon
        if (cos_theta <= 0.0) return 0.06;
        
        // Direct beam albedo (Briegleb et al. 1986 formulation)
        double alpha_dir = 0.026 / (std::pow(cos_theta, 1.7) + 0.065) + 
                           0.15 * (cos_theta - 0.1) * (cos_theta - 0.5) * (cos_theta - 1.0);
        
        // Diffuse albedo (constant for water surfaces)
        double alpha_dif = 0.06;
        
        // Fraction of direct vs diffuse radiation (cloud-dependent)
        // Clear sky: mostly direct; overcast: mostly diffuse
        double f_direct = 1.0 - 0.65 * cloud_cover * cloud_cover;
        double f_diffuse = 1.0 - f_direct;
        
        // Total albedo (weighted average of direct and diffuse components)
        double alpha_total = f_direct * alpha_dir + f_diffuse * alpha_dif;
        
        // Physical bounds for water albedo
        return std::max(0.03, std::min(alpha_total, 0.40));
    }

    //! Calculate relative humidity from current air temperature and daily minimum temperature
    //! Uses FAO-56 method: assumes T_min approximates dew point temperature
    //! T_air_current: Current air temperature (°C)
    //! T_min_daily: Daily minimum air temperature (°C)
    //! Returns: Relative humidity (0-100%)
    static inline double calc_rh_from_temp(double T_air_current, double T_min_daily) {
        // 1. Saturation vapor pressure (es) at current temperature using Magnus-Tetens formula
        double es = 0.6108 * std::exp((17.27 * T_air_current) / (T_air_current + 237.3));
        
        // 2. Actual vapor pressure (ea) using T_min as proxy for dew point
        // Assumption: T_dew ≈ T_min (air saturates at coolest moment of the day)
        double ea = 0.6108 * std::exp((17.27 * T_min_daily) / (T_min_daily + 237.3));
        
        // 3. Calculate relative humidity
        double rh = 100.0 * (ea / es);
        
        // Limit to 100% for numerical safety
        return std::min(rh, 100.0);
    }



    // === PHYSICAL PROCESS FLAGS ===
    //! Calculate water temperature?
    bool m_bDoWaterTemperature{false};

    enum class ETemperatureMode {
        Off = 0,
        Transport = 1,
        Given = 2
    };

    //! Temperature mode: off / solve transport / prescribed uniform time series
    ETemperatureMode m_eTemperatureMode{ETemperatureMode::Off};

    //! If mode=Given: file with (time_since_start_s, temperature_C)
    std::string m_strGivenTemperatureFilename;
    std::vector<double> m_vGivenTemperatureTime;
    std::vector<double> m_vGivenTemperatureValue;

    bool bIsGivenTemperature() const { return m_eTemperatureMode == ETemperatureMode::Given; }
    double getGivenTemperatureAtTime(double t) const;
    void applyGivenTemperatureAtTime(double t);
    
    //! Compute water salinity?
    bool m_bDoWaterSalinity{false};
    
    //! Compute water density?
    bool m_bDoWaterDensity{false};
    
    //! If true, Manning is calculated as a function of water level
    bool m_bManningDependsOnLevel{false};
    
    //! Compute sediment transport?
    bool m_bDoSedimentTransport{false};
    
    //! Do Dry Bed?
    bool m_bDoDryBed{false};
    
    //! Do McCormack Limiter Flux?
    bool m_bDoMcCormackLimiterFlux{false};
    
    //! Do Surface Gradient Method?
    bool m_bDoSurfaceGradientMethod{false};
    
    //! Do Source Term Balance?
    bool m_bDoSourceTermBalance{false};
    
    //! Do beta coefficient?
    bool m_bDoBetaCoefficient{false};
    
    //! Do Murillo condition?
    bool m_bDoMurilloCondition{false};
    
    //! Smooth bathymetry before simulation?
    bool m_bDoSmoothBathymetry{false};
    
    //! Smooth solution (regularization) during simulation?
    bool m_bDoSmoothSolution{false};

    //! Auto-smoothing trigger threshold computed from upstream discharge time series (|Q|, m3/s).
    //! 0 means "not available" (e.g., upstream BC is not type 3 or no time series loaded).
    double m_dAutoSmoothAbsQThreshold{0.0};

    //! Per-timestep marker: smoothing was auto-enabled this timestep (for log "*")
    bool m_bAutoSmoothAppliedThisStep{false};
    
    //! Save at every computational timestep (for debugging)?
    bool m_bSaveAllTimesteps{false};
    
    //! Flag to calculate RH from Tair and Tmin (when RH time series is not available)
    bool m_bCalculateRHFromTemperature{false};
    
    //! Simulation continuation flag
    bool m_bContinueSimulation{false};
    
    //! Read hydrograph input?
    bool m_bHydroFile{false};
    
    //! Boolean for error
    bool m_bReturnError{false};
    
    //! Is this timestep saved?
    bool m_bSaveTime{false};
    
    //! Flag for extreme flow conditions
    bool m_bExtremeFlowConditions{false};
    
    //! Use independent transport limiter (different from hydrodynamics)?
    bool m_bUseIndependentTransportLimiter{false};
    
    // === CACHE AND OPTIMIZATION ===
    //! Cache for binary search in calculateHydraulicParameters() - exploits spatial coherence
    //! Stores last interpolation index for each cross-section (speeds up ~30%)
    mutable std::vector<int> m_vLastInterpolationIndex;

    //! Initial temperature condition filename (if applicable)
    std::string m_strInitialTemperatureConditionFilename;

    //! Temperature advection-diffusion temporal term (ASt)
    std::vector<double> m_vCrossSectionTemperatureASt;

    //! Upstream temperature boundary condition type (0=free, 1=imposed, 2=time series)
    int m_nUpwardTemperatureCondition{};

    //! Downstream temperature boundary condition type
    int m_nDownwardTemperatureCondition{};

    //! Upstream boundary condition file (time series)
    std::string m_strUpwardTemperatureBoundaryConditionFilename;

    //! Upstream boundary condition time vector
    std::vector<double> m_vUpwardTemperatureBoundaryConditionTime;

    //! Upstream boundary condition value vector
    std::vector<double> m_vUpwardTemperatureBoundaryConditionValue;

    //! Upstream temperature value (if constant)
    double m_dUpwardTemperatureBoundaryValue{};


    //! Upstream surface energy balance forcing coefficients
    double m_dUpwardTemperatureOffsetBeta{};
    double m_dUpwardAtmosphericExchangekA{};
    double m_dUpwardSurfaceConcentratedHeatkR{};
    double m_dUpwardInflowWaterEffectkQ{};

    //! Downstream boundary condition file (time series)
    std::string m_strDownwardTemperatureBoundaryConditionFilename;

    //! Downstream boundary condition time vector
    std::vector<double> m_vDownwardTemperatureBoundaryConditionTime;

    //! Downstream boundary condition value vector
    std::vector<double> m_vDownwardTemperatureBoundaryConditionValue;

    //! Downstream temperature value (if constant)
    double m_dDownwardTemperatureBoundaryValue{};

    //! Thermal diffusion coefficient (m²/s)
    double m_dThermalDispersion{};

    // === Surface energy balance forcing ===


    //! Single surface energy balance forcing file (Tair, relative humidity, wind)
    std::string m_strHeatFluxFile;

    //! Coeficiente de transferencia de calor sensible (CS) - Stanton number
    double m_dHeatFlux_CS{1.3e-3};

    //! Coeficiente de transferencia de calor latente (CL) - Dalton number
    double m_dHeatFlux_CL{1.3e-3};
    
    //! Site geographic latitude (decimal degrees, North positive)
    double m_dHeatFluxLatitude{36.5};
    
    //! Cloud cover for radiation calculation (0.0 = clear sky, 1.0 = overcast)
    double m_dHeatFluxCloudCover{0.2};
    
    //! Effective reservoir depth for 0D temperature model (m)
    double m_dReservoirEffectiveDepth{5.0};

    // === Forcing time series ===

    //! Heat flux forcing times and values (Tair, relative humidity, wind, pressure)
    std::vector<double> m_vHeatFluxTime;
    std::vector<double> m_vHeatFluxAirTemp;
    std::vector<double> m_vHeatFluxRelHumidity;
    std::vector<double> m_vHeatFluxWind;
    std::vector<double> m_vHeatFluxAtmosphericPressure;  // Atmospheric pressure (Pa)
    
    //! Daily minimum temperature vector calculated outside main loop
    //! Index = day from simulation start (0, 1, 2, ...)
    std::vector<double> m_vDailyMinTemperature;
    
    // === Simulation continuation options ===
    std::string m_strContinueNetcdfPath;
    
    // === Friend classes ===
    friend class CHydrograph;
    friend class CDataWriter;
    friend class CCrossSection;
    friend class CDataReader;
    friend class CScreenPresenter;


    public:
    //! System start-simulation time
    time_t m_tSysStartTime{};

    //! Sytem loop-start-simulation time
    time_t m_tSysStartLoopTime{};
    
    // Diagnostic tracking variables
    double m_dLastLogTime{0.0};      // Last time statistics were written
    int m_nWarningCount{0};          // Counter for anomalous value warnings

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
    
    //! YAML configuration content (for NetCDF metadata)
    string m_strYAMLConfigContent;

    //! Computational timestep obtained from Courant number
    double m_dTimestep;
    
    //! Adaptive timestep control variables
    double m_dTimestepPrevious;     // Previous timestep for comparison
    double m_dTimestepMin;          // Minimum timestep reached during simulation
    double m_dTimestepMax;          // Maximum allowed timestep (= m_dSimTimestep)
    int m_nTimestepChanges;         // Counter for significant timestep changes

    //! Lambda Value
    double m_dLambda;

    //! The initial estuarine condition, IEC [0 = in calm, 1 = water flow or 2 = elevation]
    int m_nInitialEstuarineCondition;

    //! The equation for Sediment Transport
    int m_nEquationSedimentTransport;

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

    //! TVD flux limiter working vectors for passive tracers (salinity, temperature, etc.)
    //! These are reusable for any passive scalar transport to avoid code duplication
    vector<double> m_vTracer_r_ratio;         // Smoothness indicator r = ∇S[i]/∇S[i-1]
    vector<double> m_vTracer_limiter;         // Limiter function Ψ(r) 
    vector<double> m_vTracer_flux_limited;    // TVD-limited flux at cell interfaces

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
    
    //! Salt mass balance tracking (downstream boundary)
    double m_dSaltInitialMass{0.0};        // Initial total salt mass in estuary at t=0 (kg)
    double m_dSaltInflowDownstream{0.0};   // Cumulative salt entering through ocean boundary (kg)
    double m_dSaltOutflowDownstream{0.0};  // Cumulative salt leaving through ocean boundary (kg)
    double m_dSaltNetFlowDownstream{0.0};  // Net cumulative flow (inflow - outflow) (kg)
    double m_dPredictorBoundaryFlux{0.0};  // Predictor boundary flux (kg/s) - for McCormack averaging
    double m_dCorrectorBoundaryFlux{0.0};  // Corrector boundary flux (kg/s) - for McCormack averaging
    
    //! The Courant Number (user-provided, used as upper limit only)
    double m_dCourantNumber{};
    
    //! Automatic CFL control based on numerical scheme stability
    double m_dCourantNumberMaxTheoretical; // Max CFL for scheme (TVD-McCormack)
    double m_dCourantNumberOptimal;        // Optimal CFL with safety factor
    double m_dCourantNumberCurrent;        // Current effective Courant (adaptive)
    double m_dCourantNumberMin;            // Minimum Courant reached
    int m_nCourantReductions;              // Counter for Courant reductions

    //! The equation for McCormack Limiter Flux (hydrodynamics)
    int m_nEquationMcCormackLimiterFlux{};

    //! The equation for transport limiter flux (salinity/temperature)
    //! If not specified, uses same as hydrodynamics
    int m_nEquationTransportLimiterFlux{};

    //! Psi Formula
    int m_nPsiFormula{};

    //! Delta Value
    double m_dDeltaValue{};

    //! Bathymetry smoothing parameters
    int m_nBathymetrySmoothingPasses{1};
    double m_dBathymetrySmoothingAlpha{0.25};

    //! Solution smoothing parameters
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

    //! Cross-section Bed slope for predictor (forward difference) - source term balance
    vector<double> m_vCrossSectionBedSlopePredictor;

    //! Cross-section Bed slope for corrector (backward difference) - source term balance
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
    //! Previous timestep cross-section areas (for salt mass conservation)
    vector<double> m_vPreviousCrossSectionArea;

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

    //! Cross-section DhDx (water surface gradient for pressure term)
    vector<double> m_vCrossSectionDhDx;

    //! Cross-section water densities (estado actual)
    vector<double> m_vCrossSectionDensity;
    //! Cross-section water densities (predicted, for corrector/baroclinic)
    vector<double> m_vPredictedCrossSectionDensity;
    //! Density gradient dRho/dx (baroclinic)
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

    //! Previous timestep velocity (for salt mass conservation)
    vector<double> m_vPreviousCrossSectionU;

    //! Cross-section perturbation water velocities
    vector<double> m_vCrossSectionC;

    //! Cross-section salinity vector (g/kg or PSU)
    vector<double> m_vCrossSectionSalinity;

    //! Cross-section temperature vector (°C)
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
    
    //! Storage for predictor and corrector fluxes (McCormack scheme)
    vector<double> m_vCrossSectionSalinityASt_predictor;
    vector<double> m_vCrossSectionSalinityASt_corrector;
    
    //! Salinity flux terms (analogous to F0, F1 for hydrodynamics)
    vector<double> m_vCrossSectionFS;      // Total salinity flux: advection + diffusion
    vector<double> m_vCrossSectionGS;      // Salinity source terms (currently unused, for future extensions)
    
    //! Previous timestep values for mass conservation
    vector<double> m_vPreviousCrossSectionSalinity;

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
    vector<double> m_vManningNumberSquared;      // Manning² (used in friction ~4000 times/day)
    vector<double> m_vInvDX;                     // 1/ΔX (used in gradients ~8000 times/day)
    vector<double> m_vDxSum;                     // ΔX[i] + ΔX[i+1] (centered differences)
    vector<double> m_vInvDxSum;                  // 1/(ΔX[i] + ΔX[i+1])
    vector<double> m_vGtimesDX;                  // g*ΔX (constant term)

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
    [[nodiscard]] double dGetSimulationDuration() const { return m_dSimDuration; }
    //! Method for setting the simulation duration
    void dSetSimulationDuration(double simDuration) { m_dSimDuration = simDuration; }

    //! Method for getting the simulation timestep
    [[nodiscard]] double dGetSimulationTimestep() const { return m_dSimTimestep; }
    //! Method for setting the simulation timestep
    void dSetSimulationTimestep(double simTimestep) { m_dSimTimestep = simTimestep; }

    //! Method for getting the initial estuarine condition
    [[nodiscard]] int nGetInitialEstuarineCondition() const { return m_nInitialEstuarineCondition; }
    //! Method for setting the initial estuarine condition
    void nSetInitialEstuarineCondition(int initialCondition) { m_nInitialEstuarineCondition = initialCondition; }

    //! Method for getting the compute sediment transport
    [[nodiscard]] bool bGetDoSedimentTransport() const { return m_bDoSedimentTransport; }
    //! Method for setting the compute sediment transport
    void bSetDoSedimentTransport(bool doSedimentTransport) { m_bDoSedimentTransport = doSedimentTransport; }

    //! Method for getting equation of sediment transport
    [[nodiscard]] int nGetEquationSedimentTransport() const { return m_nEquationSedimentTransport; }
    //! Method for setting equation of sediment transport
    void nSetEquationSedimentTransport(int equationSedimentTransport) { m_nEquationSedimentTransport = equationSedimentTransport; }

    //! Method for getting equation of sediment transport
    [[nodiscard]] double dGetSedimentDensity() const;
    //! Method for setting equation of sediment transport
    void dSetSedimentDensity(double sedimentDensity);

    //! Method for getting the compute water temperature
    [[nodiscard]] bool bGetDoWaterTemperature() const { return m_bDoWaterTemperature; }
    //! Method for setting the compute water temperature
    void bSetDoWaterTemperature(bool doWaterTemperature) { m_bDoWaterTemperature = doWaterTemperature; }

    //! Getter for Manning level-dependence
    [[nodiscard]] bool bGetManningDependsOnLevel() const { return m_bManningDependsOnLevel; }
    //! Setter for Manning level-dependence
    void bSetManningDependsOnLevel(bool bValue) { m_bManningDependsOnLevel = bValue; }

    //! Method for getting the compute water salinity
    [[nodiscard]] bool bGetDoWaterSalinity() const { return m_bDoWaterSalinity; }
    //! Method for setting the compute water salinity
    void bSetDoWaterSalinity(bool doWaterSalinity) { m_bDoWaterSalinity = doWaterSalinity; }

    //! Method for getting the compute water salinity
    [[nodiscard]] int nGetUpwardSalinityCondition() const { return m_nUpwardSalinityCondition; }
    //! Method for setting the compute water salinity
    void nSetUpwardSalinityCondition(int nUpwardCondition) { m_nUpwardSalinityCondition = nUpwardCondition; }

    //! Method for getting the compute water salinity
    [[nodiscard]] int nGetDownwardSalinityCondition() const { return m_nDownwardSalinityCondition; }
    //! Method for setting the compute water salinity
    void nSetDownwardSalinityCondition(int nDownwardCondition) { m_nDownwardSalinityCondition = nDownwardCondition; }

    //! Method for getting the beta salinity constant
    [[nodiscard]] double dGetBetaSalinityConstant() const { return m_dBetaSalinityConstant; }
    //! Method for setting the beta salinity constant
    void dSetBetaSalinityConstant(double salinityConstant) { m_dBetaSalinityConstant = salinityConstant; }
    //! Method for getting the beta temperature constant
    [[nodiscard]] double dGetBetaTemperatureConstant() const { return m_dBetaTemperatureConstant; }
    //! Method for setting the beta temperature constant
    void dSetBetaTemperatureConstant(double tempConstant) { m_dBetaTemperatureConstant = tempConstant; }

    //! Method for getting the along channel constant elevation
    [[nodiscard]] double dGetLongitudinalDispersionConstant() const { return m_dLongitudinalDispersion; }
    //! Method for setting the along channel constant elevation
    void dSetLongitudinalDispersionConstant(double longitudinalDispersion) { m_dLongitudinalDispersion = longitudinalDispersion; }

    //! Method for getting the upward estuarine condition
    [[nodiscard]] int nGetUpwardEstuarineCondition() const { return m_nUpwardEstuarineCondition; }
    //! Method for setting the upward estuarine condition
    void nSetUpwardEstuarineCondition(int upwardEstuarineCondition) { m_nUpwardEstuarineCondition = upwardEstuarineCondition; }

    //! Method for getting the downward estuarine condition
    [[nodiscard]] int nGetDownwardEstuarineCondition() const { return m_nDownwardEstuarineCondition; }
    //! Method for setting the downward estuarine condition
    void nSetDownwardEstuarineCondition(int downwardEstuarineCondition) { m_nDownwardEstuarineCondition = downwardEstuarineCondition; }

    //! Method for getting the courant number
    [[nodiscard]] double dGetCourantNumber() const { return m_dCourantNumber; }
    //! Method for setting the courant number
    void dSetCourantNumber(double courantNumber) { m_dCourantNumber = courantNumber; }

    //! Method for getting if McComarck limiter flux is applied
    [[nodiscard]] bool bGetDoMcComarckLimiterFlux() const { return m_bDoMcCormackLimiterFlux; }
    //! Method for setting if McComarck limiter flux is applied
    void bSetDoMcComarckLimiterFlux(bool doMcComarckLimiterFlux) { m_bDoMcCormackLimiterFlux = doMcComarckLimiterFlux; }

    //! Method for getting equation limiter flux (hydrodynamics)
    [[nodiscard]] int nGetEquationLimiterFlux() const { return m_nEquationMcCormackLimiterFlux; }
    //! Method for setting equation limiter flux (hydrodynamics)
    void nSetEquationLimiterFlux(int equationLimiterFlux) { m_nEquationMcCormackLimiterFlux = equationLimiterFlux; }

    //! Method for getting transport limiter flux (salinity/temperature)
    [[nodiscard]] int nGetTransportLimiterFlux() const { 
        return m_bUseIndependentTransportLimiter ? m_nEquationTransportLimiterFlux : m_nEquationMcCormackLimiterFlux; 
    }
    //! Method for setting transport limiter flux
    void nSetTransportLimiterFlux(int transportLimiterFlux, bool independent) { 
        m_nEquationTransportLimiterFlux = transportLimiterFlux; 
        m_bUseIndependentTransportLimiter = independent;
    }

    //! Method for getting Psi formula
    [[nodiscard]] int nGetPsiFormula() const { return m_nPsiFormula; }
    //! Method for setting Psi formula
    void nSetPsiFormula(int psiFormula) { m_nPsiFormula = psiFormula; }

    //! Method for getting Delta Value
    [[nodiscard]] double dGetDeltaValue() const { return m_dDeltaValue; }
    //! Method for setting Delta Value
    void dSetDeltaValue(double deltaValue) { m_dDeltaValue = deltaValue; }

    //! Method for getting if surface gradient method is applied
    [[nodiscard]] bool bGetDoSurfaceGradientMethod() const { return m_bDoSurfaceGradientMethod; }
    //! Method for setting if surface gradient method is applied
    void bSetDoSurfaceGradientMethod(bool doSurfaceGradientMethod) { m_bDoSurfaceGradientMethod = doSurfaceGradientMethod; }

    //! Method for getting if source term balance is applied
    [[nodiscard]] bool bGetDoSurfaceTermBalance() const { return m_bDoSourceTermBalance; }
    //! Method for setting if source term balance is applied
    void bSetDoSurfaceTermBalance(bool doSourceTermBalance) { m_bDoSourceTermBalance = doSourceTermBalance; }

    //! Method for getting if beta coefficient is applied
    [[nodiscard]] bool bGetDoBetaCoefficient() const { return m_bDoBetaCoefficient; }
    //! Method for setting if beta coefficient is applied
    void bSetDoBetaCoefficient(bool doBetaCoefficient) { m_bDoBetaCoefficient = doBetaCoefficient; }

    //! Method for getting if dry bed is applied
    [[nodiscard]] bool bGetDoDryBed() const { return m_bDoDryBed; }
    //! Method for setting if dry bed is applied
    void bSetDoDryBed(bool doDryBed) { m_bDoDryBed = doDryBed; }

    //! Method for getting if Murillo condition is applied
    [[nodiscard]] bool bGetDoMurilloCondition() const { return m_bDoMurilloCondition; }
    //! Method for setting if Murillo condition is applied
    void bSetDoMurilloCondition(bool doMurilloCondition) { m_bDoMurilloCondition = doMurilloCondition; }

    //! Method for getting the compute water density
    [[nodiscard]] bool bGetDoWaterDensity() const { return m_bDoWaterDensity; }
    //! Method for setting the compute water density
    void bSetDoWaterDensity(bool doWaterDensity) { m_bDoWaterDensity = doWaterDensity; }
    
    //! Method for getting smooth bathymetry flag
    [[nodiscard]] bool bGetDoSmoothBathymetry() const { return m_bDoSmoothBathymetry; }
    //! Method for setting smooth bathymetry flag
    void bSetDoSmoothBathymetry(bool doSmoothBathymetry) { m_bDoSmoothBathymetry = doSmoothBathymetry; }
    
    //! Method for getting smooth solution flag
    [[nodiscard]] bool bGetDoSmoothSolution() const { return m_bDoSmoothSolution; }
    //! Method for setting smooth solution flag
    void bSetDoSmoothSolution(bool doSmoothSolution) { m_bDoSmoothSolution = doSmoothSolution; }

    //! Method for getting save all timesteps flag
    [[nodiscard]] bool bGetSaveAllTimesteps() const { return m_bSaveAllTimesteps; }
    //! Method for setting save all timesteps flag
    void bSetSaveAllTimesteps(bool saveAllTimesteps) { m_bSaveAllTimesteps = saveAllTimesteps; }

    //! Add output variable
    void strAddOutputVariable(const string& strItem) { m_vOutputVariables.push_back(strItem); }

    //! A vector with hydrograph objects
    vector<CHydrograph> hydrographs;

    //! Method for getting the number of hydrographs
    [[nodiscard]] int nGetHydrographsNumber() const { return m_nHydrographsNumber; }
    //! Method for setting the number of hydrographs
    void nSetHydrographsNumber(int nValue) { m_nHydrographsNumber = nValue; }

    //! Get the vector of a variable
    vector<double> vGetVariable(const string& strVariableName) const;

    void AddHydrograph();



    //! Runs the simulation
    void bDoSimulation(int, char const* []);
    void initializeVectors();
    void calculateDailyMinTemperatures();
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
    void applyBoundariesToCurrentState();
    void applyBoundariesToPredictorState();
    void calculatePredictor();
    void calculateCorrector();
    void updatePredictorBoundaries();
    void updateCorrectorBoundaries();
    void mergePredictorCorrector();
    void mergeTracerPredictorCorrector();
    void smoothSolution();
    void smoothBathymetry();
    
    // Generic TVD flux limiter for passive tracer transport
    // Computes limited fluxes to prevent spurious oscillations while maintaining 2nd order accuracy
    void compute_tracer_tvd_flux(const vector<double>& tracer, const vector<double>& discharge,
                                  vector<double>& flux_limited);
    
    // Salinity transport (modular design like hydrodynamics)
    void calculateSalinityFluxes();        // Compute advective + diffusive fluxes
    void calculateSalinitySourceTerms();   // Compute source terms (if any)
    void calculate_salinity_predictor();   // Apply McCormack predictor step
    void calculate_salinity_corrector();   // Apply McCormack corrector step

    void calculateRadiativeFluxes();
    void calculate_temperature(); // Legacy wrapper (calls predictor/corrector based on m_nPredictor)
    void calculate_temperature_predictor();
    void calculate_temperature_corrector();

    void updateReservoirTemperature0D();

    //! Calculate adaptive Manning's n coefficient based on water level
    //! Returns n(η) varying linearly between 1.0 (mean water level) and 1.2 (one meter above mwl)
    static inline double n_eta(double eta) {
        if (eta <= 0) return 1.0;
        if (eta >= 1) return 1.2;
        return 1.0 + 0.2*eta;
    }
    
    void AnnounceProgress();
    void checkAnomalousValues();
    void writePeriodicStatistics();

    //! Carries out end-of-simulation tidying (error messages etc.)
    void bDoSimulationEnd();

    private:
    // Simulation start date properties (from CDataReader)
    int m_nSimStartSec;
    int m_nSimStartMin;
    int m_nSimStartHour;
    int m_nSimStartDay;
    int m_nSimStartMonth;
    int m_nSimStartYear;

    public:
    // Simulation start date getters
    int nGetSimStartSec() const { return m_nSimStartSec; }
    int nGetSimStartMin() const { return m_nSimStartMin; }
    int nGetSimStartHour() const { return m_nSimStartHour; }
    int nGetSimStartDay() const { return m_nSimStartDay; }
    int nGetSimStartMonth() const { return m_nSimStartMonth; }
    int nGetSimStartYear() const { return m_nSimStartYear; }

    // Simulation start date setters
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

    //! Get simulation end date and time as formatted string (ISO 8601)
    string getSimulationEndDateTimeString() const {
        // Calculate end time by adding simulation duration to start time
        std::tm start_tm = {};
        start_tm.tm_year = m_nSimStartYear - 1900;
        start_tm.tm_mon = m_nSimStartMonth - 1;
        start_tm.tm_mday = m_nSimStartDay;
        start_tm.tm_hour = m_nSimStartHour;
        start_tm.tm_min = m_nSimStartMin;
        start_tm.tm_sec = m_nSimStartSec;
        
        std::time_t start_time = std::mktime(&start_tm);
        std::time_t end_time = start_time + static_cast<std::time_t>(m_dSimDuration);
        std::tm* end_tm = std::localtime(&end_time);
        
        ostringstream oss;
        oss << setfill('0') << setw(4) << (end_tm->tm_year + 1900) << "-"
            << setw(2) << (end_tm->tm_mon + 1) << "-"
            << setw(2) << end_tm->tm_mday << " "
            << setw(2) << end_tm->tm_hour << ":"
            << setw(2) << end_tm->tm_min << ":"
            << setw(2) << end_tm->tm_sec;
        return oss.str();
    }

    private:
    // === Private data members for internal precomputation ===
    vector<double> m_vBedZ;
    vector<double> m_vManningN;
    vector<double> m_vPositionX;
    vector<double> m_vBeta;
    
    vector<vector<double>> m_vWidth;
    vector<vector<double>> m_vLeftY;
    vector<vector<double>> m_vRightY;
    vector<vector<double>> m_vEstuaryAreas;
    vector<vector<double>> m_vEstuaryHydraulicRadius;
    vector<vector<double>> m_vEstuaryWaterDepths;
    vector<vector<double>> m_vPrecalculatedSecondTerm;
    vector<int> m_vElevationSectionsCount;

    public:
    // === Optimization methods ===
    void precomputeEstuaryData();
};
#endif // SIMULATION_H