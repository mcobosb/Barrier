/*!
*
 * \file cross_section.cpp
 * \brief Cross-section geometry and hydraulic properties management
 * \details Handles cross-section geometric data (width, depth, area)
 *          and computes derived hydraulic parameters (hydraulic radius,
 *          wetted perimeter, conveyance) for open channel flow calculations.
 * \author Manuel Cobos Budia

 * \date 2026
 * \copyright GNU General Public License
 *
 */

/*==============================================================================================================================

==============================================================================================================================*/
#include <fstream>
using std::ifstream;

#include <iostream>
using std::cerr;
using std::cout;
using std::endl;
using std::ios;

#include <string>
using std::to_string;

#include <algorithm>
using std::find;

#include <cmath>

#include <cross_section.h>

/**
 * @brief Construct a new CCrossSection object for channel geometry
 * 
 * Initializes all scalar members to 0:
 * - Section numbers (ID, elevation count)
 * - Position (X, Z, UTM coordinates)
 * - Hydraulic properties (Manning's n, β coefficient)
 * - Bank angles (left/right)
 * 
 * @note Vector members (width, area, hydraulic radius tables) are empty
 * @see dAppend2Vector() for populating tabulated hydraulic data
 */
CCrossSection::CCrossSection(){
    m_nSectionNumber =
    m_nElevationSectionNumber = 0;

    m_dX =
    m_dElevation =
    m_dZ =
    m_dManningNumber =
    m_dX_UTM =
    m_dY_UTM =
    m_dRightRBAngle =
    m_dLeftRBAngle =
    m_dBeta =  0.0;
}

/**
 * @brief Destructor (default implementation)
 */
CCrossSection::~CCrossSection() = default;


/**
 * @brief Append hydraulic data to tabulated vectors
 * 
 * Recognized items (case-sensitive):
 * - "elevation": Water depth (m) relative to bed
 * - "width": Channel width (m) at this elevation
 * - "area": Cross-sectional area (m²)
 * - "perimeter": Wetted perimeter (m)
 * - "hydraulic radius": Rh = A/P (m)
 * - "sigma": Width function σ(η)
 * - "left/right river bank location": Bank positions
 * - "I1", "I2": Pressure integral coefficients
 * 
 * @param strItem Column name from CSV file
 * @param dValue Numerical value to append
 * 
 * @note Called during CSV parsing (data_reader.cpp)
 * @see calculateI1() for computing I1 from geometry
 */
void CCrossSection::dAppend2Vector(const string& strItem, double dValue){
    if (strItem == "elevation")
        m_vWaterDepth.push_back(dValue);
    else if (strItem == "width")
        m_vWidth.push_back(dValue);
    else if (strItem == "area")
        m_vArea.push_back(dValue);
    else if (strItem == "perimeter")
        m_vPerimeter.push_back(dValue);
    else if (strItem == "hydraulic radius")
        m_vHydraulicRadius.push_back(dValue);
    else if (strItem == "sigma")
        m_vSigma.push_back(dValue);
    else if (strItem == "left river bank location")
        m_vLeftRBLocation.push_back(dValue);
    else if (strItem == "right river bank location")
        m_vRightRBLocation.push_back(dValue);
    else if (strItem == "I1")
        m_vI1.push_back(dValue);
    else if (strItem == "I2")
        m_vI2.push_back(dValue);
}

/**
 * @brief Setters for scalar cross-section properties
 * 
 * Simple assignment functions for:
 * - Section ID, elevation count
 * - Position (X, Z, UTM X/Y)
 * - Manning's n roughness coefficient
 * - Bank angles (degrees from north)
 * - β: Momentum correction coefficient
 * - Water depth (current simulation value)
 */
void CCrossSection::nSetSectionNumber(const int nValue){
    m_nSectionNumber = nValue;
}

void CCrossSection::dSetX(const double dValue){
    m_dX = dValue;
}
void CCrossSection::dSetZ(const double dValue){
    m_dZ = dValue;
}

void CCrossSection::dSetManningNumber(const double dValue) {
    m_dManningNumber = dValue;
}
void CCrossSection::dSetX_UTM(const double dValue) {
    m_dX_UTM = dValue;
}
void CCrossSection::dSetY_UTM(const double dValue) {
    m_dY_UTM = dValue;
}
void CCrossSection::dSetRightRBAngle(const double dValue) {
    m_dRightRBAngle = dValue;
}
void CCrossSection::dSetLeftRBAngle(const double dValue) {
    m_dLeftRBAngle = dValue;
}
void CCrossSection::dSetWaterDepth(const double dValue) {
    m_dElevation = dValue;
}
void CCrossSection::dSetBeta(const double dValue) {
    m_dBeta = dValue;
}
void CCrossSection::nSetElevationSectionsNumber(const int nValue) {
    m_nElevationSectionNumber = nValue;
}

/**
 * @brief Getters for scalar and indexed vector properties
 * 
 * Two types:
 * 1. Scalar getters: Return single values (X, Z, Manning's n, etc.)
 * 2. Indexed getters: Return value at specific elevation index
 *    e.g., dGetArea(5) returns area at 5th elevation level
 * 
 * @note Indexed getters have no bounds checking - caller must validate
 */
int CCrossSection::nGetSectionNumber() const{
  return m_nSectionNumber;
}
double CCrossSection::dGetX() const{
  return m_dX;
}
double CCrossSection::dGetZ() const {
  return m_dZ;
}
double CCrossSection::dGetManningNumber() const {
  return m_dManningNumber;
}
double CCrossSection::dGetX_UTM() const{
  return m_dX_UTM;
}
double CCrossSection::dGetY_UTM() const{
  return m_dY_UTM;
}
double CCrossSection::dGetRightRBAngle() const{
  return m_dRightRBAngle;
}
double CCrossSection::dGetLeftRBAngle() const{
  return m_dLeftRBAngle;
}
double CCrossSection::dGetBeta() const{
  return m_dBeta;
}
int CCrossSection::nGetElevationSectionsNumber() const{
    return m_nElevationSectionNumber;
}
double CCrossSection::dGetWaterDepth(const int nValue) const {
  return m_vWaterDepth[nValue];
}
double CCrossSection::dGetWidth(const int nValue) const {
  return m_vWidth[nValue];
}
double CCrossSection::dGetArea(const int nValue) const {
  return m_vArea[nValue];
}
double CCrossSection::dGetPerimeter(const int nValue) const {
  return m_vPerimeter[nValue];
}
double CCrossSection::dGetHydraulicRadius(const int nValue) const {
  return m_vHydraulicRadius[nValue];
}
double CCrossSection::dGetSigma(const int nValue) const {
  return m_vSigma[nValue];
}
double CCrossSection::dGetLeftY(const int nValue) const {
  return m_vLeftRBLocation[nValue];
}
double CCrossSection::dGetRightY(const int nValue) const {
  return m_vRightRBLocation[nValue];
}
double CCrossSection::dGetI1(const int nValue) const {
    return m_vI1[nValue];
}
double CCrossSection::dGetI2(const int nValue) const {
    return m_vI2[nValue];
}

/**
 * @brief Getters for complete hydraulic tables (returns copies)
 * 
 * Returns entire tabulated vectors:
 * - Area(elevation), Width(elevation), Hydraulic Radius(elevation)
 * - Water depth levels, Bank locations
 * - I1, I2 pressure integral coefficients
 * 
 * @note Returns by value (copy) - use carefully in tight loops
 * @see precomputeEstuaryData() in simulation.cpp for cached access
 */
vector<double> CCrossSection::vGetArea() {
    return m_vArea;
}
vector<double> CCrossSection::vGetHydraulicRadius() {
    return m_vHydraulicRadius;
}
vector<double> CCrossSection::vGetWaterDepth() {
    return m_vWaterDepth;
}
vector<double> CCrossSection::vGetWidth() {
    return m_vWidth;
}
vector<double> CCrossSection::vGetLeftRBLocation() {
    return m_vLeftRBLocation;
}
vector<double> CCrossSection::vGetRightRBLocation() {
    return m_vRightRBLocation;
}
vector<double> CCrossSection::vGetI1() {
    return m_vI1;
}
vector<double> CCrossSection::vGetI2() {
    return m_vI2;
}

/**
 * @brief Calculate I1 pressure integral for all elevation levels
 * 
 * Computes the first moment of pressure distribution:
 *   I1(h) = ∫₀ʰ (h-η)·σ(η)·dη
 * 
 * where:
 * - h: Total water depth (m)
 * - η: Elevation from bed (m, 0 ≤ η ≤ h)
 * - σ(η): Channel width at elevation η (m)
 * - (h-η): Pressure head at elevation η (m)
 * 
 * Physical interpretation:
 * - I1 represents the hydrostatic pressure moment about the water surface
 * - Used in momentum equation for non-prismatic channels
 * - For rectangular channel: I1 = 0.5·A·h
 * - For trapezoidal channel: I1 ≈ 0.4·A·h (wider at top)
 * 
 * Numerical method:
 * - Trapezoidal rule integration over tabulated σ(η)
 * - Loops over all elevation levels i = 0 to N-1
 * - For each level, integrates from bed (j=0) to current depth (j=i)
 * 
 * @note Called once during initialization (data_reader.cpp)
 * @see Chaudhry (2008): Open-Channel Flow, Section 4.3
 */
void CCrossSection::calculateI1() {
    // Clear any existing I1 values
    m_vI1.clear();
    m_vI1.resize(m_vWaterDepth.size(), 0.0);
    
    // For each elevation level (water depth h)
    for (size_t i = 0; i < m_vWaterDepth.size(); i++) {
        if (i == 0 || m_vArea[i] < 1e-6) {
            m_vI1[i] = 0.0;
            continue;
        }
        
        double h = m_vWaterDepth[i];  // Total water depth at this level
        double I1 = 0.0;
        
        // Integrate I1 = ∫₀ʰ (h-η)·σ(η)·dη using trapezoidal rule
        // Integrate from bottom (j=0) to current water level (j=i-1)
        for (size_t j = 0; j < i; j++) {
            double eta_j = m_vWaterDepth[j];
            double eta_j1 = (j+1 < m_vWaterDepth.size()) ? m_vWaterDepth[j+1] : h;
            double sigma_j = m_vWidth[j];
            double sigma_j1 = (j+1 < m_vWidth.size()) ? m_vWidth[j+1] : m_vWidth[j];
            
            // For the last segment, use current level
            if (j + 1 >= i) {
                eta_j1 = h;
                sigma_j1 = m_vWidth[i];
            }
            
            // Integrand: (h-η)·σ(η)
            double f_j = (h - eta_j) * sigma_j;
            double f_j1 = (h - eta_j1) * sigma_j1;
            
            // Trapezoidal rule: Δη/2 · (f_j + f_j+1)
            double deta = eta_j1 - eta_j;
            if (deta > 1e-10) {
                I1 += 0.5 * deta * (f_j + f_j1);
            }
        }
        
        m_vI1[i] = I1;
    }
}


/**
 * @brief Find elevation where channel width gradient exceeds threshold
 * 
 * Searches for first elevation η > maxAstronomicalTide where dB/dη > threshold.
 * 
 * Purpose:
 * - Detect transition from main channel to floodplain
 * - Identify vegetation emergence zones (salt marshes)
 * - Adjust Manning's n dynamically based on inundation
 * 
 * Algorithm:
 * 1. Compute width gradient dB/dη at each elevation
 * 2. Find first η > maxTide where gradient exceeds threshold
 * 3. Return midpoint elevation of that segment
 * 
 * @param maxAstronomicalTide Maximum tide level from tidal analysis (m)
 * @param threshold Minimum gradient to detect (m/m, e.g., 0.5 = 50% slope)
 * @param etaWidthGradientThreshold Output: elevation where gradient exceeds threshold (m)
 * 
 * @note Returns 0.0 if no gradient exceeds threshold or insufficient data
 */
void CCrossSection::calculateEtaMaxWidthGradient(double maxAstronomicalTide, double threshold, double& etaWidthGradientThreshold) const {
    double etaMaxWidthGradient = 0.0;
    if (m_vWaterDepth.size() < 2 || m_vWidth.size() < 2) {
        etaWidthGradientThreshold = 0.0;
        return;
    }
    for (size_t i = 1; i < m_vWaterDepth.size(); ++i) {
        double deta = m_vWaterDepth[i] - m_vWaterDepth[i-1];
        if (std::fabs(deta) < 1e-8) continue;
        double dB = m_vWidth[i] - m_vWidth[i-1];
        double grad = dB / deta;
        double etaMid = 0.5 * (m_vWaterDepth[i] + m_vWaterDepth[i-1]);
        if (etaMid > maxAstronomicalTide && grad > threshold) {
            etaMaxWidthGradient = etaMid;
            break;
        }
    }
    etaWidthGradientThreshold = etaMaxWidthGradient;
}