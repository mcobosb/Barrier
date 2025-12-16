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

#include <cross_section.h>

//===============================================================================================================================
//! The CCrossSection constructor
//===============================================================================================================================
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

//===============================================================================================================================
//! The CCrossSection destructor
//===============================================================================================================================
CCrossSection::~CCrossSection() = default;


//===============================================================================================================================
//! Append a cross-section object to the estuary object
//===============================================================================================================================
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

//===============================================================================================================================
//! The CCrossSection Setters
//===============================================================================================================================
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

//===============================================================================================================================
//! The CCrossSection Getters
//===============================================================================================================================
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

//! Getter for vector variables
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

//===============================================================================================================================
//! Calculate I1 pressure integral for each elevation
//! I1 = ∫₀ʰ (h-η)·σ(x,η)·dη where σ is width, η is elevation from bed, h is water depth
//! This integral represents the first moment of the pressure distribution about the water surface
//! For trapezoidal channels: I1 ≈ 0.4·A·h (compared to rectangular: I1 = 0.5·A·h)
//===============================================================================================================================
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
