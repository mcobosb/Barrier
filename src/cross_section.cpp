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
    m_dBeta = 0.0;
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
 * 
 * @param strItem Column name from CSV file
 * @param dValue Numerical value to append
 * 
 * @note Called during CSV parsing (data_reader.cpp)
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

/**
 * @brief Getters for complete hydraulic tables (no copies)
 * 
 * Returns references to entire tabulated vectors:
 * - Area(elevation), Width(elevation), Hydraulic Radius(elevation)
 * - Water depth levels, Bank locations
 * 
 * @note Returns by const reference (no copy). This is safe as long as
 *       the CCrossSection object outlives the reference.
 */
const vector<double>& CCrossSection::vGetArea() const {
  return m_vArea;
}
const vector<double>& CCrossSection::vGetHydraulicRadius() const {
  return m_vHydraulicRadius;
}
const vector<double>& CCrossSection::vGetWaterDepth() const {
  return m_vWaterDepth;
}
const vector<double>& CCrossSection::vGetWidth() const {
  return m_vWidth;
}
const vector<double>& CCrossSection::vGetLeftRBLocation() const {
  return m_vLeftRBLocation;
}
const vector<double>& CCrossSection::vGetRightRBLocation() const {
  return m_vRightRBLocation;
}