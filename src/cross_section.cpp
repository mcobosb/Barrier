/*!
*
 * \file cross_section.cpp
 * \brief
 * \details TODO 001 A more detailed description of these routines.
 * \author Manuel Cobos Budia

 * \date 28/08/2024
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
    m_dLeftRBAngle =
    m_dX_UTM =
    m_dY_UTM =
    m_dRightRBAngle =
    m_dLeftRBAngle =  0.0;
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
    else if (strItem == "beta")
       m_vBeta.push_back(dValue);
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
double CCrossSection::dGetBeta(const int nValue) const {
  return m_vBeta[nValue];
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
