/*!
*
 * \class CCrossSection
 * \brief
 * \details TODO 001 This is a more detailed description of the CCrossSection class
 * \author Manuel Cobos Budia

 * \date 2024
 * \copyright GNU General Public License
 *
 * \file cross_section.h
 * \brief Contains Cross-Section definitions
 *
 */
#ifndef CROSS_SECTION_H
#define CROSS_SECTION_H

#include <string>
using std::string;

#include <vector>
using std::vector;

class CCrossSection {
    private:
        //! Section number
        int m_nSectionNumber;

        //! The X coordinate
        double m_dX;

        //! The Z coordinate
        double m_dZ;

        //! The Manning number
        double m_dManningNo;

        //! The X coordinate in UTM
        double m_dXUTM;

        //! The Y coordinate in UTM
        double m_dYUTM;

        //! The right river bank angle
        double m_dRightRBAngle;

        //! The left river bank angle
        double m_dLeftRBAngle;

        //! Number of Elevation Sections
        int m_nElevationSectionNumber;

        //! The elevation over Z
        vector<double> m_vElevation;

        //! The channel width at Elevation
        vector<double> m_vWidth;

        //! The channel area at Elevation
        vector<double> m_vArea;

        //! The perimeter at Elevation
        vector<double> m_vPerimeter;

        //! The Hydraulic Radius at Elevation
        vector<double> m_vHydraulicRadius;

        //! The sigma value at Elevation
        vector<double> m_vSigma;

        //! The location of the left riverbank at Elevation
        vector<double> m_vLeftRBLocation;

        //! The location of the right riverbank at Elevation
        vector<double> m_vRightRBLocation;

        //! The beta value at Elevation
        vector<double> m_vBeta;

        //! The I1 value at Elevation
        vector<double> m_vI1;

        //! The I2 value at Elevation
        vector<double> m_vI2;

    public:

        CCrossSection();
        ~CCrossSection();

        void nSetSectionNumber(int nValue);
        void dSetX(double dValue);
        void dSetZ(double dValue);
        void dSetManningNo(double dValue);
        void dSetXUTM(double dValue);
        void dSetYUTM(double dValue);
        void dSetRightRBAngle(double dValue);
        void dSetLeftRBAngle(double dValue);
        void nSetElevationSectionsNumber(int nValue);
        void dAppend2Vector(string strItem, double dValue);


        // void dSetElevation(double dValue);
        // void dSetWidth(double dValue);
//        void dSetArea(double dArea);
//        void dSetPerimeter(double dPerimeter);
//        void dSetHydraulicRadius(double dHydraulicRadius);
//        void dSetSigma(double dSigma);
//        void dSetLeftY(double dLeftY);
//        void dSetRightY(double dRightY);
//        void dSetBeta(double dBeta);

        int nGetSectionNumber();
        double dGetX();
        double dGetZ();
        double dGetManningNumber();
        double dGetXUTM();
        double dGetYUTM();
        double dGetRightRBAngle();
        double dGetLeftRBAngle();
        int nGetElevationSectionsNumber();
        double dGetElevation(int nValue);
        double dGetWidth(int nValue);
        double dGetArea(int nValue);
        double dGetPerimeter(int nValue);
        double dGetHydraulicRadius(int nValue);
        double dGetSigma(int nValue);
        double dGetLeftY(int nValue);
        double dGetRightY(int nValue);
        double dGetBeta(int nValue);
        double dGetI1(int nValue);
        double dGetI2(int nValue);

        vector<double> vGetArea();
        vector<double> vGetHydraulicRadius();
        vector<double> vGetElevation();
};
#endif