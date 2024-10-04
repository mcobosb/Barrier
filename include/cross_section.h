/*!
*
 * \class CCrossSection
 * \brief
 * \details Definition of CCrossSection class that represents estuarine cross-section geometry.
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
        double m_dManningNumber;

        //! The X coordinate in UTM
        double m_dX_UTM;

        //! The Y coordinate in UTM
        double m_dY_UTM;

        //! The right river bank angle
        double m_dRightRBAngle;

        //! The left river bank angle
        double m_dLeftRBAngle;

        //! The water elevation
        double m_dElevation;

        //! Number of Elevation Sections
        int m_nElevationSectionNumber;

        //! The water depth
        vector<double> m_vWaterDepth;

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

        //! Setter
        void nSetSectionNumber(int nValue);
        void dSetX(double dValue);
        void dSetZ(double dValue);
        void dSetManningNumber(double dValue);
        void dSetX_UTM(double dValue);
        void dSetY_UTM(double dValue);
        void dSetRightRBAngle(double dValue);
        void dSetLeftRBAngle(double dValue);
        void nSetElevationSectionsNumber(int nValue);
        void dSetWaterDepth(double dValue);
        void dAppend2Vector(const string& strItem, double dValue);

        //! Getter
        [[nodiscard]] int nGetSectionNumber() const;
        [[nodiscard]] double dGetX() const;
        [[nodiscard]] double dGetZ() const;
        [[nodiscard]] double dGetManningNumber() const;
        [[nodiscard]] double dGetX_UTM() const;
        [[nodiscard]] double dGetY_UTM() const;
        [[nodiscard]] double dGetRightRBAngle() const;
        [[nodiscard]] double dGetLeftRBAngle() const;
        [[nodiscard]] int nGetElevationSectionsNumber() const;
        [[nodiscard]] double dGetWaterDepth(int nValue) const;
        [[nodiscard]] double dGetWidth(int nValue) const;
        [[nodiscard]] double dGetArea(int nValue) const;
        [[nodiscard]] double dGetPerimeter(int nValue) const;
        [[nodiscard]] double dGetHydraulicRadius(int nValue) const;
        [[nodiscard]] double dGetSigma(int nValue) const;
        [[nodiscard]] double dGetLeftY(int nValue) const;
        [[nodiscard]] double dGetRightY(int nValue) const;
        [[nodiscard]] double dGetBeta(int nValue) const;
        [[nodiscard]] double dGetI1(int nValue) const;
        [[nodiscard]] double dGetI2(int nValue) const;

        //! Vector getter
        vector<double> vGetArea();
        vector<double> vGetHydraulicRadius();
        vector<double> vGetWaterDepth();
};
#endif