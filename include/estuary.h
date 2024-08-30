/*!
*
* \class CEstuary
* \brief
* \details TODO 001 This is a more detailed description of the CSimulation class
* \author Manuel Cobos Budia

* \date 2024
* \copyright GNU General Public License
*
* \file estuary.h
* \brief Contains CEstuary definitions
*
*/

#ifndef ESTUARY_H
#define ESTUARY_H

#include <vector>
using std::vector;

#include "cross_section.h"

class CEstuary{
    //! The CCrossSection class is a friend of the CEstuary
    friend class CCrossSection;

    private:

        // void nSetSectionNumber(int nSectionNumber);
        // void dSetX(double dX);
        // void dSetZ(double dZ);
        // void dAppendVector(vector<double>, double dValue);
        //
        //
        // int nGetSectionNumber();
        // double dGetX();
        // double dGetZ();
        // double dGetElevation(int nZ);
        // double dGetWidth(int nZ);
        // double dGetArea(int nZ);
        // double dGetPerimeter(int nZ);
        // double dGetHydraulicRadius(int nZ);
        // double dGetSigma(int nZ);
        // double dGetLeftY(int nZ);
        // double dGetRightY(int nZ);
        // double dGetBeta(int nZ);



    public:
        CEstuary();
        ~CEstuary();

        //! A vector with cross-sections objects along the estuary
        vector<CCrossSection> estuary;

        void AddCrossSection();
};
#endif