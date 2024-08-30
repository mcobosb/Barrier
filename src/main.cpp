/*!
 *
 * \file main.cpp
 * \brief The start-up routine for Saint Venant solver
 * \details TODO 001 A more detailed description of this routine
 * \author Manuel Cobos Budia

 * \date 2024
 * \copyright GNU General Public License
 *
 */

/*===============================================================================================================================
===============================================================================================================================*/
#include "main.h"
#include "simulation.h"
#include "screen_presenter.h"
#include "data_reader.h"

//===============================================================================================================================
//! Saint Venant's main function
//===============================================================================================================================
int main(int argc, char const* argv[])
{
    // Enable the use of UTF-8 symbols in Saint Venant output
    // setlocale(LC_ALL, "en_GB.UTF-8");

    //! Initialize objects ---------------------------------------------------------------------------------------------
    // Create a CSimulation object
    CSimulation pSimulation;

    // Create a ScreenPresenter object
    auto* pScreenPresenter = new CScreenPresenter;
    // Create a DataReader object
    auto* pDataReader = new CDataReader();



    //! Announce Start -------------------------------------------------------------------------------------------------
    pScreenPresenter->StartingRun(argc, argv);

    bool nRtn = pDataReader->bOpenLogFile(&pSimulation);
    //! Read the .ini file and get the name of the run-data file, and path for output etc.
    nRtn = pDataReader->bReadConfigurationFile(&pSimulation);
    nRtn = pDataReader->bReadAlongChannelGeometryFile(&pSimulation);
    nRtn = pDataReader->bReadCrossSectionGeometryFile(&pSimulation);
    nRtn = pDataReader->bReadHydrographsFile(&pSimulation);


    // Run the simulation and then check how it ends
    nRtn = pSimulation.bDoSimulation(argc, argv);
    // pSimulation->DoSimulationEnd(nRtn);

    //! Announce End -------------------------------------------------------------------------------------------------
    pScreenPresenter->EndingRun();

   // Get rid of the CSimulation object and close files
   // delete pSimulation;

   // Go back to the OS
   return nRtn;
}
