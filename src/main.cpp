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

//===============================================================================================================================
//! Saint Venant main function
//===============================================================================================================================
int main(int argc, char const* argv[])
{
    // Enable the use of UTF-8 symbols in Saint Venant output
    // setlocale(LC_ALL, "en_GB.UTF-8");

    //! Initialize objects ---------------------------------------------------------------------------------------------
    // Create a CSimulation object
    CSimulation pSimulation;

    // Run the simulation and then check how it ends
    pSimulation.bDoSimulation(argc, argv);
    pSimulation.bDoSimulationEnd();

   // Go back to the OS
   return 0;
}
