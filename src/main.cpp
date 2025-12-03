/*!
 *
 * \file main.cpp
 * \brief Entry point for the Barrier estuarine hydrodynamic model
 * \details Initializes the simulation engine, handles execution flow,
 *          catches exceptions, and returns appropriate exit codes.
 * \author Manuel Cobos Budia

 * \date 2024
 * \copyright GNU General Public License
 *
 */

/*===============================================================================================================================
===============================================================================================================================*/
#include "main.h"
#include "simulation.h"
#include <iostream>
#include <exception>

//===============================================================================================================================
//! Saint Venant main function
//===============================================================================================================================
int main(int argc, char const* argv[])
{
    int returnCode = 0;
    
    try {
        //! Initialize objects ---------------------------------------------------------------------------------------------
        // Create a CSimulation object
        CSimulation pSimulation;

        // Run the simulation and then check how it ends
        pSimulation.bDoSimulation(argc, argv);
        pSimulation.bDoSimulationEnd();
    }
    catch (const std::exception& e) {
        std::cerr << "Error during simulation: " << e.what() << std::endl;
        returnCode = 1;
    }
    catch (...) {
        std::cerr << "Unknown error during simulation" << std::endl;
        returnCode = 2;
    }
    
    // Go back to the OS
    return returnCode;
}
