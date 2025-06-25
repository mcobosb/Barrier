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
#include <iostream>
#include <exception>

//===============================================================================================================================
//! Saint Venant main function
//===============================================================================================================================
int main(int argc, char const* argv[])
{
    // Enable the use of UTF-8 symbols in Saint Venant output
    // setlocale(LC_ALL, "en_GB.UTF-8");

    int returnCode = 0;
    
    try {
        //! Initialize objects ---------------------------------------------------------------------------------------------
        // Create a CSimulation object
        CSimulation pSimulation;

        // Run the simulation and then check how it ends
        pSimulation.bDoSimulation(argc, argv);
        pSimulation.bDoSimulationEnd();
        
        // std::cout << "Simulation completed successfully." << std::endl;
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
