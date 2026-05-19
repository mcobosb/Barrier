/*!
 *
 * \file main.cpp
 * \brief Entry point for the Barrier estuarine hydrodynamic model
 * \details Initializes the simulation engine, handles execution flow,
 *          catches exceptions, and returns appropriate exit codes.
 * \author Manuel Cobos Budia

 * \date 2026
 * \copyright GNU General Public License
 *
 */

/*===============================================================================================================================
===============================================================================================================================*/
#include "main.h"
#include "simulation.h"
#include <iostream>
#include <exception>

/**
 * @brief Main entry point for Barrier 1D estuarine hydrodynamic model
 * 
 * Execution flow:
 * 1. Create CSimulation object (reads config, allocates memory)
 * 2. Run simulation loop: bDoSimulation(argc, argv)
 * 3. Finalize and report results: bDoSimulationEnd()
 * 4. Handle exceptions and return appropriate exit code
 * 
 * @param argc Number of command-line arguments
 * @param argv Array of argument strings (argv[1] = config file path)
 * @return 0 on success, 1 for std::exception, 2 for unknown errors
 * 
 * @note Exception safety:
 * - std::exception catches: file I/O, memory allocation, math errors
 * - Unknown exceptions (...)  catch potential C-style errors
 * 
 * @see CSimulation::bDoSimulation() for main loop
 * @see CSimulation::bDoSimulationEnd() for cleanup/output
 */
int main(int argc, char const* argv[])
{
    int returnCode = 0;
    
    try {
        // Initialize simulation engine from config file
        CSimulation pSimulation;

        // Execute simulation loop (predictor-corrector steps)
        pSimulation.bDoSimulation(argc, argv);
        
        // Finalize: close NetCDF, print summary, cleanup
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
