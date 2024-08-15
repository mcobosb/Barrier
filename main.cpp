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

//===============================================================================================================================
//! Saint Venant's main function
//===============================================================================================================================
int main(int argc, char const* argv[])
{
   // Enable the use of UTF-8 symbols in Saint Venant output
   setlocale(LC_ALL, "en_GB.UTF-8");

   // Create a CSimulation object
   auto* pSimulation = new CSimulation;

   // Run the simulation and then check how it ends
   int nRtn = pSimulation->nDoSimulation(argc, argv);
   pSimulation->DoSimulationEnd(nRtn);

   // Get rid of the CSimulation object and close files
   delete pSimulation;

   // Go back to the OS
   return nRtn;
}
