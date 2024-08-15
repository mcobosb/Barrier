
//===============================================================================================================================
//! The nDoSimulation member function of CSimulation sets up and runs the simulation
//===============================================================================================================================
int CSimulation::nDoSimulation(int nArg, char const* pcArgv[])
{
#ifdef RANDCHECK
    CheckRand();
    return RTN_OK;
#endif

    // ================================================== initialization section ================================================
    // Hello, World!
    AnnounceStart();

    // Start the clock ticking
    StartClock();

    // Find out the folder in which the CoastalME executable sits, in order to open the .ini file (they are assumed to be in the same folder)
    if (! bFindExeDir(pcArgv[0]))
        return (RTN_ERR_CMEDIR);

    // Deal with command-line parameters
    int nRet = nHandleCommandLineParams (nArg, pcArgv);
    if (nRet != RTN_OK)
        return (nRet);

    // OK, we are off, tell the user about the licence and the start time
    AnnounceLicence();

    // Read the .ini file and get the name of the run-data file, and path for output etc.
    if (! bReadIniFile())
        return (RTN_ERR_INI);