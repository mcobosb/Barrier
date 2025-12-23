🧩 Calibration 
- [1] Prepare water level calibration routines
- [ ] Implement current (water flow) calibration
- [ ] Implement salinity calibration
- [ ] Implement temperature calibration
- [ ] Calibrate parameters: Manning coefficient n(x), longitudinal dispersion (Kx), CL and CF

🧩 Model Improvements
- [X] Add option to specify constant PSU value at the estuary mouth, or provide a time series file for mouth salinity
- [X] Test initial salinity time series input at the mouth
- [X] Add temperature transport
- [X] Test temperature transport
- [X] Add baroclinic transport terms
- [X] Add the dependence of n(eta) with the tidal range
- [X] Add the threshold for jumping to inundation flats
- [ ] Implement suspended sediment transport (affects density) and bedload transport
- [X] Include density calculation in the model
- [X] Smooth considering distante between cross-sections (non-uniform 1D Laplacian)
    - [X] Read alpha and num_passes from config.yaml
- [X] Change reading beta coefficient (F1 = β*(Q²/A + g*I1)) from along_channel file
- [ ] Add LAPE (Longitudinal Anomaly of Potential Energy) calculation
- [X] Test salinity initial condition input file
- [ ] Write log files with key information
- [X] Add upstream 0D model for well-mixed tanks inflows temperature 

🧩 Data Input
- [X] Adapt code to allow specifying a NetCDF file path for reading parameters and continue the simulation
- [X] Adapt yaml_reader for considering hydrodynamics in a block as transports

🧩 Simulation Continuation
- [X] Enable continuing the simulation

