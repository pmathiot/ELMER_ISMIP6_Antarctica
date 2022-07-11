# ELMER_ISMIP6_Antarctica

Repository to work on the Elmer/Ice IGE contribution to ISMIP6.

## TODO

### Input
- [ ] Setup a data directory and add it in the README
- [ ] SMB input file
    - [ ] Control
    - [ ] Scenario XXX
    - [ ] ... to be completed
- [ ] Ocean forcing data file
    - [ ] Control
    - [ ] Scenario XXX
    - [ ] ... to be completed
- [ ] Benoit's grid (first)
    - [ ] grid with ice tongue and instruction to do it in README or elsewhere
    - [ ] tools and instruction to convert restart from 1 grid to another that also manage change in partition numbers
    - [ ] grid refined over all the run history and instruction to do it in README or elsewhere
- [ ] Initial condition compatible with Benoit's grid (first)
    - [ ] Initial condition compatible with grid with no ice tongue
    - [ ] Initial condition compatible with grid refined
    - [ ] Build netcdf restart => this will make obsolete the two previous points.
- [ ] Clean directory for only the production file

### SIF
- [ ] Check SSA  solver
- [ ] Check H    solver
- [ ] Check XIOS solver
    - [X] Check metadata
    - [X] Check integrated value wrt to save scalar
    - [ ] Clarify with only Melt give an integrated value different in Savescalar
    - [ ] Check all netcdf varaible with paraview (need sno account)
- [X] Netcdf restart
    - [X] Tools to interpolate from a result file to a netcdf file has been done
    - [X] Differences on domain interior tiny in a comparison .result and .nc restart format
    - [X] Why is there differences on the domain edges ?
- [ ] Coulomb Regularisé
    - [ ] Convert linear friction coef from initial condition to CR friction coefficient
    - [ ] Check all blocs
- [ ] Check SMB solver
    - [ ] Read reference SMB
    - [ ] Read anomaly SMB year by year
    - [ ] Compute SMB as it is required for ISMIP6
- [ ] Check BMB solver
    - [ ] PICO melt solver blocs and parameter
    - [ ] Quadratic local melt blocs and parameter
    - [ ] Check SIF compatibility (not both melt activated ...)
    - [ ] Compute BMB as it is required for ISMIP6
- [ ] Nearestpoint => create basin mask

### Solver
- [ ] Build Quadratic melt Solver.F90
- [ ] Update Benoit's regional pp to new variable name (dh/dt for example)

### Test to do
- [ ] Check Benoit's IMITMIP diag and XIOS integrated method
   - [ ] check all ouput with paraview (min, max, pattern ...)
   - [ ] when succesful rm Benoit's diag
   - [ ] buil pp script to do it
- [ ] compare overall with Benoit's configuration

### Run to do
#### Reference run
- [ ] Control
#### Perturbation runs
- [ ] ... to be completed
#### Otional runs
- [ ]

### Post Processing
- [ ] review the ISMIP6 XIOS context file
- [ ] review the ISMIP6 XIOS output file
- [ ] build the pp script to convert XIOS output to ISMIP6 output
- [ ] review the ISMIP6 output
- [ ] process all the output

### Distribution
- [ ] push data on the server
