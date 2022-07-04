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
    - [ ] grid refined over all the run history and instruction to do it in README or elsewhere
- [ ] Initial condition compatible with Benoit's grid (first)
    - [ ] Initial condition compatible with grid with no ice tongue
    - [ ] Initial condition compatible with grid refined
- [ ] Clean directory for only the production file

### SIF
- [ ] Check SSA  solver
- [ ] Check H    solver
- [ ] Check XIOS solver
- [ ] Check SMB solver
    - [ ] Read reference SMB
    - [ ] Read anomaly SMB year by year
    - [ ] Compute SMB as it is required for ISMIP6
- [ ] Check BMB solver
    - [ ] PICO melt solver blocs
    - [ ] Quadratic local melt blocs
    - [ ] Check SIF compatibility (not both melt activated ...)
    - [ ] Compute BMB as it is required for ISMIP6
- [ ] Nearestpoint => create basin mask

### Solver
- [ ] 

### Test to do
- [ ] Check Benoit's IMITMIP diag and XIOS integrated method
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
