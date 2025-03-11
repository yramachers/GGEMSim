# GGEMSim
Simulate electron transport in a GGEM configuration.

Build depends on ROOT: easiest build environment at Warwick, source from cvmfs.

source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc13-opt/setup.sh

which delivers ROOT version 6.30, built with C++20.

The data folder contains GDML file(s) consistent with COMSOL geometries but in 3D. COMSOL field calculation
outputs in 2D are held in ROOT files, converted with a Python script from COMSOL output files.

Currently builds two executables: single_charge.exe transports a single charge from a configurable
starting point and stores all interaction points to file and returns the number of charges reaching the anode 
in the geometry, identified by GDML name.

The charge_scan executable counts events of interest, and outputs them to file. That is meant to enable
larger scale simulations. This scans, hard-coded, starting positions in steps of 0.1mm from the default start 
in x and initial energies in steps of 0.3 eV from the default start kinetic energy, 20 iterations in total
for each electron.

Another executable, gain_scan, doesn't scan positions or energies other than the initial given values but 
counts the number of charges passing the middle, z=0, of the double GGEM, in addition to the total number 
of electrons arriving at the anode.
