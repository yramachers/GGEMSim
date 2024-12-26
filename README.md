# GGEMSim
Simulate electron transport in a GGEM configuration.

Build depends on ROOT: easiest build environment at Warwick, source from cvmfs.

source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc13-opt/setup.sh

which delivers ROOT version 6.30, built with C++20.

The data folder contains GDML file(s) consistent with COMSOL geometries but in 3D. COMSOL field calculation
outputs in 2D are held in ROOT files, converted with a Python script from COMSOL output files.

Currently builds two executables: single_charge.exe transports a single charge from a configurable
starting point and stores all interactions points to file and returns the number of charges reaching the anode 
in the geometry, identified by GDML name.

The charge_scan executable counts events of interest, and outputs them to file. That is meant to enable
larger scale simulations.
