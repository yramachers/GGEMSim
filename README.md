# GGEMSim
Simulate electron transport in a GGEM configuration.

Build depends on ROOT: easiest build environment at Warwick, source from cvmfs.

source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc13-opt/setup.sh

which delivers ROOT version 6.30, built with C++20.

The present code is missing essential resource files. Needs realistic geometry GDML files, 
consistent with static electric field map from COMSOL. Both are absent. Examples for coding of GDML
files using ROOT are in folder utils. A cross section test executable is also present for electron-Argon atom
scattering at the relevant energies. The data folder
contains the output GDML files from the examples. Not realistic, and missing the field maps hence are 
pure placeholders for now.
