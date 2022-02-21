# driftvelocity
Codes for calculating electron drift velocity using anode-cathode-anode crossing tracks using ProtoDUNE-SP data
STEP1: Use the module velocity_module.cc (using velocity.fcl) to generate ntuple
STEP2: Use the ntuple to generate Spatial distortion map using data_crossing.C [very high statistics required]
use binvalues.C to generate a continuous (across cathode) distortion map combining 3D distortion map and distortion at cathode
STEP3: Use spatical distotion map to correct Y and Z distortion and ulimate calculate drift velocity
