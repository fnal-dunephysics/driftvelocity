# driftvelocity
Codes for calculating electron drift velocity using anode-cathode-anode crossing tracks using ProtoDUNE-SP data
STEP1: Use the module velocity_module.cc (using velocity.fcl) to generate ntuple
STEP2: Use the ntuple to generate Spatial distortion map using data_crossing.C [very high statistics required]
use binvalues.C to generate a continuous (across cathode) distortion map combining 3D distortion map and distortion at cathode
STEP3: Use spatial distotion map to correct Y and Z distortion and ulimate calculate drift velocity


Spatial distortion calculation requires large number of events to get a reasonable coverage (I will recommend a minimum of 1e6 events). If we don't have ~1 million events we can still estimate the drift velocity at the anode, considering there is negligible spatial distortion at the anode. Anything above 100k events will give a reasonable estimate of drift velocity at the anode. Please, remember the farther you go from the anode towards cathode there will be more spatial distortion and the drift velocity measured without correcting for spatial distortion in Z will have more and more error.
