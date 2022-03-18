# Drift times
The analysis is mainly based on Z vs drift times. We can make a selection of cathode-anode or anode-cathode-anode tracks simly based on drift time distribution. The code deltaT_plots.C makes a distribution of drift times for positive drift (deltaT_pos), negative drift(deltaT_neg) as well as anode-cathode-anode (lablelled deltaT_after_stitching). However, I have not set a strict cut on deltaT for selecting anode-cathode-anode tracks as that may limit statistics. I developed additional techniques for anode-cathode-anode selection which I have implemented on the code. 

# driftvelocity
# code for calculating the driftvelocity only at the anode
Spatial distortion calculation requires large number of events to get a reasonable coverage (I will recommend a minimum of 1e6 events). If we don't have ~1 million events we can still estimate the drift velocity at the anode, considering there is negligible spatial distortion at the anode. Anything above 100k events will give a reasonable estimate of drift velocity at the anode. Please, remember the farther you go from the anode towards cathode there will be more spatial distortion and the drift velocity measured without correcting for spatial distortion in Z will have more and more error.

For calculating drift velocity only closest to anode use, velocityonly.C and the graph med_vel will give the median velocity. I am also saving 1D distribution of velocities in differnt bin which can be used to find mean or fitted mean drift velocity easily.



# Codes for calculating electron drift velocity throughout the drift length using anode-cathode-anode crossing tracks using ProtoDUNE-SP data
STEP1: Use the module velocity_module.cc (using velocity.fcl) to generate ntuple
STEP2: Use the ntuple to generate Spatial distortion map using data_crossing.C [very high statistics required]
use binvalues.C to generate a continuous (across cathode) distortion map combining 3D distortion map and distortion at cathode
STEP3: Use spatial distotion map to correct Y and Z distortion and ulimate calculate drift velocity using velocity_code.C.




