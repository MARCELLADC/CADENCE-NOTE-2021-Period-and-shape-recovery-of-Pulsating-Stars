# CADENCE-NOTE-2021-Period-and-shape-recovery-of-Pulsating-Stars


A  notebook that use known and our metrics written for the cadence Note 2021.
The input of this metric is a theoretical light-curve or a pulsating variable template (like RRab.cvs), that is convoluted with the chosen survey strategy (OpSim) to obtain the ”observed” light curves in all the filters $ugrizY$ with the expected S/N for each phase point, 
filtered to exclude the saturated points for that specific OpSim and analyzed by the multiband gatspy package to obtain single filter and median periodograms.
the fit of the observed light curve in each band through the superposition of nHarm are given
From this analysis one can derive the number of visits, and the number of ''useful''visits, the dimension of the biggest hole in the phased light curve, the recovered  period and its accuracy,  mean magnitude (averaged in intensity) and amplitude in each filter. 
The differences of these quantities respect to those of the input light curve (Delta mean mag and Delta Amplitude) in each bands are additional very useful FoMs. 

Comments are welcome
