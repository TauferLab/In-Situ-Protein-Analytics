#In-Situ-Protein-Analytics

This repository contains tools used to perform an eigenvalue-based, 
in-situ data analytics of protein trajectories following the approach detailed in:

	T. Johnston, B. Zhang, A. Liwo, S. Crivelli, and M. Taufer.  In-Situ Data Analytics and Indexing of Protein Trajectories. Journal of Computational Chemistry (JCC), 2017.

This software targets MD simulations at the exascale and proposes a novel
technique for *in situ* data analysis and indexing of MD
trajectories. Our technique maps individual trajectories'
substructures (i.e., alpha-helices, beta-strands) to metadata
frame by frame. The metadata captures the conformational properties of
the substructures. The ensemble of metadata can be used for automatic,
strategic analysis within a trajectory or across trajectories,
without manually identify those portions of trajectories in which
critical changes take place. We demonstrate our technique's
effectiveness by applying it to 26.3k helices and 31.2k strands from
9,917 PDB proteins and by providing three empirical case studies.
