Using: process_pdb.py
	
	process_pdb.py was written for Python 2.7 and requires numpy.

	The software reads a data file and assumes that structures (locations and types)
	are annotated in the .pdb file.  It also assumes that the .pdb file contains the
	geometry of a single frame of a trajectory.  NOTE: some .pdb files from the Protein
	Data Bank come with several models in a single file.  If this is the case, then either
	the data file should be modified to contain a single protein, or the software modified
	to read only the first (for example) model in the file.

	Example usage:

	./process_pdb.py /path/to/data/file.pdb

	A sample .pdb file (1BDD) is included in the data directory.

	Running:

	./process_pdb.py ../data/1bdd.pdb

	produces:

	1bdd HELIX 10 19 10 564.681144906 207.759025 207.759025 1 2.21796231574
	1bdd HELIX 25 37 13 998.628959166 253.368147 253.368147 1 1.87528529344
	1bdd HELIX 42 55 14 1396.80643478 419.957706 419.957706 1 2.90549859231

	as output.

	Each line contains the following information:

	1bdd -> name of protein (used from last 4 characters before .pdb in file name)
	HELIX -> type of structure (from .pdb file) either HELIX or SHEET
	10 19 10 -> first, last, length: first amino acid number (1 indexed), last amino acid number, and the length of the segment (number of amino acids)
	564.681144906 -> Largest eigenvalue of the distance matrix built from CA atoms 10-19.
	207.759025 -> Maximum distance in the distance matrix built from CA atoms 10-19.
	207.759025 (2nd occurrence) -> Distance between first and last CA atoms in the segment (CA atoms 10 and 19).
	1 -> Type of helix (indicator of RH alpha helix, other options found here: http://plato.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html
		(This could include things like RH or LH, alpha, pi, omega, etc... RH alpha helices appear to be, by far, the most common helix type).
	2.21796231574 -> Angle (in radians) formed by the first, middle, and last CA atom in the segment.
		(This was used to determine the "straightness" of a segment.  I.e. the closer to 3.14159 the straighter the segment.)
