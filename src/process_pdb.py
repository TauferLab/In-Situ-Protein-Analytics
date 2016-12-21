#!/usr/bin/python

# Copyright (c)
# 	2016 by The University of Delaware
# 	Contributors: Travis Johnston
# 	Affiliation: Global Computing Laboratory, Michela Taufer PI
# 	Url: http://gcl.cis.udel.edu/, https://github.com/TauferLab
# 
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without modification,
# are permitted provided that the following conditions are met:
# 
# 	1. Redistributions of source code must retain the above copyright notice, 
# 	this list of conditions and the following disclaimer.
# 
# 	2. Redistributions in binary form must reproduce the above copyright notice,
# 	this list of conditions and the following disclaimer in the documentation
# 	and/or other materials provided with the distribution.
# 
# 	3. If this code is used to create a published work then the	following paper
# 	must be cited:
# 
# 		T. Johnston, B. Zhang, A. Liwo, S. Crivelli, and M. Taufer.
# 		In-Situ Data Analytics and Indexing of Protein Trajectories.
# 		Journal of Computational Chemistry (JCC), 2017.
# 
# 	4.  Permission of the PI must be obtained before this software is used
# 	for commercial purposes.  (Contact: taufer@acm.org)
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
# IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
# INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
# OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
# OF THE POSSIBILITY OF SUCH DAMAGE.



import sys
import numpy as np
from numpy import linalg as npALG
import math


def compute_ev_from_range(protein, structure_info, atoms, structure_type):

	X = []
	for atom in atoms:
		if atom[0]==structure_info[0] and structure_info[1]<= atom[1] and atom[1]<=structure_info[2]:
			X.append( (atom[2], atom[3], atom[4]) )
	if len(X) > 0:	
		M = np.zeros( (len(X), len(X)) )

		for i in xrange(len(X)):
			for j in xrange(i+1, len(X)):
				M[i, j] = (X[i][0]-X[j][0])**2 + (X[i][1]-X[j][1])**2 + (X[i][2]-X[j][2])**2
				M[j, i] = M[i, j]

		ev = npALG.eigvalsh(M)[-1]

		### Compute angle between first, middle, and last atoms.
		### This angle is used as a metric of "straightness" of the structure and is intended mostly for use with the sheets.
		
		if len(X) < 3:
			angle=math.pi
		else:
			u = np.zeros(3)		### Point from middle atom to first atom
			v = np.zeros(3)		### Point from middle atom to last atom
			mIdx = len(X)/2		### Integer arithmetic rounds down... this is ok.
			
			u[0] = X[0][0] - X[mIdx][0]
			u[1] = X[0][1] - X[mIdx][1]
			u[2] = X[0][2] - X[mIdx][2]

			v[0] = X[-1][0] - X[mIdx][0]
			v[1] = X[-1][1] - X[mIdx][1]
			v[2] = X[-1][2] - X[mIdx][2]

			angle = math.acos( np.dot(u, v)/(np.dot(u, u)*np.dot(v, v))**.5 )

		print protein, structure_type, structure_info[1], structure_info[2], len(X), ev, M.max(), M[0, len(X)-1], structure_info[3], angle

	return None




with open(sys.argv[1],'r') as inFile:
	Structures = [ [], [] ]		### Structures[0] is a list of info about the helices, Structures[1] is a list of info about the strands.
								### Each list will contain tuples (chain#, first_amino_acid#, last_amino_acid#)
	CA = []						### A list of information about the CA atoms.  This will be a list of tuples, (chain#, amino_acid#, x, y, z)

	for line in inFile:
		if line[:5]=="HELIX":
			assert line[19]==line[31]
			chain = line[19]
			first = int(line[21:25])
			last = int(line[33:37])
			type_of_helix = int(line[38:40])
			Structures[0].append( (chain, first, last, type_of_helix) )

		elif line[:5]=="SHEET":
			assert line[21]==line[32]	### amino acids on same chain.
			chain = line[21]			
			first = int(line[22:26])
			last = int(line[33:37])
			strand_sense = int(line[38:40])
			Structures[1].append( (chain, first, last, strand_sense) )

		elif line[:5]=="MODEL":			### Done finding structures, now we are looking at the first model in the .pdb file.
			find_structure = False
			find_atoms = True
		
		if line[:4]=="ATOM":
			if line[12:16]==" CA ":
				alternate_location = line[16]
				if alternate_location==' ' or alternate_location=='A':
					CA.append( (line[21], int(line[22:26]), float(line[30:38]), float(line[38:46]), float(line[46:54])) )
					
		elif line[:6]=="HETATM":
			if line[12:16]==" CA ":
				alternate_location = line[16]
				if alternate_location==' ' or alternate_location=='A':
					CA.append( (line[2], int(line[22:26]), float(line[30:38]), float(line[38:46]), float(line[46:54])) )
				
Deduplicated_Structures = [ [], []]

for helix in xrange(len(Structures[0])):
	### if this helix is not a superset of another structure, then add it to the Deduplicated list.
	Good = True
	for compare in xrange(len(Structures[0])):
		if compare!=helix:		### make sure we are not looking at the same index.
			a = Structures[0][helix]
			b = Structures[0][compare]
			if a[0]==b[0] and a[1]<=b[1] and b[2]<=a[2]:
				### a and b are on the same chain, and b is a subset of a
				if a[1] < b[1] or b[2] < a[2]:
					### b is a proper subset of a, DO NOT KEEP a.
					Good = False
				else:
					if helix > compare:
						Good = False	### if they have exactly the same range, only keep the one with smaller index.
	
	for compare in xrange(len(Structures[1])):
		a = Structures[0][helix]
		b = Structures[1][compare]
		if a[0]==b[0] and a[1]<=b[1] and b[2]<=a[2]:
			### a and b are on the same chain, and b is a subset of a
			Good = False
						
	if Good:
		Deduplicated_Structures[0].append( Structures[0][helix] )


for strand in xrange(len(Structures[1])):
	### if this strand is not a superset of another structure, then add it to the Deduplicated list.
	Good = True
	for compare in xrange(len(Structures[1])):
		if compare!=strand:		### make sure we are not looking at the same index.
			a = Structures[1][strand]
			b = Structures[1][compare]
			if a[0]==b[0] and a[1]<=b[1] and b[2]<=a[2]:
				### a and b are on the same chain, and b is a subset of a
				if a[1] < b[1] or b[2] < a[2]:
					### b is a proper subset of a, DO NOT KEEP a.
					Good = False
				else:
					if strand > compare:
						Good = False	### if they have exactly the same range, only keep the one with smaller index.
	
	for compare in xrange(len(Structures[0])):
		a = Structures[1][strand]
		b = Structures[0][compare]
		if a[0]==b[0] and a[1]<=b[1] and b[2]<=a[2]:
			### a and b are on the same chain, and b is a subset of a
			Good = False
						
	if Good:
		Deduplicated_Structures[1].append( Structures[1][strand] )


protein = sys.argv[1][-8:-4]

for helix in Deduplicated_Structures[0]:
	compute_ev_from_range(protein, helix, CA, "HELIX")

for strand in Deduplicated_Structures[1]:
	compute_ev_from_range(protein, strand, CA, "SHEET")
