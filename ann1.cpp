//----------------------------------------------------------------------
//		File:			ann_sample.cpp
//		Programmer:		Sunil Arya and David Mount
//		Last modified:	03/04/98 (Release 0.1)
//		Description:	Sample program for ANN
//----------------------------------------------------------------------
// Copyright (c) 1997-2005 University of Maryland and Sunil Arya and
// David Mount.  All Rights Reserved.
// 
// This software and related documentation is part of the Approximate
// Nearest Neighbor Library (ANN).  This software is provided under
// the provisions of the Lesser GNU Public License (LGPL).  See the
// file ../ReadMe.txt for further information.
// 
// The University of Maryland (U.M.) and the authors make no
// representations about the suitability or fitness of this software for
// any purpose.  It is provided "as is" without express or implied
// warranty.
//----------------------------------------------------------------------

#include <cstdlib>						// C standard library
#include <cstdio>						// C I/O (for sscanf)
#include <cstring>						// string manipulation
#include <iostream>					// file I/O
#include <fstream>
#include <ANN/ANN.h>				// ANN declarations
#include <sstream> 

using namespace std;					// make std:: accessible

//----------------------------------------------------------------------
// ann_sample
//
// This is a simple sample program for the ANN library.	 After compiling,
// it can be run as follows.
// 
// ann_sample [-d dim] [-max mpts] [-nn k] [-e eps] [-df data] [-qf query]
//
// where
//		dim				is the dimension of the space (default = 2)
//		mpts			maximum number of data points (default = 1000)
//		k				number of nearest neighbors per query (default 1)
//		eps				is the error bound (default = 0.0)
//		data			file containing data points
//		query			file containing query points
//
// Results are sent to the standard output.
//----------------------------------------------------------------------

//----------------------------------------------------------------------
//	Parameters that are set in getArgs()
//----------------------------------------------------------------------
void getArgs(int argc, char **argv);			// get command-line arguments

int				k				= 1;			// number of nearest neighbors
int				dim				= 2;			// dimension
double			eps				= 0;			// error bound
int				maxPts			= 1000;			// maximum number of data points

istream*		dataIn			= NULL;			// input for data points
istream*		queryIn			= NULL;			// input for query points

bool readPt(istream &in, ANNpoint p)			// read point (false on EOF)
{
	char line[1028];
	int check = in.getline(line, 1028);
	if (!check) {
		return false;
	}
	istringstream ss(line);
	for (int i = 0; i < dim; i++){
		if (!(ss >> p[i])){
			return false;
		}
	}
	return true;
}

void printPt(ostream &out, ANNpoint p)			// print point
{
	out << "(" << p[0];
	for (int i = 1; i < dim; i++) {
		out << ", " << p[i];
	}
	out << ")\n";
}

int main(int argc, char **argv)
{
	int					nPts;					// actual number of data points
	ANNpointArray		dataPts;				// data points
	ANNpoint			queryPt;				// query point
	ANNidxArray			nnIdx;					// near neighbor indices
	ANNdistArray		dists;					// near neighbor distances
	ANNkd_tree*			kdTree;					// search structure

	getArgs(argc, argv);						// read command-line arguments

	queryPt = annAllocPt(dim);					// allocate query point
	dataPts = annAllocPts(maxPts, dim);			// allocate data points
	nnIdx = new ANNidx[k];						// allocate near neigh indices
	dists = new ANNdist[k];						// allocate near neighbor dists

	nPts = 0;									// read data points

	cout << "Data Points:\n";
	while (nPts < maxPts && readPt(*dataIn, dataPts[nPts])) {
		printPt(cout, dataPts[nPts]);
		nPts++;
	}

	kdTree = new ANNkd_tree(					// build search structure
					dataPts,					// the data points
					nPts,						// number of points
					dim);						// dimension of space

	while (readPt(*queryIn, queryPt)) {			// read query points
		cout << "Query point: ";				// echo query point
		printPt(cout, queryPt);

		kdTree->annkSearch(						// search
				queryPt,						// query point
				k,								// number of near neighbors
				nnIdx,							// nearest neighbors (returned)
				dists,							// distance (returned)
				eps);							// error bound

		cout << "\tNN:\tIndex\tDistance\n";
		for (int i = 0; i < k; i++) {			// print summary
			dists[i] = sqrt(dists[i]);			// unsquare distance
			cout << "\t" << i << "\t" << nnIdx[i] << "\t" << dists[i] << "\n";
		}
	}
    delete [] nnIdx;							// clean things up
    delete [] dists;
    delete kdTree;
	annClose();									// done with ANN

	return EXIT_SUCCESS;
}
