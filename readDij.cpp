// readDij: prints a Dij (or LETij) influence matrix and outputs its header and the header of every spot contained in it.
//          if -extract is specified as first parameter, it'll extract the matrices.
//          The matrix format is float little endian

// Usage: readDij -extract <dij>

#include <cstdlib>
#include <cstring>
#include <stdint.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>

using namespace std;

// define the header
typedef struct {
	float gantry;
	float table;
	float collimator;
	// in mm
	float   spotSizeX; // not needed for protons
	float   spotSizeY; // not needed for protons
	float   voxSizeX;
	float   voxSizeY;
	float   voxSizeZ;
	int32_t nx;
	int32_t ny;
	int32_t nz;
	int32_t nSpots;
	float   factor;
} header_t;

// define the spotheader_t
typedef struct {
	float energy;
	float x;
	float y;
} spotheader_t;

int main(int argc, char* argv[]) {
	
	if (argc != 2 && argc != 3) {
		cerr << "Usage: readDij <dij>" << endl;
		cerr << "Usage with extracty option: readDij -extract <dij>" << endl;
		exit(EXIT_FAILURE);
	}
	
	bool ifExtract=0;
	if (argc == 3 && strcmp(argv[1],"-extract") == 0)
		ifExtract = 1;
	
	string fdij = argv[1+ifExtract];
	ifstream dij;
	dij.open(fdij.c_str(), ios::in | ios::binary);
	if (!dij.is_open()) {
		cerr << "Cannot open input file: " << fdij << endl;
		exit(EXIT_FAILURE);
	}
	
	// Write headerDose
	cout << "Reading header of Dij file ... " << endl;
	header_t header;
	dij.read(reinterpret_cast<char*>(&header), sizeof(header_t));
	cout << "\tgantry:     " << header.gantry     << endl;
	cout << "\ttable:      " << header.table      << endl;
	cout << "\tcollimator: " << header.collimator << endl;
	cout << "\tspotSizeX:  " << header.spotSizeX  << endl;
	cout << "\tspotSizeY:  " << header.spotSizeY  << endl;
	cout << "\tvoxSizeX:   " << header.voxSizeX   << endl;
	cout << "\tvoxSizeY:   " << header.voxSizeY   << endl;
	cout << "\tvoxSizeZ:   " << header.voxSizeZ   << endl;
	cout << "\tnx:         " << header.nx         << endl;
	cout << "\tny:         " << header.ny         << endl;
	cout << "\tnz:         " << header.nz         << endl;
	cout << "\tnSpots:     " << header.nSpots     << endl;
	cout << "\tfactor:     " << header.factor     << endl;
	
	uint32_t nVoxels = header.nx * header.ny * header.nz;
	
	// Write body
	cout << "Reading Dij body ..." << endl;
	
	vector<spotheader_t> spots(header.nSpots);
	vector<int32_t> nonZeroDose(header.nSpots);
	vector< vector<int32_t> > index(header.nSpots);
	vector< vector<uint16_t> > vecDose(header.nSpots);
	
	vector<float> out;
	if(ifExtract) out.resize(nVoxels, 0);
	
	for(size_t iSpot=0; iSpot<header.nSpots; iSpot++)
	{
		dij.read(reinterpret_cast<char*>(&spots.at(iSpot).energy), sizeof(float));
		dij.read(reinterpret_cast<char*>(&spots.at(iSpot).x     ), sizeof(float));
		dij.read(reinterpret_cast<char*>(&spots.at(iSpot).y     ), sizeof(float));
		dij.read(reinterpret_cast<char*>(&nonZeroDose.at(iSpot) ), sizeof(int32_t));

		// Each spot has this mini header. Then we iterate over all the non-zero voxels of the spot's dose matrix
		cout << "Spot Number: " << iSpot << endl;
		cout << "\tenergy:    " << spots.at(iSpot).energy << endl;
		cout << "\tx:         " << spots.at(iSpot).x      << endl;
		cout << "\ty:         " << spots.at(iSpot).y      << endl;
		cout << "\tnon Zero:  " << nonZeroDose.at(iSpot)  << endl;
		
		index.at(iSpot).resize( nonZeroDose.at(iSpot) );
		vecDose.at(iSpot).resize( nonZeroDose.at(iSpot) );
		for(size_t i=0; i<nonZeroDose.at(iSpot); i++)
		{
			dij.read(reinterpret_cast<char*>(&index.at(iSpot).at(i)), sizeof(int32_t));
			dij.read(reinterpret_cast<char*>(&vecDose.at(iSpot).at(i)), sizeof(uint16_t));
			// cout << indexDoseOut.at(iSpot).at(i) << "\t" << vecDose.at(iSpot).at(i) << endl;
			
			if(ifExtract)
				out.at(index.at(iSpot).at(i)) = (float) header.factor * vecDose.at(iSpot).at(i);
				// cout << outXio.at(index.at(iSpot).at(i));
		}
		
		if(ifExtract)
		{
			ostringstream it;
			it << iSpot+1;
			string file = string(5-it.str().length(), '0') + it.str();
			file = "spot." + file + ".bin";
			ofstream strout;
			strout.open(file.c_str(), ios::out | ios::binary);
			if (!strout.is_open()) {
				cerr << "Cannot open output file: " << file << endl;
				exit(EXIT_FAILURE);
			}
			strout.write(reinterpret_cast<char*>(&out[0]), nVoxels*sizeof(uint32_t));
			fill(out.begin(), out.end(), 0);
		}
	
	}
	
	exit(EXIT_SUCCESS);
}



