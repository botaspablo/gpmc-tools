// gPMC2XiO: reorientates and scales a gPMC cube file. Depending on the arguments, it can also change the endian from little to big and the type from float32 to uint32.
// The input dose should be in MeV/g and the output is in cGy. The output is also scaled by a 1.1 RBE, the number of fractions and of course the spot factor used for the simulation.
// The transformation needed to change the geometry from a gPMC file to XiO is the following:
// (x, y, z) -> (x, -y, -z)

// Usage: gPMC2XiO <output> <input> nx ny nz spotFactor nFractions

#include <cstdlib>
#include <cstring>
#include <stdint.h>
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

uint32_t SwapEndian(uint32_t word);

int main(int argc, char* argv[]) {
	
	if (argc != 8) {
		cerr << "ERROR! Wrong arguments" << endl;
		cerr << "Usage: gpmc2xio <output> <input> nx ny nz spotFactor nFractions" << endl;
		exit(EXIT_FAILURE);
	}
	
	ifstream fileIn;
	ofstream fileOut;
	
	string fout = argv[1];
	string fin  = argv[2];
	uint32_t nX = atoi( argv[3] );
	uint32_t nY = atoi( argv[4] );
	uint32_t nZ = atoi( argv[5] );
	uint32_t nVoxels = nX * nY * nZ;
	
	// Open input Cube
	fileIn.open(fin.c_str(), ios::in | ios::binary);
	if (!fileIn.is_open()) {
		cerr << "Cannot open file: " << fin << endl;
		exit(EXIT_FAILURE);
	}
	cout << "Reading input file ..." << endl;
	vector<float> Cube(nVoxels,0);
	fileIn.read((char*)&Cube[0], nVoxels * sizeof(float));
	fileIn.close();
	
	cout << "Normalization factors:" << endl;
	float spotFactor = atof( argv[6] );
	float nFractions = atof( argv[7] );
	
	float MeVg2cGy = 1.6022e-8;
	float RBE = 1.1;
	float PrescDose = nFractions * RBE * (1e9 / spotFactor);
	
	cout << "\tspotFactor = " << spotFactor << endl;
	cout << "\tnFractions = " << nFractions << endl;
	cout << "\tPrescDose  = " << PrescDose  << endl;
	
	vector<uint32_t> XiOCube(nVoxels,0);
	
	// Flip y->-y, z->-z
	uint32_t overflows=0;
	uint32_t underflows=0;
	cout << "Flipping y and z ..." << endl;
	for (int ix = 0; ix < nX; ix++)
	{
		for (int iy = 0; iy < nY; iy++)
		{
			for (int iz = 0; iz < nZ; iz++)
			{
				float value = PrescDose * MeVg2cGy * Cube.at(ix+iy*nX+iz*nX*nY);
				if (value > 0xFFFFFFFF-0.5)
				{
					overflows++;
					value = 0xFFFFFFFF-0.5;
				}
				if (value < 0x00000000)
				{
					underflows++;
					value = 0x00000000;
				}
				XiOCube[ix+(nY-iy-1)*nX+(nZ-iz-1)*nX*nY] = SwapEndian( (uint32_t)(value+0.5) );
			}
		}
	}
	
	if(overflows>0 || underflows>0)
		cerr << "\tFile: " << fin << endl;
	if(overflows>0)
		cerr << "\tOverflows of the uint32_t container: " << overflows << endl;
	if(underflows>0)
		cerr << "\tUnderflows of the uint32_t container: " << underflows << endl;
	
	// Write out
	fileOut.open(fout.c_str(), ios::out | ios::binary);
	if (!fileOut.is_open()) {
		cerr << "Cannot open output file: " << fout << endl;
		exit(EXIT_FAILURE);
	}
	cout << "Writing out ..." << endl;
	fileOut.write((char*)&XiOCube[0], nVoxels * sizeof(uint32_t));
	fileOut.close();
	
	exit(EXIT_SUCCESS);
}


uint32_t SwapEndian(uint32_t word)
{
	return ((word >> 24)&0xff) | // move byte 3 to byte 0
		((word >> 8)&0xff00)   | // move byte 2 to byte 1
		((word << 8)&0xff0000) | // move byte 1 to byte 2
		((word << 24)&0xff000000); // byte 0 to byte 3
}

