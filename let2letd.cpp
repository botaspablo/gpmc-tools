// LET cube to DOSExLET
// WARNING: dose and LET cubes are supposed to have the same endianness (big)!!

#include <cstdlib>
#include <cstring>
#include <stdint.h>
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

uint32_t SwapEndian(uint32_t word);

int main(int argc, char* argv[])
{
	
	if( argc != 4 )
	{
		cout << "ERROR! wrong arguments" << endl;
		cout << "Usage: let2letd <LETd> <LET> <Dose>" << endl;
		exit(EXIT_FAILURE);
	}
	
	string fout = argv[1];
	string fin = argv[2];
	string fdose = argv[3];
	
	ifstream strin;
	ofstream strout;
	ifstream strdose;
	
	// Open input LET
	strin.open(fin.c_str(), ios::in | ios::binary | ios::ate);
	if (!strin.is_open()) {
		cerr << "Cannot open file: " << fin << endl;
		exit(EXIT_FAILURE);
	}
	
	uint32_t nVoxels = strin.tellg() / sizeof(uint32_t);
	strin.seekg(0, ios::beg);
	
	cout << "Reading input LET ..." << endl;
	vector<uint32_t> LET(nVoxels);
	strin.read((char*)&LET[0], nVoxels * sizeof(uint32_t));
	strin.close();
	
	// Open input Dose
	strdose.open(fdose.c_str(), ios::in | ios::binary);
	if (!strdose.is_open()) {
		cerr << "Cannot open file: " << fdose << endl;
		exit(EXIT_FAILURE);
	}
	cout << "Reading input Dose ..." << endl;
	vector<uint32_t> dose(nVoxels);
	strdose.read((char*)&dose[0], nVoxels * 4);
	strdose.close();
	
	cout << "Calculating LETd ..." << endl;
	vector<uint32_t> out(nVoxels,0);
	
	int overflows=0;
	for(uint32_t i=0; i<nVoxels; i++)
	{
		// change big endian to little to perform the multiplication
		float d = SwapEndian((uint32_t)(dose[i]));
		float l = SwapEndian((uint32_t)(LET[i]));
		uint64_t value = d * l;
		if (value >= 0xFFFFFFFF-0.5)
		{
			// overflows++;
			value = 0xFFFFFFFF-0.5;
		}
		// if (value > 1000000)
		// {
		// 	// This is a possible hotspot: output indexes and coordinates
		// 	unsigned int ix = i%250;
		// 	unsigned int iy = (i/250)%250;
		// 	unsigned int iz = i/250/250;
		// 	cout << value << endl;
		// 	cout << ix << " " << iy << " " << iz << endl;
		// 	cout << "\t" << 2*ix << " " << 2*iy << " " << 2*iz << endl;
		// }
		// store the data in little endian
		out[i] = SwapEndian( (uint32_t)(value+0.5) );
	}
	
	// cout << "\tOverflows of uint32_t: " << overflows << endl;
	
	cout << "Writing output LETd ..." << endl;
	// Open output LETd
	strout.open(fout.c_str(), ios::out | ios::binary);
	if (!strout.is_open()) {
		cerr << "Cannot open file: " << fout << endl;
		exit(EXIT_FAILURE);
	}
	strout.write((char*)&out[0], nVoxels * sizeof(uint32_t));

	strout.close();
	
	exit(EXIT_SUCCESS);
}

uint32_t SwapEndian(uint32_t word) {
	return ((word >> 24)&0xff) | ((word >> 8)&0xff00) | ((word << 8)&0xff0000) | ((word << 24)&0xff000000);
}
