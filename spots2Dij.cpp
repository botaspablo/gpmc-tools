// doseSpots2DijOpt4D: reads beaminfo file of the beam and all the spots provided. With that information it creates a single array in opt4D format (Jan Unkelbach)
// The beam should be in XiO format (gPMC2XiO)
// The beam should be scored in the CT.
// the input doses should be in either cGy or 100*cGy*keV/mum
// TODO DicomRT uses Gy and Dij are in Gy/gigaproton. At the moment the number of gigaprotons is hardcoded at the beginning of the file. It is also hardcoded in gPMC when we simulate perSpot.

// Usage: doseSpots2Dij <outputDose> <geometry.dat> <beamtramp> <inputs>

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
	float   spotSizeX;  // not needed for protons
	float   spotSizeY;  // not needed for protons
	float   voxSizeX;
	float   voxSizeY;
	float   voxSizeZ;
	int32_t nx;
	int32_t ny;
	int32_t nz;
	int32_t nSpots;
	float   factor;
} header_t;

// define the header
typedef struct {
	float energy;
	float x;
	float y;
	float gp;
} spot_t;

vector<spot_t> readTrampSpots(string file);
void readGeometry(string file, header_t &header);
uint32_t SwapEndian(uint32_t word);

int main(int argc, char* argv[]) {
	
	if (argc < 5) {
		cerr << "ERROR! Wrong arguments" << endl;
		cerr << "Usage 1: spots2Dij <outputDose> <geometry.dat> <beamtramp> <inputs>" << endl;
		cerr << "Usage 2: spots2Dij -thres X <outputDose> <geometry.dat> <beamtramp> <inputs>" << endl;
		cerr << "\tInputs:  Dose or LETd inputs. Dose in cGy and LETd in 100*cGy*keV/mum." << endl;
		cerr << "\tThreshold: A single value in percentage of the peak dose. It will be checked on every voxel." << endl;
		exit(EXIT_FAILURE);
	}
	
	header_t header;
	
	int   ifOpt = 0;
	bool  ifThres = 1;
	float threshold = 0;
	if(strcmp(argv[1],"-thres") == 0)
	{
		ifOpt+=2;
		ifThres=1;
		threshold = atof(argv[2]);
	}
	string fOutDose = argv[1+ifOpt];
	string fGeo = argv[2+ifOpt];
	string fileTramp = argv[3+ifOpt];
	header.nSpots = (argc - (4+ifOpt));
	vector<string> fInDose(header.nSpots);
	for(size_t i=0; i<header.nSpots; i++)
		fInDose[i] = argv[i+(4+ifOpt)];
	
	float const cGy2Gy=0.01;
	float const nGp=100;
	
	cout << "Number of input spots: " << header.nSpots << endl;
	cout << "Output: " << fOutDose << endl;
	// cout << "Reading geometry ..." << endl;
	readGeometry(fGeo, header);
	
	uint64_t nVoxels = header.nx * header.ny * header.nz;
	header.gantry = 270 - header.gantry; // it needs the beaminfo value, not the one in PATIENT.in
	header.collimator = 0;
	header.spotSizeX = 1;
	header.spotSizeY = 1;
	
	// read tramp file
	// cout << "Reading tramp file ... " << endl;
	vector<spot_t> spots = readTrampSpots(fileTramp);
	// cout << "voxSizeX: " << header.voxSizeX << endl;
	// cout << "voxSizeY: " << header.voxSizeY << endl;
	// cout << "voxSizeZ: " << header.voxSizeZ << endl;
	// cout << "nx: " << header.nx << endl;
	// cout << "ny: " << header.ny << endl;
	// cout << "nz: " << header.nz << endl;
	// cout << "nSpots: " << header.nSpots << endl;
	
	cout << "Getting maxValue, total value per beamlet and absolute total value ... " << endl;
	double maxValue=0;
	double totalDose=0;
	vector<float> beamletDose(header.nSpots,0);
	uint32_t* XiOCube = new uint32_t[nVoxels];
	for(size_t iSpot=0; iSpot<header.nSpots; iSpot++)
	{
		ifstream strIn;
		strIn.open(fInDose[iSpot].c_str(), ios::in | ios::binary);
		if (!strIn.is_open()) {
			cerr << "Cannot open file: " << fInDose[iSpot] << endl;
			exit(EXIT_FAILURE);
		}
		// cout << "Reading input file " << fInDose.at(iSpot) << endl;
		strIn.read((char*)&XiOCube[0], nVoxels * sizeof(uint32_t));
		strIn.close();
		
		for(size_t i=0; i < nVoxels; i++)
		{
			double value = SwapEndian(XiOCube[i]) * cGy2Gy/nGp;
			totalDose += value;
			beamletDose[iSpot] += value;
			if( maxValue < value )
			{
				maxValue = value;
				if(SwapEndian(XiOCube[i]) == 0xFFFFFFFF)
				{
					cerr << "ERROR! Probable INF value (only first is detected) in data." << endl;
					cerr << "Spot number " << iSpot+1 << " " << maxValue << endl;
					cerr << "\tValue " << SwapEndian(XiOCube[i]) << endl;
				}
			}
		}
	}
	// cout << endl;
	// cout << endl;
	
	// Set dose factor
	header.factor = maxValue / 0x7FFF;
	cout << "\tgigaprotons:     " << nGp << endl;
	cout << "\tmaxValue:        " << maxValue << "  (Gy/nGp)" << endl;
	cout << "\tfactor:          " << header.factor << "  (maxValue/nGp / 0x7FFF)" << endl;
	float effthres = 0;
	if(ifThres)
	{
		effthres = threshold/100*maxValue;
		cout << "\tthreshold:       " << threshold << "  (\%)" << endl;
		cout << "\tthreshold abs:   " << effthres << "  (Gy/nGp)" << endl;
	}
	
	
	// Get useful indexes after normalization
	vector< vector<int32_t> > indexDoseOut(header.nSpots);
	vector< vector<uint16_t> > vecDose(header.nSpots);
	vector<int32_t> nonZeroDose(header.nSpots, 0);
	// cout << "Getting non-zero voxels, setting input and ouput index arrays and number of non-zero voxels ... " << endl;
	for(size_t iSpot=0; iSpot<header.nSpots; iSpot++)
	{
		// cout << 100*iSpot/header.nSpots << "\%  " << flush;
		ifstream strIn;
		strIn.open(fInDose[iSpot].c_str(), ios::in | ios::binary);
		strIn.read((char*)&XiOCube[0], nVoxels * sizeof(uint32_t));
		strIn.close();
		
		// Get non-zero indexes and the value in those indexes.
		for (int ix = 0; ix < header.nx; ix++)
		{
			for (int iy = 0; iy < header.ny; iy++)
			{
				for (int iz = 0; iz < header.nz; iz++)
				{
					int indexIn = ix+iy*header.nx+iz*header.nx*header.ny;
					
					float temp = (float)(SwapEndian(XiOCube[indexIn])) * cGy2Gy/nGp;
					if(temp > effthres)
					{
						uint16_t value = (uint16_t)( temp/header.factor + 0.5 );
						
						// OUTPUT INDEX
						// XYZ and flip everything
						// int indexOut = (header.nx-ix-1)+(header.ny-iy-1)*header.nx+(header.nz-iz-1)*header.nx*header.ny;
						// XZY Swap Y-Z and flip everything
						// int32_t indexOut = (header.nx-ix-1)+(header.nz-iz-1)*header.nx+(header.ny-iy-1)*header.nx*header.nz;
						// XZY Swap Y-Z
						int32_t indexOut = ix + iz*header.nx + iy*header.nx*header.nz;
					
						indexDoseOut[iSpot].push_back( indexOut );
						vecDose[iSpot].push_back( value );
						nonZeroDose[iSpot]+=1;
					}
				}
			}
		}
	} // spots loop
	
	delete[] XiOCube;

	ofstream strOut;
	strOut.open(fOutDose.c_str(), ios::out | ios::binary);
	if (!strOut.is_open()) {
		cerr << "Cannot open output file: " << fOutDose << endl;
		exit(EXIT_FAILURE);
	}
	
	// swap header for output:
	float tempswap = header.voxSizeY;
	header.voxSizeY = header.voxSizeZ;
	header.voxSizeZ = tempswap;
	
	int tempswap2 = header.ny;
	header.ny = header.nz;
	header.nz = tempswap2;
	
	// Write header
	// cout << "Writing header ..." << endl;
	// cout << "\tgantry: " << header.gantry << endl;
	// cout << "\ttable: " << header.table << endl;
	// cout << "\tcollimator: " << header.collimator << endl;
	// cout << "\tspotSizeX: " << header.spotSizeX << endl;
	// cout << "\tspotSizeY: " << header.spotSizeY << endl;
	// cout << "\tOuput voxSizeX: " << header.voxSizeX << endl;
	// cout << "\tOuput voxSizeY: " << header.voxSizeY << endl;
	// cout << "\tOuput voxSizeZ: " << header.voxSizeZ << endl;
	// cout << "\tOuput nx: " << header.nx << endl;
	// cout << "\tOuput ny: " << header.ny << endl;
	// cout << "\tOuput nz: " << header.nz << endl;
	// cout << "\tnSpots: " << header.nSpots << endl;
	// cout << "\tfactor: " << header.factor << endl;

	strOut.write((char*)&header, sizeof(header_t));
	// Write body
	// cout << "Writing body ..." << endl;
	unsigned int totalNonZero=0;
	for(int iSpot=0; iSpot<header.nSpots; iSpot++)
	{
		strOut.write((char*)&spots[iSpot].energy, sizeof(float));
		strOut.write((char*)&spots[iSpot].x     , sizeof(float));
		strOut.write((char*)&spots[iSpot].y     , sizeof(float));
		strOut.write((char*)&nonZeroDose[iSpot] , sizeof(int32_t));
		// Each spot has this mini header. Then we iterate over all the non-zero voxels of the spot's dose matrix
		// cout << "\tenergy:   " << spots.at(iSpot).energy << endl;
		// cout << "\tx:        " << spots.at(iSpot).x << endl;
		// cout << "\ty:        " << spots.at(iSpot).y << endl;
		// cout << "\tnon Zero: " << nonZeroDose.at(iSpot) << endl;
		totalNonZero += nonZeroDose.at(iSpot);
		
		for(int i=0; i<nonZeroDose[iSpot]; i++)
		{
			strOut.write((char*)&indexDoseOut[iSpot][i], sizeof(int32_t));
			strOut.write((char*)&vecDose[iSpot][i]     , sizeof(uint16_t));
			// cout << indexDoseOut.at(iSpot).at(i) << "\t" << vecDose.at(iSpot).at(i) << endl;
		}
	}
	
	strOut.close();
	
	cout << "\tNon-zero voxels: " << totalNonZero << "  (" << 100.0f*(float)(totalNonZero)/(nVoxels*header.nSpots) << " \%)" << endl;
	
	exit(EXIT_SUCCESS);
}








vector<spot_t> readTrampSpots(string file)
{
	vector<spot_t> spots;
	spot_t thisSpot;
	
	ifstream ifs(file.c_str(), fstream::in);
	if(!ifs.good())
	{
		cerr << "Tramp file " << file << " was not found or not correctly opened." << endl;
		exit(EXIT_FAILURE);
	}
	string line;
	int i=0;
	while ( getline( ifs, line, '\n' ) )
	{
		if( line.empty() || line[0] == '#' ) continue;

		istringstream parser(line);
		string s;
		vector<string> token;
		while ( parser >> s ) token.push_back(s);

		thisSpot.energy = atof( token[0].c_str() );
		thisSpot.x      = atof( token[1].c_str() );
		thisSpot.y      = atof( token[2].c_str() );
		thisSpot.gp     = atof( token[3].c_str() );
		spots.push_back(thisSpot);
		i++;
	}
	ifs.close();
	
	return spots;
}






void readGeometry(string file, header_t &header) {
	
	ifstream ifs(file.c_str(), fstream::in);
	if(!ifs.good())
	{
		cerr << "Configuration file " << file << " was not found or not correctly opened." << endl;
		exit(EXIT_FAILURE);
	}
	
	string line;
	while ( getline( ifs, line, '\n' ) )
	{
		istringstream parser(line);

		string s;
		vector<string> token;
		while ( parser >> s ) token.push_back(s);
		
		if( line.find("gantry") != std::string::npos )
			header.gantry = atof( token.back().c_str() );
		if( line.find("collimator") != std::string::npos )
			header.collimator = atof( token.back().c_str() );
		if( line.find("couch") != std::string::npos )
			header.table = atof( token.back().c_str() );
		if( line.find("sizeVoxelsX") != std::string::npos )
			header.voxSizeX = atof( token.back().c_str() )*10;
		if( line.find("sizeVoxelsY") != std::string::npos )
			header.voxSizeY = atof( token.back().c_str() )*10;
		if( line.find("sizeVoxelsZ") != std::string::npos )
			header.voxSizeZ = atof( token.back().c_str() )*10;
		if( line.find("nVoxelsX") != std::string::npos )
			header.nx = atoi( token.back().c_str() );
		if( line.find("nVoxelsY") != std::string::npos )
			header.ny = atoi( token.back().c_str() );
		if( line.find("nVoxelsZ") != std::string::npos )
			header.nz = atoi( token.back().c_str() );
	}
	
	ifs.close();
}





uint32_t SwapEndian(uint32_t word)
{
	return ((word >> 24)&0xff) | // move byte 3 to byte 0
		((word >> 8)&0xff00)   | // move byte 2 to byte 1
		((word << 8)&0xff0000) | // move byte 1 to byte 2
		((word << 24)&0xff000000); // byte 0 to byte 3
}


