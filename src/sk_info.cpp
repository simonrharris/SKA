#include <string> //std::string
#include <vector> //std::vector
#include <iostream> //std::cout
#include <fstream> //std::ofstream
#include <math.h>
#include "kmers.hpp"
#include "general.hpp"

using namespace std;



int getKmerFileInfo(const vector<string> & kmerfiles, const bool tabulated){

	vector < string > sampleNames;
	int kmersize;
	ifstream fileStream;

	if (tabulated){
		cout << "File\tKmer size\t# samples\t#kmers" << endl;
	}

	for (int s = 0; s < kmerfiles.size(); ++s){

		
		if (tabulated){
			cout << splitFileName(kmerfiles[s]);
		}
		else {
			cout << endl << splitFileName(kmerfiles[s]) << endl;
			cout << string(kmerfiles[s].size(), '=') << endl;
		}

		fileStream.open(kmerfiles[s], ios::in);

		if (fileStream.fail()) {
			cout << "Failed to open " << kmerfiles[s] << "\n" << endl;
			return 0;
		}
		
		try {
			int returnval = readKmerHeader(fileStream, kmersize, sampleNames);
		}
		catch (int e){
			cout << "An exception occurred when reading file " << kmerfiles[s] << ". Please check the format. Exception Nr. " << e << '\n';
			return 1;
		}

		if (tabulated){
			cout << "\t" << kmersize << "\t" << sampleNames.size();
		}
		else {
			cout << "Kmer size: " << kmersize << endl;
			cout << "Number of samples: " << sampleNames.size() << endl;
		}

		int kmercount=0;
		char basebuffer[1];
		char kmerbuffer[(kmersize*2/3)+1];
		char asciibuffer[int(ceil(float(sampleNames.size())/6))];

		while (fileStream.read(asciibuffer, sizeof(asciibuffer))){//read the ascii representation of the taxa
			string asciibits (asciibuffer, sizeof(asciibuffer));
			vector < bool > mybits;
			vectorbool_from_ascii(asciibits, mybits);//read the ascii representation of the taxa to a vector of bools

			while (fileStream.peek()!='\n' && fileStream.read(kmerbuffer, sizeof(kmerbuffer))){
				kmercount++;
			}
		}
		if (tabulated){
			cout << "\t" << kmercount << endl;
		}
		else {
			cout << "Number of kmers: " << kmercount << endl << endl;
		}
		
		fileStream.close();
	}
	return 0;
}