#include <string> //std::string
#include <vector> //std::vector
#include <iostream> //std::cout
#include <fstream> //std::ofstream
#include <math.h>
#include "kmers.hpp"
#include "general.hpp"
#include "DNA.hpp"

using namespace std;

int getKmerFileInfo(const vector<string> & kmerfiles, const bool tabulated){

	vector < string > sampleNames;
	int kmersize;
	int namesPerLine=8;
	ifstream fileStream;

	if (tabulated){
		cout << "File\tKmer size\t# samples\t# kmers\t# sample patterns" << endl;
	}

	for (int s = 0; s < kmerfiles.size(); ++s){
		
		if (tabulated){
			cout << splitFileName(kmerfiles[s]);
		}
		else {
			cout << endl << splitFileName(kmerfiles[s]) << endl;
			cout << string(kmerfiles[s].size(), '=') << endl;
		}

		if (openFileStream(kmerfiles[s], fileStream, false)){return 1;};
		
		readKmerHeader(fileStream, kmersize, sampleNames);

		if (tabulated){
			cout << "\t" << kmersize << "\t" << sampleNames.size();
		}
		else {
			cout << "Kmer size: " << kmersize << endl;
			cout << "Number of samples: " << sampleNames.size() << endl;
			cout << "Sample names:"<< endl;
			int j=1;
			for (int i=0; i<sampleNames.size(); ++i, ++j){
				cout << sampleNames[i];
				if (i!=sampleNames.size()-1){
					cout << ", ";
				}
				if (j==namesPerLine || i==sampleNames.size()-1){
					cout << endl;
					j=0;
				}
			}
		}

		int kmercount=0;
		int patterncount=0;
		char kmerbuffer[(kmersize*2/3)+1];
		char asciibuffer[int(ceil(float(sampleNames.size())/6))];

		while (fileStream.read(asciibuffer, sizeof(asciibuffer))){//read the ascii representation of the taxa
			patterncount++;
			while (fileStream.peek()!='\n' && fileStream.read(kmerbuffer, sizeof(kmerbuffer))){//read the kmers on the line
				kmercount++;
			}
			fileStream.ignore(256,'\n');//skip the end ofline character
		}
		if (tabulated){
			cout << "\t" << kmercount << "\t" << patterncount << endl;
		}
		else {
			cout << "Number of kmers: " << kmercount << endl;
			cout << "Number of sample patterns: " << patterncount << endl << endl;
		}
		sampleNames.clear();
		fileStream.close();
	}
	return 0;
}