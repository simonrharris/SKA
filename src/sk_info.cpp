#include <string> //std::string
#include <vector> //std::vector
#include <iostream> //std::cout
#include <cmath> //std::ceil
#include "kmers.hpp"
#include "general.hpp"

int getKmerFileInfo(const std::vector< std::string > & kmerfiles, const bool tabulated){

	std::vector < std::string > sampleNames;
	int kmersize;
	int namesPerLine=8;
	std::ifstream fileStream;

	if (tabulated){
		std::cout << "File\tKmer size\t# samples\t# kmers\t# sample patterns" << std::endl;
	}

	for (int s = 0; s < kmerfiles.size(); ++s){//for each input kmer or kmerge file
		
		if (tabulated){//print the file name
			std::cout << splitFileName(kmerfiles[s]);
		}
		else {
			std::cout << std::endl << splitFileName(kmerfiles[s]) << std::endl;
			std::cout << std::string(kmerfiles[s].size(), '=') << std::endl;
		}

		if (openFileStream(kmerfiles[s], fileStream, false)){return 1;};//open the file
		
		readKmerHeader(fileStream, kmersize, sampleNames);//read the header to get the kmer size and sample names

		if (tabulated){//print the kmer size and sample information
			std::cout << "\t" << kmersize << "\t" << sampleNames.size();
		}
		else {
			std::cout << "Kmer size: " << kmersize << std::endl;
			std::cout << "Number of samples: " << sampleNames.size() << std::endl;
			std::cout << "Sample names:"<< std::endl;
			int j=1;
			for (int i=0; i<sampleNames.size(); ++i, ++j){//print the sample names with namesPerLine sample names per line
				std::cout << sampleNames[i];
				if (i!=sampleNames.size()-1){
					std::cout << ", ";
				}
				if (j==namesPerLine || i==sampleNames.size()-1){
					std::cout << std::endl;
					j=0;
				}
			}
		}

		int kmercount=0;
		int patterncount=0;
		char kmerbuffer[(kmersize*2/3)+1];
		char asciibuffer[int(std::ceil(float(sampleNames.size())/6))];

		while (fileStream.read(asciibuffer, sizeof(asciibuffer))){//read the ascii representation of the samples
			patterncount++;//increment the number of patterns
			while (fileStream.peek()!='\n' && fileStream.read(kmerbuffer, sizeof(kmerbuffer))){//read the kmers on the line
				kmercount++;//increment the number of kmers
			}
			fileStream.ignore(256,'\n');//skip the end ofline character
		}
		if (tabulated){//print the number of kmers and patterns
			std::cout << "\t" << kmercount << "\t" << patterncount << std::endl;
		}
		else {
			std::cout << "Number of kmers: " << kmercount << std::endl;
			std::cout << "Number of sample patterns: " << patterncount << std::endl << std::endl;
		}
		sampleNames.clear();//clear the sample names vector
		fileStream.close();//close the file
	}
	return 0;
}