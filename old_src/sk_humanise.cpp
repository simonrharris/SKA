#include <string> //std::string
#include <vector> //std::vector
#include <iostream> //std::cout
#include <cmath> //std::ceil
#include <chrono> //std::chrono
#include "kmers.hpp"
#include "general.hpp"
#include "DNA.hpp"

int humaniseKmers(const std::string & kmerfile, const std::string & outputfile){
	
	const std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();//start the clock

	std::ofstream humanfile(outputfile);

	std::vector < std::string > sampleNames;
	int kmersize;
	std::ifstream fileStream;

	if (openFileStream(kmerfile, fileStream, false)){return 1;};//open the file
	
	readKmerHeader(fileStream, kmersize, sampleNames);//read the header to get the kmer size and sample names

	int numSamples=sampleNames.size();

	std::vector < int > fileInclude;
	for (int i=0; i<numSamples; ++i){ //put the index of all sample names that are to be included into a vector
		fileInclude.push_back(i);
	}

	for (int i=0; i<sampleNames.size(); ++i){//print the sample names with namesPerLine sample names per line
		humanfile << "\t" << sampleNames[i];
	}
	humanfile << std::endl;

	std::unordered_map < std::string, std::string > kmerMap;

	char base;
	char kmerbuffer[kmersize*2/3];
	char asciibuffer[int(ceil(float(sampleNames.size())/6))];

	while (fileStream.read(asciibuffer, sizeof(asciibuffer))){//read the next ascii bitstring
		std::string asciibits (asciibuffer, sizeof(asciibuffer));
		std::vector < bool > mybits;
		vectorbool_from_ascii(asciibits, mybits);//convert the ascii representation of the taxa to a vector of bools

		while (fileStream.peek()!='\n' && fileStream.get(base)){
			fileStream.read(kmerbuffer, sizeof(kmerbuffer));
			std::string kmer (kmerbuffer, kmersize*2/3);
			
			addKmerToStringMap(kmerMap, kmer, base, mybits, fileInclude, 0, numSamples, numSamples);
			
    	}
    	fileStream.ignore(256,'\n');//skip the end ofline character
    }

    for (std::unordered_map < std::string, std::string >::iterator it=kmerMap.begin(); it!=kmerMap.end(); ++it){
    	std::string asciikmer=it->first;
		std::string kmer=codons_from_ascii(asciikmer);
		std::string allbases=it->second;
		humanfile << kmer;

		for (int i=0; i<allbases.length(); ++i){
			humanfile << "\t" << allbases[i];
		}
		humanfile << std::endl;

	}

	printDuration(start);

	return 0;
}