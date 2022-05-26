#include <iostream> //std::cout std::cerr
#include <fstream> //std::istream
#include <string> //std::string
#include <vector> //std::vector
#include <map> //std::map
#include <cmath> //std::ceil
#include "kmers.hpp"
#include "general.hpp"

int summariseKmerFiles(const std::vector < std::string > & kmerfiles)
{
	int kmersize;
	std::vector < std::string > fileSampleNames;
	std::ifstream fileStream;
	
	std::vector < char > baseorder = { 'A', 'C', 'G', 'T', 'N', 'O'};
	std::cout << "Sample\tKmer size\tTotal kmers\tAs\tCs\tGs\tTs\tNs\tOthers\tGC Content" << std::endl;
	for (int s = 0; s < kmerfiles.size(); ++s){
		
		if (openFileStream(kmerfiles[s], fileStream, false)){return 1;};
		
		readKmerHeader(fileStream, kmersize, fileSampleNames);

		std::map < char, std::vector < int > > basecounts{ {'A', std::vector < int > (fileSampleNames.size(), 0)}, {'C', std::vector < int > (fileSampleNames.size(), 0)}, {'G', std::vector < int > (fileSampleNames.size(), 0)}, {'T', std::vector < int > (fileSampleNames.size(), 0)}, {'N', std::vector < int > (fileSampleNames.size(), 0)}, {'O', std::vector < int > (fileSampleNames.size(), 0)}};

		std::vector < int > kmers(fileSampleNames.size(), 0);
		
		char base;
		char kmerbuffer[kmersize*2/3];
		char asciibuffer[int(ceil(float(fileSampleNames.size())/6))];

		while (fileStream.read(asciibuffer, sizeof(asciibuffer))){//read the ascii representation of the taxa
			std::string asciibits (asciibuffer, sizeof(asciibuffer));
			std::vector < bool > mybits;
			vectorbool_from_ascii(asciibits, mybits);//read the ascii representation of the taxa to a vector of bools

			while (fileStream.peek()!='\n' && fileStream.get(base)){
				base=toupper(base);
				fileStream.read(kmerbuffer, sizeof(kmerbuffer));
				std::string kmer (kmerbuffer, kmersize*2/3);
				
				for (int i=0; i<fileSampleNames.size(); ++i){ //add the kmer count to all samples that are true in the bitset
					if (mybits[i]){
						kmers[i]++;
						std::map < char, std::vector < int > >::iterator it = basecounts.find(base);//check if the kmer is in the map
						if ( it != basecounts.end() ){//if the kmer is in the map
							it->second[i]++;//increment the count for the base
						}
						else {//if the kmer isn't in the map
							basecounts['O'][i]++;//increment the count for other
						}
					}
				}
			}
			fileStream.ignore(256,'\n');//skip the end ofline character
    	}
		fileStream.close();
		
		for (int i=0; i<fileSampleNames.size(); ++i){ //print the results	
			std::cout << fileSampleNames[i] << "\t" << kmersize << "\t" << kmers[i] << "\t";
			for (std::vector < char >::iterator it=baseorder.begin(); it!=baseorder.end(); ++it){
				std::cout << basecounts[*it][i] << "\t";
			}
			float gccontent=(float(basecounts['C'][i])+float(basecounts['G'][i]))/(float(basecounts['A'][i])+float(basecounts['C'][i])+float(basecounts['G'][i])+float(basecounts['T'][i]));
			std::cout << gccontent << std::endl;
		}
    	fileSampleNames.clear();
	}

	return 0;
		
}


