//g++ -O3 -std=c++0x src/sk_fastq -lz -o bin/sk_fastq
#include <unordered_map> //std::unordered_map
#include <iostream> //std::cout
#include <algorithm> //std::reverse std::transform
#include <string> //std::string
#include <sstream> //std::stringstream
#include <vector> //std::vector
#include <array> //std::array
#include <chrono> //std::chrono
#include "general.hpp"
#include "kmers.hpp"
#include "gzstream.h"
#include "DNA.hpp"
//#include <string_view>//Bear in mind for future that string_view allows 'in place' substrings!


int fastqToKmers(const std::vector < std::string > & fastqs, const std::string & outfilename, const int & kmerlen, const int & userminquality, const int & userfilecutoff, const int & usercovcutoff, const float & userminmaf, const bool printAlleles)
{

	const std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();//start the clock
	
	std::unordered_map < std::string, std::array < int, 8 > > kmerMap;//create a map to store the kmers for each base
	int substringlength=(kmerlen*2)+1;
	int numseqs=0;
	int numbases=0;
	int numreadbases=0;

	std::string sequence;
	sequence.reserve(500);
	std::string quality;
	quality.reserve(500);
	std::string subsequence;
	subsequence.reserve(500);
	std::string kmer;
	kmer.reserve((kmerlen*2)+1);
	char base;
	bool isFirst=true;

	for (int s = 0; s < fastqs.size(); ++s){
		
		std::cout << "Reading " << fastqs[s] << std::endl;
		igzstream gzfileStream;
		gzfileStream.open(fastqs[s].c_str());//open the fastq file
		if (gzfileStream.fail()){
			std::cerr << std::endl << "Error: Failed to open " << fastqs[s] << std::endl << std::endl;
			return 1;
		}

		while (gzfileStream.peek()!=EOF){
			if (readNextFastqSequence(gzfileStream, fastqs[s], sequence, quality)){return 1;};//read the next fastq sequence

			numseqs++;
			numreadbases+=sequence.length();

			lowqualitytoN(sequence, quality, userminquality);
			
			std::stringstream sequencestream;
			sequencestream << sequence;//convert the sequence to stringstream
			
			while (std::getline(sequencestream, subsequence, 'N')){//for each subsequence separated by Ns
				
				if (subsequence.length()<substringlength){//if the subsequence is too small then continue
					continue;
				}
				
				int i=0;
				
				for (i=0; i<subsequence.length()-(substringlength-1); ++i){//for each base in the subsequence
					
					kmer=subsequence.substr(i,substringlength);//extract the kmer from the subsequence
					
					reverseComplementIfMin(kmer);//If the rev comp of the kmer is 'smaller' then rev comp the kmer

					if(extractMiddleBase(kmer, base)){return 1;}//extract the middle base from the kmer

					addKmerToBaseArrayMap(kmerMap, kmer, base, isFirst);//add the kmer and base to the map
					
				}
				numbases+=i;
			}
    	}
		gzfileStream.close();//close the file
		isFirst=false;

		std::cout << "Added " << numbases << " kmers from " << numseqs << " sequences" << std::endl;
		std::cout << kmerMap.size() << " unique kmers in map" << std::endl;

		applyFileKmerArrayMapFilters(kmerMap, userfilecutoff);//Filter file kmers to remove those below the maf and file cov cutoff thresholds
		
	}

	applyFinalKmerArrayMapFilters(kmerMap, usercovcutoff);//Filter the kmers to remove those below the maf and cov cutoff thresholds
	
	if(printAlleles){
		if(printKmerAlleleFrequencies(kmerMap, outfilename+"_alleles.tsv")){return 1;}//print the kmers to file
	}

	if(printKmerFile(kmerMap, outfilename+".skf", kmerlen, userminmaf)){return 1;}//print the kmers to file

	printDuration(start);//stop the clock
	
	return 0;
	
}


