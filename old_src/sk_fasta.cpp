#include <unordered_map> //std::unordered_map
#include <iostream> //std::cout std::cerr
#include <algorithm> //std::reverse std::transform
#include <string> //std::string
#include <sstream> //std::stringstream
#include <vector> //std::vector
#include <array> //std::array
#include <chrono> //timing
#include "general.hpp"
#include "kmers.hpp"
#include "gzstream.h"
#include "DNA.hpp"
//#include <string_view>//Bear in mind for future that string_view allows 'in place' substrings!

int fastaToKmers(const std::vector < std::string > & fastas, const std::string & outfilename, const long & kmerlen, const bool circularContigs)
{

	const std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();//start the clock

	std::unordered_map < std::string, std::array < int, 8 > > kmerMap;//create a map to store the kmers for each base
	int substringlength=(kmerlen*2)+1;
	int numseqs=0;
	int numbases=0;
	int numreadbases=0;
	char base;

	for (std::vector < std::string >::const_iterator fastait=fastas.begin(); fastait!=fastas.end(); ++fastait){//for each fasta file
		
		std::string filename = *fastait;

		std::cout << "Reading " << filename << std::endl;
		igzstream gzfileStream;

		if(openGzFileStream(filename, gzfileStream)){return 1;}//open the fasta file. May be gzipped.

		while (gzfileStream.peek() != EOF){

			std::string sequence;
			std::string name;

			if(readNextFastaSequence(gzfileStream, filename, name, sequence)){return 1;}//read the next sequence from the file

			if(IUPACToN(sequence)){return 1;}//convert any IUPAC characters to N and reject any unrecognised characters

			if (circularContigs){
				if(circulariseSequence(sequence, kmerlen)){return 1;}//if the contigs are circular, add an extra kmer lenght of sequence from the opposite end to each end
			}

			numseqs++;
			numreadbases += sequence.length();
			
			std::stringstream sequencestream;
			sequencestream << sequence;//convert the sequence to stringstream
			std::string subsequence;
			
			while (std::getline(sequencestream, subsequence, 'N')){//read each section of sequence between Ns
				
				if (subsequence.length() < substringlength){//if the subsequence is too short continue
					continue;
				}
				
				int i;

				for (i=0; i < subsequence.length() - (substringlength-1); ++i){//for each base in the subsequence

					std::string kmer = subsequence.substr(i,substringlength);//extract the kmer
					
					reverseComplementIfMin(kmer);//if the revcomp of the kmer is 'smaller' then revcomp it

					if(extractMiddleBase(kmer, base)){return 1;}//extract the middle base
					
					addKmerToBaseArrayMap(kmerMap, kmer, base, true);//add the kmer and base to the map

				}
				numbases += i;
    		}
		}
 		gzfileStream.close();//close the file
	}
	
	std::cout << "Added " << numbases << " kmers from " << numseqs << " sequences of total length " << numreadbases << std::endl;
	std::cout << kmerMap.size() << " unique kmers in map" << std::endl;
	
	printKmerFile(kmerMap, outfilename, kmerlen);//print the kmer file

	printDuration(start);//stop the clock
	
	return 0;
	
}


