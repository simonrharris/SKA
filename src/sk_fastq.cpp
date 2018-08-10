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
using namespace std;


int fastqToKmers(const vector<string> & fastqs, const string & outfilename, const int & kmerlen, const int & userminquality, const int & userfilecutoff, const int & usercovcutoff, const float & userminmaf)
{

	const chrono::high_resolution_clock::time_point start = chrono::high_resolution_clock::now();//start the clock
	
	unordered_map <string, array < int, 8 > > kmerMap;//create a map to store the kmers for each base
	int substringlength=(kmerlen*2)+1;
	int numseqs=0;
	int numbases=0;
	int numreadbases=0;

	string sequence;
	sequence.reserve(500);
	string quality;
	quality.reserve(500);
	string subsequence;
	subsequence.reserve(500);
	string kmer;
	kmer.reserve((kmerlen*2)+1);
	char base;
	bool isFirst=true;

	for (int s = 0; s < fastqs.size(); ++s){
		
		cout << "Reading " << fastqs[s] << endl;
		igzstream gzfileStream;
		gzfileStream.open(fastqs[s].c_str());//open the fastq file
		if (gzfileStream.fail()){
			cerr << endl << "Error: Failed to open " << fastqs[s] << endl << endl;
			return 1;
		}
		

		while (true){
			if (readNextFastqSequence(gzfileStream, fastqs[s], sequence, quality)){return 1;};//read the next fastq sequence

			numseqs++;
			numreadbases+=sequence.length();

			lowqualitytoN(sequence, quality, userminquality);
			
			stringstream sequencestream;
			sequencestream << sequence;//convert the sequence to stringstream
			
			while (getline(sequencestream, subsequence, 'N')){//for each subsequence separated by Ns
				
				if (subsequence.length()<substringlength){//if the subsequence is too small then continue
					continue;
				}
				
				transform(subsequence.begin(), subsequence.end(), subsequence.begin(), ::toupper);//change all letters in the string to upper case	
				
				int i=0;
				
				for (i=0; i<subsequence.length()-(substringlength-1); ++i){//for each base in the subsequence
					
					kmer=subsequence.substr(i,substringlength);//extract the kmer from the subsequence
					
					reverseComplementIfMin(kmer);//If the rev comp of the kmer is 'smaller' then rev comp the kmer

					if(extractMiddleBase(kmer, base)){return 1;}//extract the middle base from the kmer

					addKmerToBaseArrayMap(kmerMap, kmer, base, isFirst);//add the kmer and base to the map
					
				}
				numbases+=i;
			}
			if (gzfileStream.peek()==EOF){
				break;
			}
    	}
		gzfileStream.close();//close the file
		isFirst=false;

		cout << "Added " << numbases << " kmers from " << numseqs << " sequences" << endl;
		cout << kmerMap.size() << " unique kmers in map" << endl;

		int ret = applyFileKmerArrayMapFilters(kmerMap, userfilecutoff, userminmaf);//Filter file kmers to remove those below the maf and file cov cutoff thresholds
		
	}

	applyFinalKmerArrayMapFilters(kmerMap, usercovcutoff, userminmaf);//Filter the kmers to remove those below the maf and cov cutoff thresholds
	
	if(printKmerFile(kmerMap, outfilename, kmerlen)){return 1;}//print the kmers to file

	printDuration(start);//stop the clock
	
	return 0;
	
}


