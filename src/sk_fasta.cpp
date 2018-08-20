#include <unordered_map> //std::unordered_map
#include <iostream> //std::cout
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
using namespace std;

int fastaToKmers(const vector <string> & fastas, const string & outfilename, const long & kmerlen)
{

	const chrono::high_resolution_clock::time_point start = chrono::high_resolution_clock::now();//start the clock

	unordered_map < string, array < int, 8 > > kmerMap;//create a map to store the kmers for each base
	int substringlength=(kmerlen*2)+1;
	int numseqs=0;
	int numbases=0;
	int numreadbases=0;

	string sequence;
	string header;
	char base;

	for (vector <string>::const_iterator fastait=fastas.begin(); fastait!=fastas.end(); ++fastait){//for each fasta file
		
		string filename=*fastait;

		cout << "Reading " << filename << endl;
		igzstream gzfileStream;
		gzfileStream.open(filename.c_str());//open the file as a gzFileStream
		if (gzfileStream.fail()){
			cerr << endl << "Error: Failed to open " << filename << endl << endl;
			return 1;
		}

		if (gzfileStream.get()!='>'){//read the first character of the fasta file to check it is a >
			cout << filename << " is not in the correct format. Expecting header to start with >." << endl << endl;
			return 1;
		}

		while (getline(gzfileStream, header)){//read the next header line
			getline(gzfileStream, sequence, '>');//read until the next > i.e. read the sequence

			sequence.erase(std::remove_if( sequence.begin(), sequence.end(), ::isspace ), sequence.end() );//remove whitespace from the sequence

			if(IUPACToN(sequence)){return 1;}//convert any IUPAC characters to N and reject any unrecognised characters

			numseqs++;
			numreadbases+=sequence.length();
			
			stringstream sequencestream;
			sequencestream << sequence;//convert the sequence to stringstream
			string subsequence;
			
			while (getline(sequencestream, subsequence, 'N')){//read each section of sequence between Ns
				
				if (subsequence.length()<substringlength){//if the subsequence is too short continue
					continue;
				}
				
				transform(subsequence.begin(), subsequence.end(), subsequence.begin(), ::toupper);//change all 	letters in the string to upper case. Could put this into IUPACToN
				
				int i;

				for (i=0; i<subsequence.length()-(substringlength-1); ++i){//for each base in the subsequence

					string kmer=subsequence.substr(i,substringlength);//extract the kmer
					
					reverseComplementIfMin(kmer);//if the revcomp of the kmer is 'smaller' then revcomp it

					if(extractMiddleBase(kmer, base)){return 1;}//extract the middle base
					
					addKmerToBaseArrayMap(kmerMap, kmer, base, true);//add the kmer and base to the map

				}
				numbases+=i;
    		}
		}
 		gzfileStream.close();//close the file
	}
	
	cout << "Added " << numbases << " kmers from " << numseqs << " sequences of total length " << numreadbases << endl;
	cout << kmerMap.size() << " unique kmers in map" << endl;
	
	printKmerFile(kmerMap, outfilename, kmerlen);//print the kmer file

	printDuration(start);//stop the clock
	
	return 0;
	
}


