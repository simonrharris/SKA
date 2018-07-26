#include <unordered_map> //std::unordered_map
#include <iostream> //std::cout
#include <fstream> //std::ofstream
#include <algorithm> //std::reverse std::transform
#include <cassert> //std::assert
#include <string> //std::string
#include <sstream> //std::stringstream
#include <vector> //std::vector
#include <array> //std::array
#include <chrono> //timing
#include "kmers.hpp"
#include "gzstream.h"
//#include <string_view>//Bear in mind for future that string_view allows 'in place' substrings!
using namespace std;

int fastaToKmers(const vector<string> & fastas, const string & outfilename, const long & kmerlen)
{


	auto start = chrono::high_resolution_clock::now();

	// Create the kmer map
	unordered_map<string, array<int, 5> > kmerMap;

	int substringlength=(kmerlen*2)+1;
	int numseqs=0;
	int numbases=0;
	int lastnumbases=0;

	string sequence;
	string header;

	for (int s = 0; s < fastas.size(); ++s){


		
		cout << "Reading " << fastas[s] << endl;
		igzstream gzfileStream;
		gzfileStream.open(fastas[s].c_str());

		if (gzfileStream.get()!='>'){
			cout << fastas[s] << " is not in the correct format. Expecting header to start with >." << endl << endl;
		}

		while (getline(gzfileStream, header)){
			getline(gzfileStream, sequence, '>');
			sequence.erase(std::remove_if( sequence.begin(), sequence.end(), ::isspace ), sequence.end() );

			numseqs++;
			
			stringstream sequencestream;
			sequencestream << sequence;//convert the sequence to stringstream
			string subsequence;
			
			
			while (getline(sequencestream, subsequence, 'N')){
				
				if (subsequence.length()<substringlength){
					continue;
				}
				
				transform(subsequence.begin(), subsequence.end(), subsequence.begin(), ::toupper);//change all 	letters in the string to upper case
				
				int i=0;


				for (auto iti = subsequence.cbegin(), end = subsequence.cend()-(substringlength-1); iti != end; ++iti, ++i){
					numbases+=1;
					string kmer=subsequence.substr(i,substringlength);
					
					if (reverse_is_min(kmer, kmerlen+1)){
						reverse(kmer.begin(), kmer.end());
						transform(kmer.begin(),kmer.end(),kmer.begin(),complement);
					}
					
					char base=kmer[kmerlen];
					kmer.erase(kmer.begin()+kmerlen);
					
					auto it = kmerMap.find(kmer);//check if the kmer is in the hash
					
					if ( it != kmerMap.end() ){//if the kmer is in the hash
						it->second[base_score[int(base)]]++;//increment the count of the kmer
					
					}
					else {//if the kmers isn't in the hash
						auto ret = kmerMap.insert(make_pair(kmer, array< int, 5>()));//insert the kmer into the hash
						ret.first->second[base_score[int(base)]]=1;//increment the count of the kmer
					}
				}
		
    		}
		}
 		gzfileStream.close();
	}
	
	cout << kmerMap.size() << " unique kmers in map" << endl;
	
	cout << "Writing kmers to " << outfilename << endl;
	
	printkmerfile(kmerMap, outfilename, kmerlen);

	auto finish = std::chrono::high_resolution_clock::now();
	chrono::duration<double> elapsed = finish - start;
	cout << "Done" << endl;
	cout << "Total time required: " << elapsed.count() << "s" << endl << endl;
	
	return 0;
	
	
}


