//g++ -O3 -std=c++0x src/sk_fastq -lz -o bin/sk_fastq
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

//int main(int argc, char *argv[])
int fastqToKmers(const vector<string> & fastqs, const string & outfilename, const int & kmerlen, const int & userminquality, const int & userfilecutoff, const int & usercovcutoff, const float & userminmaf)
{

	auto start = chrono::high_resolution_clock::now();
	
	int minquality=userminquality+33;

	
	// Create the kmer map
	unordered_map<string, array<int,5> > kmerMap;
	int substringlength=(kmerlen*2)+1;
	int numseqs=0;
	int numbases=0;
	int numreadbases=0;
	int lastnumbases=0;

	string sequence;
	sequence.reserve(1000);
	string quality;
	quality.reserve(1000);
	string line;
	line.reserve(1000);
	string subsequence;
	subsequence.reserve(1000);
	int i;
	string kmer;
	kmer.reserve((kmerlen*2)+1);
	char base;

	for (int s = 0; s < fastqs.size(); ++s){
		
		cout << "Reading " << fastqs[s] << endl;
		igzstream gzfileStream;
		gzfileStream.open(fastqs[s].c_str());

		int count=0;
		while (true){
			if (gzfileStream.peek()!='@'){
				cout << fastqs[s] << " is not in the correct format. Expecting header line to start with @." << endl << endl;
				return 1;
			}
		  	getline(gzfileStream, line);
		  	getline(gzfileStream, sequence);
		  	if (gzfileStream.peek()!='+'){
				cout << fastqs[s] << " is not in the correct format. Expecting separator line to start with +." << endl << endl;
				return 1;
			}
		  	getline(gzfileStream, line);
		  	getline(gzfileStream, quality);
		  	if (quality.length()!=sequence.length()){
		  		cout << fastqs[s] << " is not in the correct format. Sequence and quality lines must be of equal length." << endl << endl;
				return 1;
		  	}

	  		if (not gzfileStream.good()){
	  			cout << fastqs[s] << " is not in the correct format. Expecting a multiple of 4 lines." << endl << endl;
	  			return 1;
			}


			numseqs++;
			numreadbases+=sequence.length();

			lowqualitytoN(sequence, quality, minquality);
			
			stringstream sequencestream;
			sequencestream << sequence;//convert the sequence to stringstream

			
			while (getline(sequencestream, subsequence, 'N')){
				
				if (subsequence.length()<substringlength){
					continue;
					}
				
				transform(subsequence.begin(), subsequence.end(), subsequence.begin(), ::toupper);//change all letters in the string to upper case	
				
				i=0;
				
				for (auto iti = subsequence.cbegin(), end = subsequence.cend()-(substringlength-1); iti != end; ++iti, ++i){
					numbases+=1;
					kmer=subsequence.substr(i,substringlength);
					
					if (reverse_is_min(kmer, kmerlen)){
						reverse(kmer.begin(), kmer.end());
						transform(kmer.begin(),kmer.end(),kmer.begin(),complement);
					}
					
					base=kmer[kmerlen];
					kmer.erase(kmer.begin()+kmerlen);
					
					auto it = kmerMap.find(kmer);//check if the kmer is in the hash
					
					if ( it != kmerMap.end() ){//if the kmer is in the hash
						it->second[base_score[int(base)]]++;//increment the count of the base for the kmer
						if (s>0){//if this is the second file
							it->second[4]++;//increment the count of the kmer
						}
					}
					else if (s==0) {//if the kmers isn't in the hash
						auto ret = kmerMap.insert(make_pair(kmer, array< int, 5>()));//insert the kmer to the hash
						ret.first->second[base_score[int(base)]]=1;//increment the count of the base for the kmer
					}
				}
			}
			if (gzfileStream.peek()==EOF){
				break;
			}
    	}
		gzfileStream.close();

		
		if (s==0 && userfilecutoff>0){
			cout << "Added " << numbases << " kmers from " << numseqs << " sequences\n";
			cout << kmerMap.size() << " unique kmers in map\n";
			cout << "Filtering kmers for file coverage\n";

			auto it = kmerMap.begin();
			auto endIter = kmerMap.end();

			for (; it!=endIter; ){
				int basecoverage=0;
				for (auto it2=it->second.begin(); it2!=it->second.end()-1; ++it2){
					basecoverage+=*it2;
				}
				//filter sinlgletons here
				for (auto it2=it->second.begin(); it2!=it->second.end()-1; ++it2){
					if (*it2<userfilecutoff){
						basecoverage-=*it2;
						*it2=0;
					}
				}
				
				if (basecoverage==0){
						kmerMap.erase(it++);
					}
				else{
					++it;
				}
				
			}
			cout << kmerMap.size() << " unique kmers in map after filtering\n";
		}
		
		
	}
	cout << "Added " << numbases << " kmers from " << numseqs << " sequences of total length " << numreadbases << "\n";
	cout << "Running final filtering of kmers\n";

	int totalcoverage=0;

	auto it = kmerMap.begin();
	auto endIter = kmerMap.end();
	for (; it!=endIter; ){

		//remove kmers with file coverage lower than the user defined cutoff if there were two fastq files supplied
		if (fastqs.size()>1 && it->second[4]<userfilecutoff){
			kmerMap.erase(it++);
			continue;
		}
		else{
			it->second[4]=0;
		}

		//calculate the total kmer coverage
		int basecoverage=0;
		for (auto it2=it->second.begin(); it2!=it->second.end()-1; ++it2){
			basecoverage+=*it2;
		}

		float covcutoff=float(basecoverage)*userminmaf;
		if (covcutoff<usercovcutoff){//change 8 here to change minimum acceptable coverage
			covcutoff=usercovcutoff;
		}

		//filter bases that don't meet the coverage cutoff
		for (auto it2=it->second.begin(); it2!=it->second.end()-1; ++it2){
			if (*it2<covcutoff){
				basecoverage-=*it2;
				*it2=0;
			}
		}

		//remove kmers with coverage lower than the user defined cutoff
		if (basecoverage<usercovcutoff){
				kmerMap.erase(it++);
			}
		else{
			totalcoverage+=basecoverage;
			++it;
		}
		
	}
	cout << kmerMap.size() << " unique kmers in map after final filtering\n";
	cout << "Mean kmer coverage is " << float(totalcoverage)/kmerMap.size() << "\n";
	
	cout << "Writing kmers to " << outfilename << "\n";
	
	printkmerfile(kmerMap, outfilename, kmerlen);


	auto finish = std::chrono::high_resolution_clock::now();
	chrono::duration<double> elapsed = finish - start;
	cout << "Done\n";
	cout << "Total time required: " << elapsed.count() << "s\n";
	
	return 0;
	
}


