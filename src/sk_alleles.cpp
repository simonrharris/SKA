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
#include <set> //set
#include "kmers.hpp"
#include "general.hpp"
#include "gzstream.h"
#include "DNA.hpp"
//#include <string_view>//Bear in mind for future that string_view allows 'in place' substrings!
using namespace std;

int allelesToKmers(const vector<string> & alleles, const long & kmerlen)
{

	const chrono::high_resolution_clock::time_point start = chrono::high_resolution_clock::now();//start the clock

	int substringlength=(kmerlen*2)+1;
	int numseqs=0;
	int numbases=0;
	int lastnumbases=0;

	string sequence;
	string header;
	char base;

	for (int s = 0; s < alleles.size(); ++s){

		cout << "Creating split kmers for " << alleles[s] << ".";

		int numAlleles=0;
		vector <string> alleleNames;

		//for (int s = 0; s < alleles.size(); ++s){
		igzstream gzfileStream;
		gzfileStream.open(alleles[s].c_str());

		if (gzfileStream.get()!='>'){
			cout << endl << "Error: " << alleles[s] << " is not in the correct format. Expecting header to start with >." << endl << endl;
			return 1;
		}

		while (getline(gzfileStream, header)){
			getline(gzfileStream, sequence, '>');
			numAlleles++;
		}
 		gzfileStream.close();//annoyingly have to close and reopen the file as seek/rewind aren't implemented for gzstream

		// Create the kmer map
		unordered_map < string, string > kmerMap;
		string emptySequence (numAlleles , '-');

		int alleleNumber=0;

		gzfileStream.clear();
		gzfileStream.open(alleles[s].c_str());


		//cout << gzfileStream.get();

		if (gzfileStream.get()!='>'){
			cout << endl << "Error: " << alleles[s] << " is not in the correct format. Expecting header to start with >." << endl << endl;
			continue;
		}

		while (getline(gzfileStream, header)){

			stringstream headerstream;
			headerstream << header;//convert the sequence to stringstream
			string alleleName;
			getline(headerstream, alleleName, ' ');
			alleleNames.push_back(alleleName);

			getline(gzfileStream, sequence, '>');
			sequence.erase(std::remove_if( sequence.begin(), sequence.end(), ::isspace ), sequence.end() );

			if(IUPACToN(sequence)){return 1;}//convert any IUPAC charaxcters to N and reject any unrecognised characters

			numseqs++;
			
			stringstream sequencestream;
			sequencestream << sequence;//convert the sequence to stringstream
			string subsequence;
			
			while (getline(sequencestream, subsequence, 'N')){
				
				if (subsequence.length()<substringlength){
					continue;
				}
				
				transform(subsequence.begin(), subsequence.end(), subsequence.begin(), ::toupper);//change all 	letters in the string to upper case

				for (int i=0; i<subsequence.length()-(substringlength-1); ++i){
					numbases+=1;
					string kmer=subsequence.substr(i,substringlength);
					
					reverseComplementIfMin(kmer);
					
					extractMiddleBase(kmer, base);
					
					auto it = kmerMap.find(kmer);//check if the kmer is in the hash
					
					if ( it != kmerMap.end() ){//if the kmer is in the hash
						if (it->second[alleleNumber]=='-'){
							it->second[alleleNumber]=base;//add the base to the correct sample
						}
						else if (it->second[alleleNumber]!=base){
							it->second[alleleNumber]='N';//change to N if it's a repeat
						}
					
					}
					else {//if the kmers isn't in the hash
						auto ret = kmerMap.insert(make_pair(kmer, emptySequence));//insert the kmer into the hash
						ret.first->second[alleleNumber]=base;//add the base to the correct sample
					}
				}
		
    		}
    		alleleNumber++;
		}
 		gzfileStream.close();//close the file
	

		unordered_map < vector < bool >,  vector < string > > revKmerMap;

		auto it = kmerMap.begin();
		auto endIter = kmerMap.end();
		int kmercount = 0;

		for (; it!=endIter; ){

			string kmer=it->first;
			ascii_codons(kmer);

			set < char > bases;
			for (auto it2=it->second.begin(); it2!=it->second.end(); ++it2){
				if (*it2!='-'){
					bases.insert(*it2);
				}
			}

			for (auto it2=bases.begin(); it2!=bases.end(); ++it2){

				kmercount++;

				vector < bool > baseBitString(numAlleles, false);
				int i=0;
				for (auto it3=it->second.begin(); it3!=it->second.end(); ++it3, ++i){
					if (*it2==*it3){
						baseBitString[i]=true;
					}
				}

				string kmerString = *it2 + kmer;

				auto it3 = revKmerMap.find(baseBitString);//check if the bitset is in the map
				if ( it3 != revKmerMap.end() ){//if the bitset is in the map
					it3->second.push_back(kmerString); //add the kmer
				}
				else {//if the bitset isn't in the map
					vector <string> myvector; //create a new vector of strings
					myvector.push_back(kmerString); //add the bitset to the vector
					revKmerMap.insert(make_pair(baseBitString, myvector)); //insert the vector into the map
				}
			}
			kmerMap.erase(it++); //remove kmer from kmerMap to save space
		}

		if(printMergedKmerFile(revKmerMap, alleles[s], alleleNames, kmerlen)){return 1;}

	}

	printDuration(start);//stop the clock
	
	return 0;

}


