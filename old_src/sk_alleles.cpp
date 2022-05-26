#include <unordered_map> //std::unordered_map
#include <iostream> //std::cout
#include <string> //std::string
#include <sstream> //std::stringstream
#include <vector> //std::vector
#include <chrono> //timing
#include <set> //set
#include "kmers.hpp"
#include "general.hpp"
#include "gzstream.h"
#include "DNA.hpp"


int reverseKmerMap(std::unordered_map < std::string, std::string > & myKmerMap, std::unordered_map < std::vector < bool >, std::vector < std::string > > & myRevKmerMap, const int alleleCount){
	std::unordered_map < std::string, std::string >::iterator kmit = myKmerMap.begin();
	std::unordered_map < std::string, std::string >::iterator endKmit = myKmerMap.end();
	int kmercount = 0;
	for (; kmit!=endKmit; ){
		std::string kmer=kmit->first;
		if(ascii_codons(kmer)){return 1;}
		std::set < char > mybases;
		for (std::string::iterator it2=kmit->second.begin(); it2!=kmit->second.end(); ++it2){
			if (*it2!='-'){
				mybases.insert(*it2);
			}
		}
		for (std::set < char >::iterator it2=mybases.begin(); it2!=mybases.end(); ++it2){
			kmercount++;
			std::vector < bool > baseBitString(alleleCount, false);
			int i=0;
			for (std::string::iterator it3=kmit->second.begin(); it3!=kmit->second.end(); ++it3, ++i){
				if (*it2==*it3){
					baseBitString[i]=true;
				}
			}
			std::string kmerString = *it2 + kmer;
			std::unordered_map < std::vector < bool >,  std::vector < std::string > >::iterator it3 = myRevKmerMap.find(baseBitString);//check if the bitset is in the map
			if ( it3 != myRevKmerMap.end() ){//if the bitset is in the map
				it3->second.push_back(kmerString); //add the kmer
			}
			else {//if the bitset isn't in the map
				std::vector < std::string > myvector; //create a new vector of strings
				myvector.push_back(kmerString); //add the bitset to the vector
				myRevKmerMap.insert(std::make_pair(baseBitString, myvector)); //insert the vector into the map
			}
		}
		myKmerMap.erase(kmit++); //remove kmer from kmerMap to free up memory
	}
	return 0;
}


void addKmerToStringMap(std::unordered_map < std::string, std::string > & myKmerMap, const std::string & myKmer, const char myBase, const int alleleNumber, const int alleleCount){
	std::unordered_map < std::string, std::string >::iterator kmit = myKmerMap.find(myKmer);//check if the kmer is in the hash				
	if ( kmit != myKmerMap.end() ){//if the kmer is in the hash
		if (kmit->second[alleleNumber]=='-'){
			kmit->second[alleleNumber]=myBase;//add the base to the correct sample
		}
		else if (kmit->second[alleleNumber]!=myBase){
			kmit->second[alleleNumber]='N';//change to N if it's a repeat
		}
	}
	else {//if the kmers isn't in the hash
		std::pair < std::unordered_map < std::string, std::string >::iterator, bool > ret = myKmerMap.insert(std::make_pair(myKmer, std::string(alleleCount , '-')));//insert the kmer into the hash
		ret.first->second[alleleNumber]=myBase;//add the base to the correct sample
	}
}


int allelesToKmers(const std::vector < std::string > & alleles, const long & kmerlen)
{

	const std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();//start the clock

	int substringlength=(kmerlen*2)+1;

	for (int s = 0; s < alleles.size(); ++s){

		std::cout << "Creating split kmers for " << alleles[s] << std::endl;

		int numAlleles=0;
		std::vector < std::string > alleleNames;	

		if(countSequencesinFasta(alleles[s], numAlleles)){return 1;}//count the number of sequences in the file

		// Create the kmer map
		std::unordered_map < std::string, std::string > kmerMap;

		int alleleNumber=0;
		
		igzstream gzfileStream;
		if(openGzFileStream(alleles[s], gzfileStream)){return 1;}

		while (gzfileStream.peek()!=EOF){

			std::string alleleName;
			std::string sequence;
			sequence.reserve(10000);

			if(readNextFastaSequence(gzfileStream, alleles[s], alleleName, sequence)){return 1;}

			alleleNames.push_back(alleleName);

			if(IUPACToN(sequence)){return 1;}//convert any IUPAC charaxcters to N and reject any unrecognised characters
			
			std::stringstream sequencestream;
			sequencestream << sequence;//convert the sequence to stringstream
			std::string subsequence;
			
			while (std::getline(sequencestream, subsequence, 'N')){
				
				if (subsequence.length()<substringlength){
					continue;
				}

				for (int i=0; i<subsequence.length()-(substringlength-1); ++i){
					std::string kmer=subsequence.substr(i,substringlength);
					
					reverseComplementIfMin(kmer);
					
					char base;
					extractMiddleBase(kmer, base);
					
					addKmerToStringMap(kmerMap, kmer, base, alleleNumber, numAlleles);
				}
		
    		}
    		alleleNumber++;
		}
 		gzfileStream.close();//close the file
	
 		std::unordered_map < std::vector < bool >,  std::vector < std::string > > revKmerMap;

		if(reverseKmerMap(kmerMap, revKmerMap, numAlleles)){return 1;}

		if(printMergedKmerFile(revKmerMap, alleles[s]+".skf", alleleNames, kmerlen)){return 1;}

	}

	printDuration(start);//stop the clock
	
	return 0;

}


