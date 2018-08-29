#include <unordered_map> //std::unordered_map
#include <iostream> //std::cout std::cerr
#include <fstream> //std::ofstream std::ifstream
#include <string> //std::string
#include <vector> //std::vector
#include <set> //std::set
#include <chrono> //std::chrono
#include <cmath>  //std::ceil
#include <algorithm> //std::count
#include "kmers.hpp"
#include "general.hpp"

int uniqueKmers(const std::vector < std::string > & ingroupsamples, const std::vector < std::string > & kmerfiles, const float & minproportion, const std::string & outputfile, const bool incNs)
{

	const std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
	// Create the kmer map
	std::unordered_map < std::string , std::vector < bool > > kmerMap;
	int oldkmersize=0;
	int totalingroupkmers=0;

	std::vector < std::string > sampleNames;

	if (collectSampleNames(kmerfiles, sampleNames)!=0){return 1;}
	int numSamples=sampleNames.size();

	float maxmissing=(1.0-minproportion)*ingroupsamples.size();

	float minrequired=ingroupsamples.size()-maxmissing;

	std::set < std::string > ingroupSet(ingroupsamples.begin(), ingroupsamples.end());

	if (ingroupSet.size()!=ingroupsamples.size()){
		std::cout << "Warning: Duplicate samples in your ingroup sample file have been removed" << std::endl;
	}

	std::vector < int > ingroupSamplePositions;

	for (int i=0; i<sampleNames.size(); ++i){
		std::set < std::string >::iterator it = ingroupSet.find(sampleNames[i]);//check if the kmer is in the map
			if ( it != ingroupSet.end() ){//if the kmer is in the map
				ingroupSamplePositions.push_back(i);
			}
	}
	
	if (ingroupsamples.size()!=ingroupSamplePositions.size()){
		std::cout << ingroupsamples.size() << " " << ingroupSamplePositions.size() << std::endl;
		std::cout << "Error: Some of your ingroup samples are not in your input files" << std::endl << std::endl;
		return 1;
	}

	std::vector < std::string > fileSampleNames;	
	int sampleNum=0;
	int numfiles=kmerfiles.size();
	std::ifstream fileStream;
	int kmersize;
	char * kmerbuffer;


	for (int s = 0; s < kmerfiles.size(); ++s){

		if (openFileStream(kmerfiles[s], fileStream)){return 1;};

		int newkmersize;
		
		readKmerHeader(fileStream, newkmersize, fileSampleNames);
		
		if (s==0){
			kmersize=newkmersize;
		}

		if (newkmersize!=kmersize){
			std::cout << "kmer files have different kmer sizes" << std::endl << std::endl;
			return 0;
		}


		int numFileSamples=fileSampleNames.size();

		char basebuffer[1];
		char kmerbuffer[(kmersize*2/3)+1];
		char asciibuffer[int(std::ceil(float(fileSampleNames.size())/6))];

		while (fileStream.read(asciibuffer, sizeof(asciibuffer))){//read the ascii representation of the taxa
			std::string asciibits (asciibuffer, sizeof(asciibuffer));
			std::vector < bool > mybits;
			vectorbool_from_ascii(asciibits, mybits);//read the ascii representation of the taxa to a vector of bools

			while (fileStream.peek()!='\n' && fileStream.read(kmerbuffer, sizeof(kmerbuffer))){
				std::string kmer (kmerbuffer, (kmersize*2/3)+1);
				
				std::unordered_map < std::string , std::vector < bool > >::iterator it = kmerMap.find(kmer);//check if the kmer is in the map
				if ( it != kmerMap.end() ){//if the kmer is in the map
					for (int i=0; i<mybits.size(); ++i){
						it->second[i+sampleNum]=mybits[i];
					}
				}
				else {//if the kmer is not in the hash
					std::vector < bool > taxonBitset(numSamples, false);//create a new vector of bools for all samples
					for (int i=0; i<mybits.size(); ++i){ //change the bools to true based for taxa in the file containing the kmer
						taxonBitset[i+sampleNum]=mybits[i];
					}
					kmerMap.insert(std::make_pair(kmer, taxonBitset));//add the new vector to the map
				}
	    	}
	    	fileStream.ignore(256,'\n');//skip the end ofline character
    	}
    	sampleNum+=numFileSamples; //add the number of samples in the file to the count of total samples
    	fileSampleNames.clear();
		fileStream.close();
	}

	std::cout << kmerMap.size() << " kmers identified from " << sampleNum << " samples in " << numfiles << " files" << std::endl;

	int uniquecount=0;
	std::unordered_map < std::string , std::vector < bool > >::iterator it = kmerMap.begin();
	std::unordered_map < std::string , std::vector < bool > >::iterator endIter = kmerMap.end();
	std::unordered_map < std::vector < bool >,  std::vector < std::string > > revKmerMap;

	for (; it!=endIter; ){

		if (not incNs && it->first[0]=='N'){
			kmerMap.erase(it++);
			continue;
		}

		int allcount=std::count(it->second.begin(), it->second.end(), true);
		int ingroupcount=0;

		for (std::vector < int >::iterator it2=ingroupSamplePositions.begin(); it2!=ingroupSamplePositions.end(); ++it2){
			if (it->second[*it2]){
				ingroupcount++;
			}
		}
		if (ingroupcount>=minrequired && (allcount-ingroupcount)==0){
			uniquecount++;

			std::vector < bool > ingroupBitString (ingroupSamplePositions.size(), false);
			int i=0;
			for (std::vector < int >::iterator it2=ingroupSamplePositions.begin(); it2!=ingroupSamplePositions.end(); ++it2, ++i){
				ingroupBitString[i]=it->second[*it2];
			}

			std::unordered_map < std::vector < bool >,  std::vector < std::string > >::iterator it2 = revKmerMap.find(ingroupBitString);//check if the bitset is in the map
			if ( it2 != revKmerMap.end() ){//if the bitset is in the map
				it2->second.push_back(it->first); //add the kmer
			}
			else {//if the bitset isn't in the map
				//std::vector < std::string > myvector; //create a new vector of strings
				//myvector.push_back(it->first); //add the bitset to the vector
				revKmerMap.insert(std::make_pair(ingroupBitString, std::vector < std::string > {it->first})); //insert the vector into the map
			}
		}
		kmerMap.erase(it++); //remove kmer from kmerMap to save space
	}

	std::cout << uniquecount << " unique shared kmers will be written to " << outputfile << std::endl;

	std::vector < std::string > orderedingroupsamples;
	for ( std::vector < int >::iterator it=ingroupSamplePositions.begin(); it!=ingroupSamplePositions.end(); ++it){ //print each sample name to output file stream
		orderedingroupsamples.push_back(sampleNames[*it]); 
	}

	if(printMergedKmerFile(revKmerMap, outputfile, orderedingroupsamples, kmersize)){return 1;}

	printDuration(start);
	
	return 0;
	
	
}


