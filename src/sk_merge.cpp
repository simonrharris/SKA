#include <unordered_map> //std::unordered_map
#include <iostream> //std::cout std::cerr
#include <fstream> //std::ofstream std::ifstream
#include <string> //std::string
#include <vector> //std::vector
#include <cmath>  //std::ceil
#include <chrono> //std::chrono
#include <algorithm> //std::count
#include "general.hpp"
#include "kmers.hpp"


int mergeKmerFiles(const std::string & outfile, const std::vector < std::string > & kmerfiles, const std::vector < std::string > & sample)
{

	const std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();//start the clock

	int numfiles=kmerfiles.size();
	
	// Create the kmer map
	std::unordered_map < std::string, std::vector < bool > > kmerMap;
	
	std::vector< int > kmerCounts;
	int kmersize=0;
	std::vector < std::string > sampleNames;

	if (collectSampleNames(kmerfiles, sampleNames)!=0){return 1;}

	std::vector < bool > include (sampleNames.size());
	getSubsample(sample, sampleNames, include);//get a vector of bools representing which samples to include based on a provided sample file. This also removed duplicate samples in the input files

	std::vector < std::string > includedSampleNames;
	for (std::vector < bool >::iterator it=include.begin(); it!=include.end(); ++it){//make a vector of all included sample names to help printing output files later
		if (*it){
			includedSampleNames.push_back(sampleNames[distance(include.begin(), it)]);
		}
	}

	//int numSamples=sampleNames.size();
	int numSamples = std::count(include.begin(), include.end(), true);

	std::vector < std::string > fileSampleNames;

	int sampleNum=0;
	int includedSampleNum=0;
	std::ifstream fileStream;

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

		std::vector < int > fileInclude;

		for (int i=0; i<fileSampleNames.size(); ++i){ //put the index of all sample names that are to be included into a vector
			if (include[sampleNum+i]){
				fileInclude.push_back(i);
			}
		}
		if (fileInclude.size()==0){ // if no sample names in the file are going to be included then don't read the file
			fileStream.close();
			sampleNum+=int(fileSampleNames.size());
			fileSampleNames.clear();
			continue;
		}

		//int numFileSamples=fileSampleNames.size();
		int numFileSamples = fileInclude.size();

		char basebuffer[1];
		char kmerbuffer[(kmersize*2/3)+1];
		char asciibuffer[int(std::ceil(float(fileSampleNames.size())/6))];

		while (fileStream.read(asciibuffer, sizeof(asciibuffer))){//read the ascii representation of the taxa
			std::string asciibits (asciibuffer, sizeof(asciibuffer));
			std::vector < bool > mybits;
			vectorbool_from_ascii(asciibits, mybits);//read the ascii representation of the taxa to a vector of bools

			int mybitcount=0;
			for (int i=0; i< fileInclude.size(); ++i){
				if (mybits[fileInclude[i]]){
					mybitcount++;
				}
			}

			if (mybitcount==0){
				std::string line;
				std::getline(fileStream, line);
				continue;
			}

			while (fileStream.peek()!='\n' && fileStream.read(kmerbuffer, sizeof(kmerbuffer))){
				std::string kmer (kmerbuffer, (kmersize*2/3)+1);
				
				std::unordered_map < std::string, std::vector < bool > >::iterator it = kmerMap.find(kmer);//check if the kmer is in the map
				if ( it != kmerMap.end() ){//if the kmer is in the map
					for (int i=0; i<fileInclude.size(); ++i){
						it->second[i+includedSampleNum]=mybits[fileInclude[i]];
					}
				}
				else {//if the kmer is not in the hash
					std::vector < bool > taxonBitset(numSamples, false);//create a new vector of bools for all samples
					for (int i=0; i<fileInclude.size(); ++i){ //change the bools to true based for taxa in the file containing the kmer
						taxonBitset[i+includedSampleNum]=mybits[fileInclude[i]];
					}
					kmerMap.insert(std::make_pair(kmer, taxonBitset));//add the new vector to the map
				}
	    	}
	    	fileStream.ignore(256,'\n');//skip the end of line character
    	}
    	sampleNum+=numFileSamples; //add the number of samples in the file to the count of total samples
    	includedSampleNum+=fileInclude.size();
    	fileSampleNames.clear();
		fileStream.close();
	}

	std::cout << kmerMap.size() << " kmers identified from " << includedSampleNum << " samples in " << numfiles << " files" << std::endl;
	
	std::cout << "Merging..." << std::endl;

	std::unordered_map < std::vector < bool >,  std::vector < std::string > > revKmerMap;

	reverseVectorBoolKmerMap(kmerMap, revKmerMap);//reverse the map so that we have a new map of kmers for each combination of samples

	std::cout << revKmerMap.size() << " unique taxon combinations in map" << std::endl;
	
	if(printMergedKmerFile(revKmerMap, outfile, includedSampleNames, kmersize)){return 1;}
	
	printDuration(start);//stop the clock
	
	return 0;
		
}


