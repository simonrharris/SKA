#include <unordered_map> //std::unordered_map
#include <iostream> //std::cout std::cerr
#include <fstream> //std::ifstream
#include <string> //std::string
#include <vector> //std::vecto
#include <cmath> //std::ceil 
#include <iterator> // std::distance
#include "general.hpp"
#include "kmers.hpp"
#include "DNA.hpp"
#include <chrono> //std::chrono


void filterAlignment(std::unordered_map < std::string, std::string > & myKmerMap, std::vector < int > & constantBaseVector, const int numSamples, const int minrequired, const bool variantOnly){

	std::unordered_map < std::string, std::string >::iterator kmit = myKmerMap.begin();
	std::unordered_map < std::string, std::string >::iterator kmitend = myKmerMap.end();

	for (; kmit!=kmitend; ){
		std::vector < int > baseVector (6, 0);
		for (int i=0; i<numSamples; ++i){//count the number of each base in the samples
			baseVector[base_score[kmit->second[i]]]++;
		}

		int acgtTotal=0;
		int acgtCount=0;
		int constbase=-1;

		for (int i=0; i<4; ++i){//count the total number of ACGTs and the number of different bases at the site
			if (baseVector[i]>0){
				acgtCount++;
				constbase=i;//record the site as the constant base if it is present in at least one smaple. This iwll only be used if the site is constant
			}
			acgtTotal+=baseVector[i];
		}

		if (acgtTotal<minrequired || acgtTotal==0){ //if the total number of samples matching the kmer with a base is below the minimum, exclude the site from the alignment
			myKmerMap.erase(kmit++);
		}
		else if (acgtCount==1 && variantOnly){//if the site is constant and variantonlly has been selected, record the variant base and exclude the site from the alignment
			constantBaseVector[constbase]++;
			myKmerMap.erase(kmit++);
		}
		else {
			++kmit;
		}
		
	}
}

void printAlignment(const std::string & outputfilename, const std::unordered_map < std::string, std::string > & myKmerMap,std::vector < int > & constantBaseVector, const std::vector < std::string > sampleNames, const bool variantOnly){
	
	std::cout << "Printing alignment of "<< myKmerMap.size() << " sites" << endl;
	
	std::ofstream alignfile(outputfilename);

	for (int i=0; i<sampleNames.size(); ++i){
		float nonns=0.0;
		std::string mysequence (myKmerMap.size(), '-');
		int j=0;
		for (std::unordered_map < std::string, std::string >::const_iterator kmit=myKmerMap.begin(); kmit!=myKmerMap.end(); ++kmit, ++j){
			mysequence[j] = kmit->second[i];
			if(base_score[kmit->second[i]]<4){
				nonns++;
			}
		}
		alignfile << ">" << sampleNames[i] << endl << mysequence << endl;
		if ((nonns/myKmerMap.size())<0.5){
			std::cerr << "Warning: " << sampleNames[i] << " only matches " << nonns/myKmerMap.size()*100 << "% of kmers in the alignment" << endl;
		}
	}
	alignfile.close();

	if (variantOnly){//if the variant only option was selected, print the constant site counts
		std::cout << "Constant sites matching filters (a c g t):" << endl;
		std::cout << constantBaseVector[0] << " " << constantBaseVector[1]  << " " << constantBaseVector[2]  << " " << constantBaseVector[3]  << endl;
	}
}


void addKmerToStringMap(std::unordered_map < std::string, std::string > & myKmerMap, const std::string & myKmer, const char myBase, const std::vector < bool > & myBits, const std::vector < int > & mySamples, int currentSampleNumber, int totalSamples, int maxMissing){
	std::unordered_map < std::string, std::string >::iterator kmit = myKmerMap.find(myKmer);//check if the kmer is in the map
	if ( kmit != myKmerMap.end() ){//if the kmer is in the map
		for (int i=0; i<mySamples.size(); ++i){ //add the base to all samples that are true in the bitset
			if (myBits[mySamples[i]]){
				kmit->second[i+currentSampleNumber]=myBase;
			}
		}
	}
	else {//if the kmer isn't in the map
		if ((currentSampleNumber)<=maxMissing){
			std::string newsequence (totalSamples , '-');//create an empty sequence
			for (int i=0; i<mySamples.size(); ++i){ //add the base to all included samples that are true in the bitset
				if (myBits[mySamples[i]]){
					newsequence[i+currentSampleNumber]=myBase;
				}
			}
			myKmerMap.insert(std::make_pair(myKmer, newsequence));
		}
	}
}


int alignKmers(const float & minproportion, const std::string & outputfile, const std::vector < std::string > & kmerfiles, const bool & variantonly, const std::vector < std::string > & sample)
{

	const std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();//start the clock

	int numfiles=kmerfiles.size();//record the number of files

	std::vector < std::string > sampleNames;
	if (collectSampleNames(kmerfiles, sampleNames)!=0){return 1;}//get the total number of samples in all input files

	std::vector < bool > include (sampleNames.size());
	getSubsample(sample, sampleNames, include);//get a vector of bools representing which samples to include based on a provided sample file. This also removed duplicate samples in the input files

	std::vector < std::string > includedSampleNames;
	for (std::vector < bool >::iterator vbit=include.begin(); vbit!=include.end(); ++vbit){//make a vector of all included sample names to help printing output files later
		if (*vbit){
			includedSampleNames.push_back(sampleNames[std::distance(include.begin(), vbit)]);
		}
	}

	int numSamples = std::count(include.begin(), include.end(), true);//count the number of included samples

	std::cout << numSamples << " samples will be included in the alignment" << endl;
	
	float maxmissing=round((1.0-minproportion)*numSamples);//calculate the maximum number of samples that can be missing for the site to be included in the alignment

	float minrequired=numSamples-maxmissing;//calculate the minimum number of samples that must contain a kmer for the site to be included in the alignment
	
	std::cout << "Keeping variants for which at least " << int(minrequired) << " samples include kmer matches" << endl << endl;
	
	// Create the kmer map
	std::unordered_map < std::string, std::string > kmerMap;
	int oldkmersize=0;
	int sampleNum=0;
	int includedSampleNum=0;

	char base;

	std::ifstream fileStream;
	for (int s = 0; s < kmerfiles.size(); ++s){//read each file and make a map of the kmers which stores the bases for each included sample
		
		if (openFileStream(kmerfiles[s], fileStream)){return 1;};

		int kmersize;
		std::vector < std::string > names;
		
		readKmerHeader(fileStream, kmersize, names);//read the header from the kmer file to get the kmer size and sample names
		
		if (s==0){ //If it's the first file, set the oldkmersize
			oldkmersize=kmersize;
		}

		if (kmersize!=oldkmersize){ //if the file kmer size isn't the same as the other files then ditch out
			cerr << "kmer files have different kmer sizes" << endl << endl;
			return 1;
		}

		std::vector < int > fileInclude;

		for (int i=0; i<names.size(); ++i){ //put the index of all sample names that are to be included into a vector
			if (include[sampleNum+i]){
				fileInclude.push_back(i);
			}
		}

		int fileIncludeSize=fileInclude.size();

		if (fileIncludeSize==0){ // if no sample names in the file are going to be included then don't read the file
			fileStream.close();
			sampleNum+=int(names.size());
			names.clear();
			continue;
		}

		char kmerbuffer[kmersize*2/3];
		char asciibuffer[int(ceil(float(names.size())/6))];

		while (fileStream.read(asciibuffer, sizeof(asciibuffer))){//read the next ascii bitstring
			string asciibits (asciibuffer, sizeof(asciibuffer));
			std::vector < bool > mybits;
			vectorbool_from_ascii(asciibits, mybits);//convert the ascii representation of the taxa to a vector of bools

			while (fileStream.peek()!='\n' && fileStream.get(base)){
				base=toupper(base);
				fileStream.read(kmerbuffer, sizeof(kmerbuffer));
				std::string kmer (kmerbuffer, kmersize*2/3);
				
				addKmerToStringMap(kmerMap, kmer, base, mybits, fileInclude, includedSampleNum, numSamples, maxmissing);
				
	    	}
	    	fileStream.ignore(256,'\n');//skip the end ofline character
	    }
	    sampleNum+=int(names.size()); //add the number of samples in the file to the count of total samples
	    includedSampleNum+=fileIncludeSize;
		fileStream.close();
	}
	
	std::cout << kmerMap.size() << " kmers identified from " << numSamples << " samples in " << numfiles << " files" << endl;
	
	std::vector < int > constantBases (4,0);

	filterAlignment(kmerMap, constantBases, numSamples, minrequired, variantonly);
	
	printAlignment(outputfile, kmerMap, constantBases, includedSampleNames, variantonly);

	printDuration(start);
	
	return 0;
	
}


