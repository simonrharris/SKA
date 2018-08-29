#include <unordered_map> //std::unordered_map
#include <iostream> //std::cout
#include <fstream> //std::ofstream
#include <string> //std::string
#include <vector> //std::vector
#include <set> //std::set
#include <chrono> //std::chrono
#include <cmath> //std::ceil
#include <algorithm> //std::count
#include "kmers.hpp"
#include "general.hpp"
#include "DNA.hpp"


int addKmersByFilters (std::set < std::string > & myKmerSet, const std::vector < std::string > & myFiles, const int minSamples, const int maxSamples, int & myKmerSize, const int numSamples) {
	std::unordered_map < std::string, int > kmerCountMap;
	std::ifstream myFileStream;

	if (minSamples>0){
		std::cout << "Weeding kmers that are present in fewer than " << minSamples << " samples" << std::endl;
	}
	if (maxSamples<numSamples){
		std::cout << "Weeding kmers that are present in more than " << maxSamples << " samples" << std::endl;
	}

	for (std::vector < std::string >::const_iterator it = myFiles.begin(); it != myFiles.end(); ++it){

		if (openFileStream(*it, myFileStream)){return 1;};

		int kmersize;
		std::vector < std::string > names;
		
		readKmerHeader(myFileStream, kmersize, names);
		
		if (myKmerSize==0){
			myKmerSize=kmersize;
		}
		if (kmersize!=myKmerSize){
			std::cerr << "Files have different kmer sizes" << std::endl << std::endl;
			return 1;
		}

		bool firstkmer=true;
		char base;
		char kmerbuffer[kmersize*2/3];
		char asciibuffer[int(ceil(float(names.size())/6))];
		while (myFileStream.read(asciibuffer, sizeof(asciibuffer))){
			std::string asciibits (asciibuffer, sizeof(asciibuffer));
			std::vector < bool > mybits;
			vectorbool_from_ascii(asciibits, mybits);//read the ascii representation of the taxa to a vector of bools
			int matchcount=std::count(mybits.begin(), mybits.end(), true);
			
			while (myFileStream.peek()!='\n' && myFileStream.get(base)){
				myFileStream.read(kmerbuffer, sizeof(kmerbuffer));
				std::string kmer (kmerbuffer, kmersize*2/3);
				
				std::unordered_map < std::string, int >::iterator kcmit = kmerCountMap.find(kmer);//check if the kmer is in the hash
				if ( kcmit != kmerCountMap.end() ){//if the kmer is not in the hash
					kcmit->second+=matchcount;
				}
				else {
					kmerCountMap.insert(std::make_pair(kmer, matchcount));
				}
			}
			myFileStream.ignore(256,'\n');//skip the end ofline character
    	}
		myFileStream.close();
	}

	std::unordered_map < std::string, int >::iterator kcmit = kmerCountMap.begin();
	std::unordered_map < std::string, int >::iterator kcmend = kmerCountMap.end();

	for (; kcmit!=kcmend; ){
		if (kcmit->second<minSamples || kcmit->second>maxSamples){
			//std::cout << kcmit->first << " " << kcmit->second << std::endl;
			myKmerSet.insert(kcmit->first);
		}
		kmerCountMap.erase(kcmit++);
	}
	return 0;
}


int weedKmers(const std::vector < std::string > & weedfiles, const std::string & kmerfile, const float minproportion, const float maxproportion, int minsamples, int maxsamples)
{

	const std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
	// Create the kmer map
	std::set < std::string > kmerSet;
	
	std::ifstream fileStream;

	

	int kmersize=0;
	char base;
	
	

	if (kmerfile!=""){
		if (openFileStream(kmerfile, fileStream)){return 1;};
		std::vector < std::string > weednames;
		
		readKmerHeader(fileStream, kmersize, weednames);
		char kmerbuffer[kmersize*2/3];
		char asciibuffer[int(ceil(float(weednames.size())/6))];
		
		while (fileStream.read(asciibuffer, sizeof(asciibuffer))){
			std::string asciibits (asciibuffer, sizeof(asciibuffer));
			
			while (fileStream.peek()!='\n' && fileStream.get(base)){
				fileStream.read(kmerbuffer, sizeof(kmerbuffer));
				std::string kmer (kmerbuffer, sizeof(kmerbuffer));
				kmerSet.insert(kmer);
			}
			fileStream.ignore(256,'\n');//skip the end ofline character
	    }
		fileStream.close();
		
		std::cout << kmerSet.size() << " unique kmers in weed set" << std::endl;

		
	}


	std:: vector < std::string > sampleNames;
	if (collectSampleNames(weedfiles, sampleNames)!=0){return 1;}

	if (maxsamples==0){
		maxsamples=sampleNames.size();
	}
	else if (maxsamples>sampleNames.size()){
		std::cerr << "Warning: Maxsamples is greater than the number of samples in your files." << std::endl;
	}
	if (minsamples>sampleNames.size()){
		std::cerr << "Error: Minsamples is greater than the number of samples in your files." << std::endl;
		return 1;
	}

	int myminsamples=std::ceil(minproportion*sampleNames.size());
	int mymaxsamples;
	if (maxproportion==1.0){
		mymaxsamples=sampleNames.size();
	}
	else{
		mymaxsamples=std::floor(maxproportion*sampleNames.size());
	}

	if (myminsamples>minsamples){
		minsamples=myminsamples;
	}
	if (mymaxsamples<maxsamples){
		maxsamples=mymaxsamples;
	}

	if (minsamples>0 || maxsamples<sampleNames.size()){
		if (addKmersByFilters (kmerSet, weedfiles, minsamples, maxsamples, kmersize, sampleNames.size())){return 1;};
		std::cout << kmerSet.size() << " unique kmers in weed set" << std::endl;
	}
	
	int weeded=0;
	int totalweedkmers=0;
	
	for (std::vector < std::string >::const_iterator it = weedfiles.begin(); it != weedfiles.end(); ++it){

		if (openFileStream(*it, fileStream)){return 1;};

		int weedkmersize;
		std::vector < std::string > names;
		
		readKmerHeader(fileStream, weedkmersize, names);
		char kmerbuffer[weedkmersize*2/3];
		
		if (weedkmersize!=kmersize){
			std::cerr << "Files have different kmer sizes" << std::endl << std::endl;
			return 1;
		}
		std::string filename=std::string(*it);
		std::string outputfile=filename.substr(0, filename.find_last_of("."))+".weeded"+filename.substr(filename.find_last_of("."), filename.length());

		std::cout << "Writing weeded kmers to " << outputfile << std::endl;
	
		std::ofstream weededfile(outputfile);
		weededfile << "SKA v" << versionNumber << std::endl;
		weededfile << kmersize << std::endl;
		for ( std::vector < std::string >::iterator it=names.begin(); it!=names.end(); ++it){ //print each sample name to output file stream
			weededfile << *it << " "; 
		}
		weededfile << std::endl;

		bool firstkmer=true;

		char asciibuffer[int(ceil(float(names.size())/6))];
		while (fileStream.read(asciibuffer, sizeof(asciibuffer))){
			std::string asciibits (asciibuffer, sizeof(asciibuffer));
			
			while (fileStream.peek()!='\n' && fileStream.get(base)){
				fileStream.read(kmerbuffer, sizeof(kmerbuffer));
				std::string kmer (kmerbuffer, kmersize*2/3);
				
				std::set < std::string >::iterator it2 = kmerSet.find(kmer);//check if the kmer is in the hash
				if ( it2 == kmerSet.end() ){//if the kmer is not in the hash
					if (firstkmer){//if it's the first kmer, print the ascii bitstring
						weededfile << asciibits;
						firstkmer=false;
					}
					weededfile << base << kmer;//print the kmer
					weeded++;
				}
				totalweedkmers++;
			}
			if (firstkmer==false){
				weededfile << std::endl;
			}
			firstkmer=true;
			fileStream.ignore(256,'\n');//skip the end ofline character
    	}
		fileStream.close();
		weededfile.close();
	}
	std::cout << totalweedkmers << " kmers prior to weeding" << std::endl;
	std::cout << weeded << " kmers remaining after weeding" << std::endl;
	std::cout << totalweedkmers-weeded << " kmers weeded (Note this may be more than the number of kmers in the weed set due to mulitple variants of a kmer being in the file)" << std::endl;

	printDuration(start);
	
	return 0;
	
	
}


