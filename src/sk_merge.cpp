//g++ -O3 -std=c++0x splitkmer.cpp -lz -o splitkmer
#include <unordered_map>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <set>
#include <vector>
#include <math.h>
#include <algorithm>
#include "general.hpp"
#include "kmers.hpp"
#include "DNA.hpp"
#include <chrono> //timing
using namespace std;


//int main(int argc, char *argv[])
int mergeKmerFiles(const string & outfile, const vector<string> & kmerfiles, const vector <string> & sample)
{

	const chrono::high_resolution_clock::time_point start = chrono::high_resolution_clock::now();//start the clock

	int numfiles=kmerfiles.size();
	
	// Create the kmer map
	unordered_map < string, vector < bool > > kmerMap;
	
	vector< int > kmerCounts;
	int kmersize=0;
	vector < string > sampleNames;

	if (collectSampleNames(kmerfiles, sampleNames)!=0){return 1;}

	vector < bool > include (sampleNames.size());
	getSubsample(sample, sampleNames, include);//get a vector of bools representing which samples to include based on a provided sample file. This also removed duplicate samples in the input files

	vector < string > includedSampleNames;
	for (auto it=include.begin(); it!=include.end(); ++it){//make a vector of all included sample names to help printing output files later
		if (*it){
			includedSampleNames.push_back(sampleNames[distance(include.begin(), it)]);
		}
	}

	//int numSamples=sampleNames.size();
	int numSamples = count(include.begin(), include.end(), true);

	vector < string > fileSampleNames;

	int sampleNum=0;
	int includedSampleNum=0;
	ifstream fileStream;

	for (int s = 0; s < kmerfiles.size(); ++s){

		if (openFileStream(kmerfiles[s], fileStream)){return 1;};

		int newkmersize;
		
		readKmerHeader(fileStream, newkmersize, fileSampleNames);
		
		if (s==0){
			kmersize=newkmersize;
		}

		if (newkmersize!=kmersize){
			cout << "kmer files have different kmer sizes" << endl << endl;
			return 0;
		}

		vector < int > fileInclude;

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
		char asciibuffer[int(ceil(float(fileSampleNames.size())/6))];

		while (fileStream.read(asciibuffer, sizeof(asciibuffer))){//read the ascii representation of the taxa
			string asciibits (asciibuffer, sizeof(asciibuffer));
			vector < bool > mybits;
			vectorbool_from_ascii(asciibits, mybits);//read the ascii representation of the taxa to a vector of bools

			int mybitcount=0;
			for (int i=0; i< fileInclude.size(); ++i){
				if (mybits[fileInclude[i]]){
					mybitcount++;
				}
			}

			if (mybitcount==0){
				string line;
				getline(fileStream, line);
				continue;
			}

			while (fileStream.peek()!='\n' && fileStream.read(kmerbuffer, sizeof(kmerbuffer))){
				string kmer (kmerbuffer, (kmersize*2/3)+1);
				
				auto it = kmerMap.find(kmer);//check if the kmer is in the map
				if ( it != kmerMap.end() ){//if the kmer is in the map
					for (int i=0; i<fileInclude.size(); ++i){
						it->second[i+includedSampleNum]=mybits[fileInclude[i]];
					}
				}
				else {//if the kmer is not in the hash
					vector < bool > taxonBitset(numSamples, false);//create a new vector of bools for all samples
					for (int i=0; i<fileInclude.size(); ++i){ //change the bools to true based for taxa in the file containing the kmer
						taxonBitset[i+includedSampleNum]=mybits[fileInclude[i]];
					}
					kmerMap.insert(make_pair(kmer, taxonBitset));//add the new vector to the map
				}
	    	}
	    	fileStream.ignore(256,'\n');//skip the end of line character
    	}
    	sampleNum+=numFileSamples; //add the number of samples in the file to the count of total samples
    	includedSampleNum+=fileInclude.size();
    	fileSampleNames.clear();
		fileStream.close();
	}

	cout << kmerMap.size() << " kmers identified from " << includedSampleNum << " samples in " << numfiles << " files" << endl;
	
	cout << "Merging..." << endl;

	unordered_map < vector < bool >,  vector < string > > revKmerMap;

	reverseVectorBoolKmerMap(kmerMap, revKmerMap);//reverse the map so that we have a new map of kmers for each combination of samples

	cout << revKmerMap.size() << " unique taxon combinations in map" << endl;
	
	if(printMergedKmerFile(revKmerMap, outfile, includedSampleNames, kmersize)){return 1;}
	
	printDuration(start);//stop the clock
	
	return 0;
		
}


