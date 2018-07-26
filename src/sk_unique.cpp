//g++ -O3 -std=c++0x src/sk_weed.cpp -lz -o bin/sk_weed
#include <unordered_map> //std::unordered_map
#include <iostream> //std::cout
#include <fstream> //std::ofstream
#include <string> //std::string
#include <vector> //std::vector
#include <set> //std::set
#include <sstream>
#include <chrono> //timing
#include <algorithm> //std::transform
#include "kmers.hpp"

using namespace std;

//bin/sk_weed outputfile 

//int main(int argc, char *argv[])
int uniqueKmers(const vector <string> & ingroupsamples, const vector<string> & kmerfiles, const float & minproportion, const string & outputfile)
{

	auto start = chrono::high_resolution_clock::now();
	// Create the kmer map
	unordered_map<string, vector < bool > > kmerMap;
	int oldkmersize=0;
	int totalingroupkmers=0;

	vector < string > sampleNames;

	if (collectSampleNames(kmerfiles, sampleNames)!=0){return 1;}
	int numSamples=sampleNames.size();

	float maxmissing=(1.0-minproportion)*ingroupsamples.size();

	float minrequired=ingroupsamples.size()-maxmissing;

	set < string > ingroupSet(ingroupsamples.begin(), ingroupsamples.end());

	if (ingroupSet.size()!=ingroupsamples.size()){
		cout << "Warning: Duplicate samples in your ingroup sample file have been removed" << endl;
	}

	vector < int > ingroupSamplePositions;

	for (int i=0; i<sampleNames.size(); ++i){
		auto it = ingroupSet.find(sampleNames[i]);//check if the kmer is in the map
			if ( it != ingroupSet.end() ){//if the kmer is in the map
				ingroupSamplePositions.push_back(i);
			}
	}
	
	if (ingroupsamples.size()!=ingroupSamplePositions.size()){
		cout << ingroupsamples.size() << " " << ingroupSamplePositions.size() << endl;
		cout << "Error: Some of your ingroup samples are not in your input files" << endl << endl;
		return 1;
	}

	vector < string > fileSampleNames;	
	int sampleNum=0;
	int numfiles=kmerfiles.size();
	ifstream fileStream;
	int kmersize;
	char * kmerbuffer;


	for (int s = 0; s < kmerfiles.size(); ++s){

		cout << "Reading " << kmerfiles[s] << endl;
		fileStream.open(kmerfiles[s], ios::in);

		if (fileStream.fail()) {
			cout << "Failed to open " << kmerfiles[s] << "\n" << endl;
			return 0;
		}
		int newkmersize;
		try {
			int returnval = readKmerHeader(fileStream, newkmersize, fileSampleNames);
		}
		catch (int e){
			cout << "An exception occurred when reading file " << kmerfiles[s] << ". Please check the format. Exception Nr. " << e << '\n';
			return 1;
		}
		if (s==0){
			kmersize=newkmersize;
		}

		if (newkmersize!=kmersize){
			cout << "kmer files have different kmer sizes\n" << endl;
			return 0;
		}


		int numFileSamples=fileSampleNames.size();

		char basebuffer[1];
		char kmerbuffer[(kmersize*2/3)+1];
		char asciibuffer[int(ceil(float(fileSampleNames.size())/6))];

		while (fileStream.read(asciibuffer, sizeof(asciibuffer))){//read the ascii representation of the taxa
			string asciibits (asciibuffer, sizeof(asciibuffer));
			vector < bool > mybits;
			vectorbool_from_ascii(asciibits, mybits);//read the ascii representation of the taxa to a vector of bools

			while (fileStream.peek()!='\n' && fileStream.read(kmerbuffer, sizeof(kmerbuffer))){
				string kmer (kmerbuffer, (kmersize*2/3)+1);
				

				auto it = kmerMap.find(kmer);//check if the kmer is in the map
				if ( it != kmerMap.end() ){//if the kmer is in the map
					for (int i=0; i<mybits.size(); ++i){
						it->second[i+sampleNum]=mybits[i];
					}
				}
				else {//if the kmer is not in the hash
					vector < bool > taxonBitset(numSamples, false);//create a new vector of bools for all samples
					for (int i=0; i<mybits.size(); ++i){ //change the bools to true based for taxa in the file containing the kmer
						taxonBitset[i+sampleNum]=mybits[i];
					}
					kmerMap.insert(make_pair(kmer, taxonBitset));//add the new vector to the map
				}
	    	}
	    	fileStream.ignore(256,'\n');//skip the end ofline character
    	}
    	sampleNum+=numFileSamples; //add the number of samples in the file to the count of total samples
    	fileSampleNames.clear();
		fileStream.close();
	}

	cout << kmerMap.size() << " kmers identified from " << sampleNum << " samples in " << numfiles << " files" << endl;

	int uniquecount=0;
	auto it = kmerMap.begin();
	auto endIter = kmerMap.end();
	unordered_map < vector < bool >,  vector < string > > revKmerMap;

	for (; it!=endIter; ){
		int allcount=count(it->second.begin(), it->second.end(), true);
		int ingroupcount=0;

		for (auto it2=ingroupSamplePositions.begin(); it2!=ingroupSamplePositions.end(); ++it2){
			if (it->second[*it2]){
				ingroupcount++;
			}
		}
		if (ingroupcount>=minrequired && (allcount-ingroupcount)==0){
			uniquecount++;

			vector < bool > ingroupBitString (ingroupSamplePositions.size(), false);
			int i=0;
			for (auto it2=ingroupSamplePositions.begin(); it2!=ingroupSamplePositions.end(); ++it2, ++i){
				ingroupBitString[i]=it->second[*it2];
			}

			auto it2 = revKmerMap.find(ingroupBitString);//check if the bitset is in the map
			if ( it2 != revKmerMap.end() ){//if the bitset is in the map
				it2->second.push_back(it->first); //add the kmer
			}
			else {//if the bitset isn't in the map
				vector <string> myvector; //create a new vector of strings
				myvector.push_back(it->first); //add the bitset to the vector
				revKmerMap.insert(make_pair(ingroupBitString, myvector)); //insert the vector into the map
			}
		}
		kmerMap.erase(it++); //remove kmer from kmerMap to save space
	}
	if (kmerMap.size()!=0){
		kmerMap.clear(); //clear the kmerMap. This should already be empty.
		cout << "Why wasn't the kmerMap empty?" << endl; 
	}

	cout << uniquecount << " unique shared kmers will be written to " << outputfile << endl;

	//This is identical to code in sk_merge. Move to kmers.cpp?

	ofstream kmerout(outputfile); //open output file stream 
	kmerout << kmersize << endl; // print kmer size to output file stream
	for ( auto it=ingroupSamplePositions.begin(); it!=ingroupSamplePositions.end(); ++it){ //print each sample name to output file stream
		kmerout << sampleNames[*it] << " "; 
	}
	kmerout << endl;

	for ( auto it=revKmerMap.begin(); it!=revKmerMap.end(); ++it){
		stringstream bitstringstream;
		for (auto it2=it->first.begin(); it2!=it->first.end(); ++it2){
			bitstringstream << *it2;
		}
		string bitstring = bitstringstream.str();
		int myremainder=::fmod(int(bitstring.length()),6);
		if (myremainder>0){
			for (int i = 0; i<(6-myremainder); ++i){
				bitstring.push_back('0');
			}
		}
		ascii_bitstring(bitstring);
		kmerout << bitstring;
		for (auto it2=it->second.begin(); it2!=it->second.end(); ++it2){
			kmerout << *it2;
		}
		kmerout << endl;
	}
	
	kmerout.close();


	auto finish = std::chrono::high_resolution_clock::now();
	chrono::duration<double> elapsed = finish - start;
	cout << "Done\n";
	cout << "Total time required: " << elapsed.count() << "s\n\n";
	
	return 0;
	
	
}


