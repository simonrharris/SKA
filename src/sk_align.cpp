//g++ -O3 -std=c++0x splitkmer.cpp -lz -o splitkmer
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <cmath>       /* ceil */
#include <iterator>     // std::distance
#include "general.hpp"
#include "kmers.hpp"
#include <chrono> //timing
#include <algorithm>    // std::count
using namespace std;


//int main(int argc, char *argv[])
int alignKmers(const float & minproportion, const string & outputfile, const vector<string> & kmerfiles, const bool & variantonly, const vector <string> & sample)
{

	auto start = chrono::high_resolution_clock::now();//start the clock


	int numfiles=kmerfiles.size();//record the number of files

	vector < string > sampleNames;
	if (collectSampleNames(kmerfiles, sampleNames)!=0){return 1;}//get the total number of samples in all input files

	vector < bool > include (sampleNames.size());
	getSubsample(sample, sampleNames, include);//get a vector of bools representing which samples to include based on a provided sample file. This also removed duplicate samples in the input files

	vector < string > includedSampleNames;
	for (auto it=include.begin(); it!=include.end(); ++it){//make a vector off all included sample names to help printing output files later
		if (*it){
			includedSampleNames.push_back(sampleNames[distance(include.begin(), it)]);
		}
	}

	int numSamples = count(include.begin(), include.end(), true);//count the number of included samples

	cout << numSamples << " samples will be included in the alignment" << endl;

	
	float maxmissing=round((1.0-minproportion)*numSamples);//calculate the maximum number of samples that can be missing for the site to be included in the alignment

	/*if (maxmissing<1){
		maxmissing=1;
	}*/

	float minrequired=numSamples-maxmissing;//calculate the minimum number of samples that must contain a kmer for the site to be included in the alignment
	
	cout << "Keeping variants for which at least " << int(minrequired) << " samples include kmer matches" << endl << endl;
	
	// Create the kmer map
	unordered_map<string, string> kmerMap;
	string emptySequence (numSamples , '-');
	int oldkmersize=0;
	int sampleNum=0;
	int includedSampleNum=0;

	char basebuffer[1];

	ifstream fileStream;
	for (int s = 0; s < kmerfiles.size(); ++s){//read each file and make a map of the kmers which stores the bases for each included sample
		cout << "Reading " << kmerfiles[s] << endl;
		fileStream.open(kmerfiles[s], ios::in);

		if (fileStream.fail()) {
			cout << "Failed to open " << kmerfiles[s] << "\n\n";
			return 0;
		}
		int kmersize;
		vector < string > names;
		try {
			int returnval = readKmerHeader(fileStream, kmersize, names);//read the header from the kmer file to get the kmer size and sample names
		}
		catch (int e){
			cout << "An exception occurred when reading file " << kmerfiles[s] << ". Please check the format. Exception Nr. " << e << '\n';
			return 1;
		}
		
		if (s==0){ //If it's the first file, set the oldkmersize
			oldkmersize=kmersize;
		}

		if (kmersize!=oldkmersize){ //if the file kmer size isn't the same as the other files then ditch out
			cout << "kmer files have different kmer sizes\n\n";
			return 0;
		}

		vector < int > fileInclude;

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
			vector < bool > mybits;
			vectorbool_from_ascii(asciibits, mybits);//convert the ascii representation of the taxa to a vector of bools

			while (fileStream.peek()!='\n' && fileStream.read(basebuffer, sizeof(basebuffer))){
				string base (basebuffer, 1);
				base[0]=toupper(base[0]);
				fileStream.read(kmerbuffer, sizeof(kmerbuffer));
				string kmer (kmerbuffer, kmersize*2/3);
				
				auto it = kmerMap.find(kmer);//check if the kmer is in the map
				if ( it != kmerMap.end() ){//if the kmer is in the map
					for (int i=0; i<fileIncludeSize; ++i){ //add the base to all samples that are true in the bitset
						if (mybits[fileInclude[i]]){
							it->second[i+includedSampleNum]=base[0];
						}
					}
				}
				else {//if the kmer isn't in the map
					if ((includedSampleNum)<=maxmissing){
						string newsequence = emptySequence;//create an empty sequence
						for (int i=0; i<fileInclude.size(); ++i){ //add the base to all included samples that are true in the bitset
							if (mybits[fileInclude[i]]){
								newsequence[i+includedSampleNum]=base[0];
							}
						}
						kmerMap.insert(make_pair(kmer, newsequence));
					}
				}
	    	}
	    	fileStream.ignore(256,'\n');//skip the end ofline character
	    }
	    sampleNum+=int(names.size()); //add the number of samples in the file to the count of total samples
	    includedSampleNum+=fileIncludeSize;
		fileStream.close();
	}
	
	cout << kmerMap.size() << " kmers identified from " << numSamples << " samples in " << numfiles << " files" << endl;
	
	int consta = 0;
	int constc = 0;
	int constg = 0;
	int constt = 0;

	auto it = kmerMap.begin();
	auto endIter = kmerMap.end();

	for (; it!=endIter; ){
		int acgt=0;
		set < char > bases;
		for (int i=0; i<numSamples; ++i){
			switch (it->second[i])
			{
				case 'A':
					acgt++;
					bases.insert('A');
					break;
				case 'C':
					acgt++;
					bases.insert('C');
					break;
				case 'G':
					acgt++;
					bases.insert('G');
					break;
				case 'T':
					acgt++;
					bases.insert('T');
					break;
			}
		}
		if (acgt<minrequired || bases.size()==0){ // (afound+cfound+gfound+tfound)==0){
			kmerMap.erase(it++);
		}
		else if (bases.size()==1 && variantonly){
			for (auto it2=bases.begin(); it2!= bases.end(); ++it2){
				switch (*it2){
					case 'A':
						consta++;
					case 'C':
						constc++;
						break;
					case 'G':
						constg++;
						break;
					case 'T':
						constt++;
						break;
				}
			}
			kmerMap.erase(it++);
		}
		else {
			++it;
		}
		
		
	}
	
	
	cout << "Printing alignment of "<< kmerMap.size() << " sites" << endl;
	
	ofstream alignfile(outputfile);
	
	for (int i=0; i<numSamples; ++i){
		string sampleName =includedSampleNames[i];
		alignfile << ">" << sampleName << "\n";
		float nonns=0.0;
		for (auto it=kmerMap.begin(); it!=kmerMap.end(); ++it){
			alignfile << it->second[i];
			if (it->second[i]=='A' || it->second[i]=='C' || it->second[i]=='G' || it->second[i]=='T'){
				nonns++;
			}
		}
		alignfile << "\n";
		if ((nonns/kmerMap.size())<0.5){
			cout << "Warning: " << sampleName << " only matches " << nonns/kmerMap.size()*100 << "% of kmers in the alignment\n";
		}
	}
	alignfile.close();

	if (variantonly){
		cout << "Constant sites matching filters (a c g t):" << endl;
		cout << consta << " " << constc << " " << constg << " " << constt << endl;
	}

	auto finish = std::chrono::high_resolution_clock::now();
	chrono::duration<double> elapsed = finish - start;
	cout << "Done\n";
	cout << "Total time required: " << elapsed.count() << "s" << endl << endl;
	
	return 0;
	
	
}


