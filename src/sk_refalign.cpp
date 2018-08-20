#include <unordered_map>
#include <set>
#include <iostream>
#include <fstream>
#include <algorithm>
//#include <stdio.h>
#include <string>
#include <vector>
#include <cmath>       // ceil
#include "kmers.hpp"
#include "general.hpp"
#include "DNA.hpp"
#include <chrono> //timing
#include "gzstream.h"
using namespace std;


int alignKmersToReference(const string & reference, const string & outputfile, const vector<string> & kmerfiles, const int & kmerlen, const bool & includeref, const bool & maprepeats, const bool &fillall, const bool &variantonly, const vector <string> & sample)
{

	const chrono::high_resolution_clock::time_point start = chrono::high_resolution_clock::now();
	
	ofstream alignfile(outputfile);
	if (alignfile.fail()){
		cerr << endl << "Error: Failed to open " << outputfile << endl << endl;
		return 1;
	}
	
	// Create the kmer map
	unordered_map < string, vector <int> > kmerMap;
	set <int> revSet;
	int substringlength=(kmerlen*2)+1;
	int numseqs=0;
	string filename=splitFileName(reference);
	cout << "Reading " << filename << endl;
	if (includeref){
		alignfile << ">" << filename << endl;
	}

	string sequence;
	string header;
	int basenum=0;
	string refseq;
	char base;

	igzstream gzfileStream;
	gzfileStream.open(reference.c_str());
	if (gzfileStream.fail()){
		cerr << endl << "Error: Failed to open " << reference << endl << endl;
		return 1;
	}

	if (gzfileStream.get()!='>'){
		cout << reference << " is not in the correct format. Expecting header to start with >." << endl << endl;
	}

	while (getline(gzfileStream, header)){
		getline(gzfileStream, sequence, '>');
		sequence.erase(remove_if( sequence.begin(), sequence.end(), ::isspace ), sequence.end() );
		transform(sequence.begin(), sequence.end(), sequence.begin(), ::toupper);//change all letters in the string to upper case

		refseq.append(sequence);
		if (includeref){
			alignfile << sequence;
		}
		numseqs++;
		
		int i=0;
		
		for (string::iterator iti = sequence.begin(), end = sequence.end()-(substringlength-1); iti != end; ++iti, ++i){
			string kmer=sequence.substr(i,substringlength);
			
			if (reverseComplementIfMin(kmer)){
				revSet.insert(basenum+i+kmerlen);
			}
			extractMiddleBase(kmer, base);
			ascii_codons(kmer);
			
			unordered_map < string, vector <int> >::iterator it = kmerMap.find(kmer);//check if the kmer is in the hash
			if ( it != kmerMap.end() ){//if the kmer is in the hash
				it->second.push_back(basenum+i+kmerlen);//add the location of the match to the map
			}
			else {
				kmerMap.insert(make_pair(kmer, vector<int> {basenum+i+kmerlen}));
			}
		}
		basenum+=sequence.length();
			
    }
    gzfileStream.close();

	if (includeref){
		alignfile << endl;
	}
	
	cout << kmerMap.size() << " kmers read" << endl;
	
	vector < string > sampleNames;
	if (collectSampleNames(kmerfiles, sampleNames)!=0){return 1;}//get sample names from all kmerfiles

	vector < bool > include (sampleNames.size());
	getSubsample(sample, sampleNames, include);//get a vector of bools representing which samples to include based on a provided sample file. This also removed duplicate samples in the input files

	vector < string > includedSampleNames;
	for (vector < bool >::iterator it=include.begin(); it!=include.end(); ++it){//make a vector off all included sample names to help printing output files later
		if (*it){
			includedSampleNames.push_back(sampleNames[distance(include.begin(), it)]);
		}
	}
	
	int numSamples = count(include.begin(), include.end(), true);//count the number of included samples

	cout << numSamples << " samples will be aligned to " << reference << endl;

	vector < string > sequences(numSamples, string(basenum , '-'));//create a vector to store the new sequences;

	ifstream fileStream;

	char kmerbuffer[kmerlen*2/3];
	int sampleNum=0;
	int includedSampleNum=0;
	vector < string > names;
	
	for (int s = 0; s < kmerfiles.size(); ++s){

		string filename=splitFileName(kmerfiles[s]);
		if (openFileStream(kmerfiles[s], fileStream)){return 1;};

		int kmersize;
		
		readKmerHeader(fileStream, kmersize, names);//read the header from the kmer file to get the kmer size and sample names
		
		if (kmersize!=kmerlen){
			cerr << endl << "kmer size in " << filename << " is not " << kmerlen << endl << endl;
			return 1;
		}

		vector < int > fileInclude;

		for (int i=0; i<names.size(); ++i){ //put the index of all sample names that are to be included into a vector
			if (include[sampleNum+i]){
				fileInclude.push_back(i);
			}
		}

		if (fileInclude.size()==0){ // if no sample names in the file are going to be included then don't read the file
			fileStream.close();
			cerr << "Nothing to align" << endl;
			sampleNum+=int(names.size());
			names.clear();
			continue;
		}

		char asciibuffer[int(ceil(float(names.size())/6))];

		while (fileStream.read(asciibuffer, sizeof(asciibuffer))){
			string asciibits (asciibuffer, sizeof(asciibuffer));
			vector < bool > mybits;
			vectorbool_from_ascii(asciibits, mybits);//read the ascii representation of the taxa to a vector of bools
			
			while (fileStream.peek()!='\n' && fileStream.get(base)){
				base=toupper(base);
				fileStream.read(kmerbuffer, sizeof(kmerbuffer));
				string kmer (kmerbuffer, kmerlen*2/3);
				unordered_map < string, vector <int> >::iterator it = kmerMap.find(kmer);//check if the kmer is in the hash
				if ( it != kmerMap.end() ){//if the kmer is in the hash
					if (it->second.size()==1 || maprepeats){
						for (vector < int >::iterator itb = it->second.begin(); itb != it->second.end(); ++itb) {
							set < int >::iterator itc = revSet.find(*itb);
							if ( itc != revSet.end() ){
								base=complement(base);
							}
							for (int j=0; j<fileInclude.size(); ++j){ //add the base to all samples that are true in the bitset
								if (mybits[fileInclude[j]]){
									sequences[includedSampleNum+j][*itb]=base;
									if (fillall || base!=refseq[*itb]){//} || *itb==15 || *itb==basenum-16){
										for (int i=*itb-kmerlen; i<=*itb+kmerlen; ++i){
											if (i==*itb){
												continue;
											}
											if (sequences[includedSampleNum+j][i]=='-'){
												sequences[includedSampleNum+j][i]=tolower(refseq[i]);
											}
										}
									}
								}
							}
						}
					}
				}
			}
			fileStream.ignore(256,'\n');//skip the end ofline character
    	}

    	for (int i=0; i<fileInclude.size(); ++i){
	    	int mappedbases=basenum;
	    	int mappedNs=0;
	    	for (int j=0; j<basenum; ++j){
	    		if (sequences[includedSampleNum+i][j]=='-'){
	    			mappedbases--;
	    		}
	    	}
			cout << includedSampleNames[i] << ": " << float(mappedbases)/(basenum)*100 << "% of reference bases mapped" << endl;
    	}
    	sampleNum+=int(names.size()); //add the number of samples in the file to the count of total samples
    	includedSampleNum+=fileInclude.size();
    	names.clear();
		fileStream.close();
	}

	
	if (variantonly){
		vector < int > sitestoprint;
		sitestoprint.reserve(10000);
		vector < int > constantBases (4,0);
		for (int i=0; i<basenum; ++i){

			vector < int > baseVector (6, 0);

			for (int j=0; j<includedSampleNames.size(); ++j){
				char base=toupper(sequences[j][i]);
				baseVector[base_score[base]]++;
			}

			int acgtCount=0;
			int constbase=-1;

			for (int i=0; i<4; ++i){
				if (baseVector[i]>0){
					acgtCount++;
					constbase=i;
				}
			}

			if (acgtCount<2){
				constantBases[constbase]++;
			}
			else {
				sitestoprint.push_back(i);
			}
			
		}

		cout << "Printing alignment of " << sitestoprint.size() << " variant sites" << endl;

		for (int i=0; i<includedSampleNames.size(); ++i){
			string mySequence (sitestoprint.size(), '-');
			for (int j=0; j<sitestoprint.size(); ++j){
				mySequence[j]=sequences[i][sitestoprint[j]];
			}
			alignfile << ">" << includedSampleNames[i] << endl << mySequence << endl;
		}

		cout << "Constant sites (a c g t):" << endl;
		cout << constantBases[0] << " " << constantBases[1]  << " " << constantBases[2]  << " " << constantBases[3]  << endl;
	}
	else{
		cout << "Printing alignment" << endl;
		for (int i=0; i<includedSampleNames.size(); ++i){
			alignfile << ">" << includedSampleNames[i] << endl << sequences[i] << endl;
		}
	}

	alignfile.close();
	
	printDuration(start);

	return 0;
	
	
}


