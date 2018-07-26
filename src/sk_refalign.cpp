//g++ -O3 -std=c++0x splitkmer.cpp -lz -o splitkmer
//usage: sk_refalign outputfilename reference.fasta reference.kmers <kmer files>
#include <unordered_map>
#include <set>
#include <iostream>
#include <fstream>
#include <cassert>
#include <algorithm>
#include <zlib.h>
#include <stdio.h>
//#include "kseq.h"
#include <string>
#include <tuple>
#include <vector>
#include <cmath>       /* ceil */
#include "kmers.hpp"
#include "general.hpp"
#include <chrono> //timing
#include "gzstream.h"
//KSEQ_INIT(gzFile, gzread)
using namespace std;


//int main(int argc, char *argv[])
int alignKmersToReference(const string & reference, const string & outputfile, const vector<string> & kmerfiles, const int & kmerlen, const bool & includeref, const bool & maprepeats, const bool &fillall, const bool &variantonly, const vector <string> & sample)
{

	auto start = chrono::high_resolution_clock::now();
	
	ofstream alignfile(outputfile);
	
	// Create the kmer map
	unordered_map<string, vector<int> > kmerMap;
	set<int> revSet;
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

	igzstream gzfileStream;
	gzfileStream.open(reference.c_str());

	if (gzfileStream.get()!='>'){
		cout << reference << " is not in the correct format. Expecting header to start with >." << endl << endl;
	}

	while (getline(gzfileStream, header)){
		getline(gzfileStream, sequence, '>');
		sequence.erase(std::remove_if( sequence.begin(), sequence.end(), ::isspace ), sequence.end() );
		transform(sequence.begin(), sequence.end(), sequence.begin(), ::toupper);//change all letters in the string to upper case

		refseq.append(sequence);
		if (includeref){
			alignfile << sequence;
		}
		numseqs++;
		
		int i=0;
		
		for (auto iti = sequence.cbegin(), end = sequence.cend()-(substringlength-1); iti != end; ++iti, ++i){
			string kmer=sequence.substr(i,substringlength);
			
			if (reverse_is_min(kmer, kmerlen+1)){
				reverse(kmer.begin(), kmer.end());
				transform(kmer.begin(),kmer.end(),kmer.begin(),complement);
				revSet.insert(basenum+i+kmerlen);
			}
			
			char base=kmer[kmerlen];
			kmer.erase(kmer.begin()+kmerlen);
			ascii_codons(kmer);
			
			auto it = kmerMap.find(kmer);//check if the kmer is in the hash
			if ( it != kmerMap.end() ){//if the kmer is in the hash
				it->second.push_back(basenum+i+kmerlen);//add the location of the match to the map
			}
			else {
				auto ret = kmerMap.insert(make_pair(kmer, vector<int>()));
				ret.first->second.push_back(basenum+i+kmerlen);
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
	for (auto it=include.begin(); it!=include.end(); ++it){//make a vector off all included sample names to help printing output files later
		if (*it){
			includedSampleNames.push_back(sampleNames[distance(include.begin(), it)]);
		}
	}
	
	int numSamples = count(include.begin(), include.end(), true);//count the number of included samples

	cout << numSamples << " samples will be aligned to " << reference << endl;

	vector < string > sequences(numSamples, string(basenum , '-'));//create a vector to store the new sequences;

	ifstream fileStream;
	char basebuffer[1];
	char kmerbuffer[kmerlen*2/3];
	int sampleNum=0;
	int includedSampleNum=0;
	vector < string > names;
	
	for (int s = 0; s < kmerfiles.size(); ++s){

		string filename=splitFileName(kmerfiles[s]);
		cout << "Aligning " << filename << endl;
		fileStream.open(kmerfiles[s], ios::in);

		if (fileStream.fail()) {
			cout << "\nFailed to open" << kmerfiles[s] << endl << endl;
			return 0;
		}
		int kmersize;
		
		
		try {
			int returnval = readKmerHeader(fileStream, kmersize, names);//read the header from the kmer file to get the kmer size and sample names
		}
		catch (int e){
			cout << "An exception occurred when reading file " << kmerfiles[s] << ". Please check the format. Exception Nr. " << e << '\n';
			return 1;
		}
		if (kmersize!=kmerlen){
			cout << "\nkmer size in " << filename << " is not " << kmerlen << endl << endl;
			return 0;
		}

		vector < int > fileInclude;

		for (int i=0; i<names.size(); ++i){ //put the index of all sample names that are to be included into a vector
			if (include[sampleNum+i]){
				fileInclude.push_back(i);
			}
		}

		if (fileInclude.size()==0){ // if no sample names in the file are going to be included then don't read the file
			fileStream.close();
			cout << "Nothing to align" << endl;
			sampleNum+=int(names.size());
			names.clear();
			continue;
		}


		char asciibuffer[int(ceil(float(names.size())/6))];

		while (fileStream.read(asciibuffer, sizeof(asciibuffer))){
			string asciibits (asciibuffer, sizeof(asciibuffer));
			vector < bool > mybits;
			vectorbool_from_ascii(asciibits, mybits);//read the ascii representation of the taxa to a vector of bools
			
			while (fileStream.peek()!='\n' && fileStream.read(basebuffer, sizeof(basebuffer))){
				string base (basebuffer, 1);
				base[0]=toupper(base[0]);
				fileStream.read(kmerbuffer, sizeof(kmerbuffer));
				string kmer (kmerbuffer, kmerlen*2/3);
				auto it = kmerMap.find(kmer);//check if the kmer is in the hash
				if ( it != kmerMap.end() ){//if the kmer is in the hash
					if (it->second.size()==1 || maprepeats){
						for (auto itb = it->second.begin(); itb != it->second.end(); ++itb) {
							auto itc = revSet.find(*itb);
							if ( itc != revSet.end() ){
								transform(base.begin(),base.end(),base.begin(),complement);
							}
							for (int j=0; j<fileInclude.size(); ++j){ //add the base to all samples that are true in the bitset
								if (mybits[fileInclude[j]]){
									sequences[includedSampleNum+j][*itb]=base[0];
									if (fillall || base[0]!=refseq[*itb]){//} || *itb==15 || *itb==basenum-16){
										for (int i=*itb-kmerlen; i<=*itb+kmerlen; ++i){
											if (i==*itb){
												continue;
											}
											if (sequences[includedSampleNum+j][i]=='-'){
												sequences[includedSampleNum+j][i]=tolower(refseq[i]);
											}
											/*else if ((sequences[sampleNum+j][i]=='a' || sequences[sampleNum+j][i]=='c' || sequences[sampleNum+j][i]=='g' ||sequences[sampleNum+j][i]=='t') && sequences[sampleNum+j][i]!=tolower(refseq[i])){
												sequences[sampleNum+j][i]='n';
											}*/
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
		int consta = 0;
		int constc = 0;
		int constg = 0;
		int constt = 0;
		for (int i=0; i<basenum; ++i){
			int a = 0;
			int c = 0;
			int g = 0;
			int t = 0;
			int afound = 0;
			int cfound = 0;
			int gfound = 0;
			int tfound = 0;
			for (int j=0; j<includedSampleNames.size(); ++j){
				switch (toupper(sequences[j][i]))
				{
					case 'A':
						a++;
						afound=1;
						break;
					case 'C':
						c++;
						cfound=1;
						break;
					case 'G':
						g++;
						gfound=1;
						break;
					case 'T':
						t++;
						tfound=1;
						break;
				}
			}
			if ((afound+cfound+gfound+tfound)==1){
				if (a>0){
					consta++;
				}
				else if (c>0){
					constc++;
				}
				else if (g>0){
					constg++;
				}
				else if (t>0){
					constt++;
				}
			}
			else if ((afound+cfound+gfound+tfound)>1) {
				sitestoprint.push_back(i);
			}
		}

		cout << "Printing alignment of " << sitestoprint.size() << " variant sites" << endl;

		for (int i=0; i<includedSampleNames.size(); ++i){
			alignfile << ">" << includedSampleNames[i] << endl;
			for (int j=0; j<sitestoprint.size(); ++j){
				alignfile << sequences[i][sitestoprint[j]];
			}
			alignfile << endl;
		}

		cout << "Constant sites (a c g t):" << endl;
		cout << consta << " " << constc << " " << constg << " " << constt << endl;

	}
	else{
		cout << "Printing alignment" << endl;
		for (int i=0; i<includedSampleNames.size(); ++i){
			alignfile << ">" << includedSampleNames[i] << endl;
			alignfile << sequences[i] << endl;
		}
	}
	


	alignfile.close();
	
	auto finish = std::chrono::high_resolution_clock::now();
	chrono::duration<double> elapsed = finish - start;
	cout << "Done" << endl;
	cout << "Total time required: " << elapsed.count() << "s" << endl << endl;

	return 0;
	
	
}


