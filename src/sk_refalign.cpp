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
#include "kseq.h"
#include <string>
#include <tuple>
#include <vector>
#include "kmers.hpp"
#include "general.hpp"
#include <chrono> //timing
KSEQ_INIT(gzFile, gzread)
using namespace std;


//int main(int argc, char *argv[])
int alignKmersToReference(const string & reference, const string & outputfile, const vector<string> & kmerfiles, const int & kmerlen, const bool & includeref, const bool & maprepeats)
{

	auto start = chrono::high_resolution_clock::now();
	
	ofstream alignfile(outputfile);
	
	// Create the kmer map
	unordered_map<string, vector<int> > kmerMap;
	set<int> revSet;
	int substringlength=(kmerlen*2)+1;
	int numfiles=0;
	++numfiles;
	int numseqs=0;
	string filename=splitFileName(reference);
	cout << "Reading " << filename << "\n";
	if (includeref){
		alignfile << ">" << filename << "\n";
	}

	gzFile f = gzopen(reference.c_str(), "r");
	assert(f);
	kseq_t *seq;
	seq = kseq_init(f); // initialize seq 
	int basenum=0;
	string refseq;
	while (kseq_read(seq) >= 0) { // read the sequence  
       		//seq->name.s = sequence name
        	//seq->comment.s = comment
        	//seq->seq.s = sequence 
        	//seq->qual.s = quality if present
		
		string sequence(seq->seq.s);//convert the char * sequence to string
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
	kseq_destroy(seq);
 	gzclose(f);
	if (includeref){
		alignfile << "\n";
	}
	
	cout << kmerMap.size() << " kmers read\n";
	
	ifstream fileStream;
	char basebuffer[1];
	char kmerbuffer[kmerlen*2/3];
	
	for (int s = 0; s < kmerfiles.size(); ++s){
		++numfiles;
		string newSequence (basenum , '-');
		int mappedkmers=0;
		int unmappedkmers=0;

		string filename=splitFileName(kmerfiles[s]);
		cout << "Aligning " << filename;
		alignfile << ">" << filename << "\n";
		fileStream.open(kmerfiles[s], ios::in);

		if (fileStream.fail()) {
			cout << "\nFailed to open" << kmerfiles[s] << "\n\n";;
			return 0;
		}
		int kmersize=readKmerHeader(fileStream);
		if (kmersize!=kmerlen){
			cout << "\nkmer size in " << filename << " is not " << kmerlen << "\n\n";
			return 0;
		}

		while (fileStream.read(basebuffer, sizeof(basebuffer))){
			string base (basebuffer, 1);
			base[0]=toupper(base[0]);
			fileStream.read(kmerbuffer, sizeof(kmerbuffer));
			string kmer (kmerbuffer, kmerlen*2/3);
			auto it = kmerMap.find(kmer);//check if the kmer is in the hash
			if ( it != kmerMap.end() ){//if the kmer is in the hash
				mappedkmers++;
				if (it->second.size()==1 || maprepeats){
					for (auto itb = it->second.begin(); itb != it->second.end(); ++itb) {
						auto itc = revSet.find(*itb);
						if ( itc != revSet.end() ){
							transform(base.begin(),base.end(),base.begin(),complement);
							newSequence[*itb]=base[0];
						}
						else {
							newSequence[*itb]=base[0];//
						}
						if (newSequence[*itb]!=refseq[*itb] || *itb==15 || *itb==basenum-16){
							for (int i=*itb-15; i<=*itb+15; ++i){
								if (i==*itb){
									continue;
								}
								newSequence[i]=refseq[i];
							}
						}
						/*if (newSequence[*itb]!=refseq[*itb]){
							cout << *itb+1 << "\t" << refseq[*itb] << "\t" << newSequence[*itb] << "\n";
						}*/
						
					}
				}
				else{
					unmappedkmers++;
				}
			}
				
		
    	}

    	int mappedbases=basenum;
    	int mappedNs=0;
    	for (int i=0; i<basenum; ++i){
    		if (newSequence[i]=='-'){
    			mappedbases--;
    		}
    		else if (newSequence[i]=='N'){
    			mappedbases++;
    		}
    	}

		fileStream.close();
		//cout << mappedkmers << "/" << unmappedkmers+mappedkmers << " (" << float(mappedkmers)/(unmappedkmers+mappedkmers)*100 << "%) kmers mapped";
		//cout << mappedbases << "/" << basenum << " (" << float(mappedbases)/(basenum)*100 << "%) reference bases mapped";
		cout << " ... " << float(mappedkmers)/(unmappedkmers+mappedkmers)*100 << "% of kmers mapped to ";
		cout << float(mappedbases)/(basenum)*100 << "% of reference bases\n";
		alignfile << newSequence << "\n";
	}
	alignfile.close();
	
	auto finish = std::chrono::high_resolution_clock::now();
	chrono::duration<double> elapsed = finish - start;
	cout << "Done\n";
	cout << "Total time required: " << elapsed.count() << "s\n\n";

	return 0;
	
	
}


