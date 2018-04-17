//g++ -O3 -std=c++0x src/sk_fastq -lz -o bin/sk_fastq
#include <unordered_map> //std::unordered_map
#include <iostream> //std::cout
#include <fstream> //std::ofstream
#include <algorithm> //std::reverse std::transform
#include <cassert> //std::assert
#include <zlib.h> //required for kseq
#include "kseq.h" //reading fasta/fastq
#include <string> //std::string
#include <sstream> //std::stringstream
#include <vector> //std::vector
#include <array> //std::array
#include <chrono> //timing
#include "kmers.hpp"
//#include <string_view>//Bear in mind for future that string_view allows 'in place' substrings!
KSEQ_INIT(gzFile, gzread)
using namespace std;

//int main(int argc, char *argv[])
int fastqToKmers(const vector<string> & fastqs, const string & outfilename, const int & kmerlen, const int & userminquality, const int & userfilecutoff, const int & usercovcutoff, const float & userminmaf)
{

	auto start = chrono::high_resolution_clock::now();
	
	int minquality=userminquality+33;

	
	// Create the kmer map
	unordered_map<string, array<int,5> > kmerMap;
	int substringlength=(kmerlen*2)+1;
	int numseqs=0;
	int numbases=0;
	int numreadbases=0;
	int lastnumbases=0;
	for (int s = 0; s < fastqs.size(); ++s){ 
		gzFile f = gzopen(fastqs[s].c_str(), "r");
		assert(f);
		kseq_t *seq;
		seq = kseq_init(f); // initialize seq 
		while (kseq_read(seq) >= 0) { // read the sequence  
       		//seq->name.s = sequence name
        	//seq->comment.s = comment
        	//seq->seq.s = sequence 
        	//seq->qual.s = quality if present
			
			numseqs++;
			string sequence(seq->seq.s);//convert the char * sequence to string
			string quality(seq->qual.s);//convert the char * quality to string
			numreadbases+=sequence.length();

			lowqualitytoN(sequence, quality, minquality);
			
			stringstream sequencestream;
			sequencestream << sequence;//convert the sequence to stringstream
			string subsequence;
			
			while (getline(sequencestream, subsequence, 'N')){
				
				if (subsequence.length()<substringlength){
					continue;
					}
				
				transform(subsequence.begin(), subsequence.end(), subsequence.begin(), ::toupper);//change all letters in the string to upper case	
				
				int i=0;
				
				for (auto iti = subsequence.cbegin(), end = subsequence.cend()-(substringlength-1); iti != end; ++iti, ++i){
					numbases+=1;
					string kmer=subsequence.substr(i,substringlength);
					
					if (reverse_is_min(kmer, kmerlen)){
						reverse(kmer.begin(), kmer.end());
						transform(kmer.begin(),kmer.end(),kmer.begin(),complement);
					}
					
					char base=kmer[kmerlen];
					kmer.erase(kmer.begin()+kmerlen);
					
					auto it = kmerMap.find(kmer);//check if the kmer is in the hash
					
					if ( it != kmerMap.end() ){//if the kmer is in the hash
						it->second[base_score[int(base)]]++;//increment the count of the base for the kmer
						if (s>0){//if this is the second file
							it->second[4]++;//increment the count of the kmer
						}
					}
					else if (s==0) {//if the kmers isn't in the hash
						auto ret = kmerMap.insert(make_pair(kmer, array< int, 5>()));//insert the kmer to the hash
						ret.first->second[base_score[int(base)]]=1;//increment the count of the base for the kmer
					}
				}
			}		
    	}
		kseq_destroy(seq);
 		gzclose(f);
		
		if (s==0 && userfilecutoff>0){
			cout << "Added " << numbases << " kmers from " << numseqs << " sequences\n";
			cout << kmerMap.size() << " unique kmers in map\n";
			cout << "Filtering kmers for file coverage\n";

			auto it = kmerMap.begin();
			auto endIter = kmerMap.end();

			for (; it!=endIter; ){
				int basecoverage=0;
				for (auto it2=it->second.begin(); it2!=it->second.end()-1; ++it2){
					basecoverage+=*it2;
				}
				//filter sinlgletons here
				for (auto it2=it->second.begin(); it2!=it->second.end()-1; ++it2){
					if (*it2<userfilecutoff){
						basecoverage-=*it2;
						*it2=0;
					}
				}
				
				if (basecoverage==0){
						kmerMap.erase(it++);
					}
				else{
					++it;
				}
				
			}
			cout << kmerMap.size() << " unique kmers in map after filtering\n";
		}
		
		
	}
	cout << "Added " << numbases << " kmers from " << numseqs << " sequences of total length " << numreadbases << "\n";
	cout << "Running final filtering of kmers\n";

	int totalcoverage=0;

	auto it = kmerMap.begin();
	auto endIter = kmerMap.end();
	for (; it!=endIter; ){

		//remove kmers with file coverage lower than the user defined cutoff if there were two fastq files supplied
		if (fastqs.size()>1 && it->second[4]<userfilecutoff){
			kmerMap.erase(it++);
			continue;
		}
		else{
			it->second[4]=0;
		}

		//calculate the total kmer coverage
		int basecoverage=0;
		for (auto it2=it->second.begin(); it2!=it->second.end()-1; ++it2){
			basecoverage+=*it2;
		}

		float covcutoff=float(basecoverage)*userminmaf;
		if (covcutoff<usercovcutoff){//change 8 here to change minimum acceptable coverage
			covcutoff=usercovcutoff;
		}

		//filter bases that don't meet the coverage cutoff
		for (auto it2=it->second.begin(); it2!=it->second.end()-1; ++it2){
			if (*it2<covcutoff){
				basecoverage-=*it2;
				*it2=0;
			}
		}

		//remove kmers with coverage lower than the user defined cutoff
		if (basecoverage<usercovcutoff){
				kmerMap.erase(it++);
			}
		else{
			totalcoverage+=basecoverage;
			++it;
		}
		
	}
	cout << kmerMap.size() << " unique kmers in map after final filtering\n";
	cout << "Mean kmer coverage is " << float(totalcoverage)/kmerMap.size() << "\n";
	
	cout << "Writing kmers to fasta file\n";
	
	printkmerfile(kmerMap, outfilename, kmerlen);


	auto finish = std::chrono::high_resolution_clock::now();
	chrono::duration<double> elapsed = finish - start;
	cout << "Done\n";
	cout << "Total time required: " << elapsed.count() << "s\n";
	
	return 0;
	
	/*
	bool gzoutput=false;
	
	if (gzoutput){
	gzFile *fi = (gzFile *)gzopen(argv[1],"wb");
	for (auto it=kmerMap.begin(); it!=kmerMap.end(); ++it){
		
		string kmer=it->first;
		gzwrite(fi,">",strlen(">"));
		if (it->second.size()>1){
			gzwrite(fi, "N", strlen("N"));
		}
		else{
			string base={it->second.begin()->first};
			gzwrite(fi, base.c_str(), strlen(base.c_str()));
		}
		gzwrite(fi,"\n",strlen("\n"));
		ascii_codons(kmer);
		gzwrite(fi,kmer.c_str(),strlen(kmer.c_str()));
		gzwrite(fi,"\n",strlen("\n"));
		
	}
	gzclose(fi);
	}
	else {
	FILE *f = fopen(argv[1],"wb");
	for (auto it=kmerMap.begin(); it!=kmerMap.end(); ++it){
		
		string kmer=it->first;
			
		fwrite(">", sizeof(char), strlen(">"), f);
		if (it->second.size()>1){
			fwrite("N", sizeof(char), strlen("N"), f);//Could as IUPAC codes here?
		}
		else{
			string base={it->second.begin()->first};
			fwrite(base.c_str(), sizeof(char), strlen(base.c_str()), f);
		}
		fwrite("\n", sizeof(char), strlen("\n"), f);
		ascii_codons(kmer);
		fwrite(kmer.c_str(), sizeof(char), strlen(kmer.c_str()), f);
		fwrite("\n", sizeof(char), strlen("\n"), f);
		
	}
	fclose(f);
	}
	
	ofstream kmerfile(argv[1]);
	for (auto it=kmerMap.begin(); it!=kmerMap.end(); ++it){
	
		string kmer=it->first;
		string base;
		int nonzero=0;
		int i=0;
		for (auto it2=it->second.begin(); it2!=it->second.end(); ++it2, ++i){
		
			if (*it2>0){
				if (nonzero>0){
					base='N';
					break;
				}
				else{
					base=bases[i];
				}
			}
		
		}
		kmerfile << base << " " << kmer <<"\n";
		
	
	}
	kmerfile.close();
	
	return 0;*/
	
	
}


