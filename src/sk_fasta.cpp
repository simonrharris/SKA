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

int fastaToKmers(const vector<string> & fastas, const string & outfilename, const long & kmerlen)
{


	auto start = chrono::high_resolution_clock::now();

	// Create the kmer map
	unordered_map<string, array<int, 5> > kmerMap;

	int substringlength=(kmerlen*2)+1;
	int numseqs=0;
	int numbases=0;
	int lastnumbases=0;
	for (int s = 0; s < fastas.size(); ++s){
		cout << "Reading " << fastas[s] << "\n";
		gzFile f = gzopen(fastas[s].c_str(), "r");
		assert(f);
		kseq_t *seq;
		seq = kseq_init(f); // initialize seq 
		while (kseq_read(seq) >= 0) { // read the sequence  
       		//seq->name.s = sequence name
        	//seq->comment.s = comment
        	//seq->seq.s = sequence 
        	//seq->qual.s = quality if present
			
			string sequence(seq->seq.s);//convert the char * sequence to string
			numseqs++;
			
			stringstream sequencestream;
			sequencestream << sequence;//convert the sequence to stringstream
			string subsequence;
			
			
			while (getline(sequencestream, subsequence, 'N')){
				
				if (subsequence.length()<substringlength){
					continue;
				}
				
				transform(subsequence.begin(), subsequence.end(), subsequence.begin(), ::toupper);//change all 	letters in the string to upper case
				
				int i=0;


				for (auto iti = subsequence.cbegin(), end = subsequence.cend()-(substringlength-1); iti != end; ++iti, ++i){
					numbases+=1;
					string kmer=subsequence.substr(i,substringlength);
					
					if (reverse_is_min(kmer, kmerlen+1)){
						reverse(kmer.begin(), kmer.end());
						transform(kmer.begin(),kmer.end(),kmer.begin(),complement);
					}
					
					char base=kmer[kmerlen];
					kmer.erase(kmer.begin()+kmerlen);
					
					auto it = kmerMap.find(kmer);//check if the kmer is in the hash
					
					if ( it != kmerMap.end() ){//if the kmer is in the hash
						it->second[base_score[int(base)]]++;//increment the count of the kmer
					
					}
					else {//if the kmers isn't in the hash
						auto ret = kmerMap.insert(make_pair(kmer, array< int, 5>()));//insert the kmer into the hash
						ret.first->second[base_score[int(base)]]=1;//increment the count of the kmer
					}
				}
		
    			}
		}
		kseq_destroy(seq);
 		gzclose(f);
	}
	
	cout << kmerMap.size() << " unique kmers in map\n";
	
	cout << "Writing kmers to " << outfilename << "\n";
	
	printkmerfile(kmerMap, outfilename, kmerlen);

	auto finish = std::chrono::high_resolution_clock::now();
	chrono::duration<double> elapsed = finish - start;
	cout << "Done\n";
	cout << "Total time required: " << elapsed.count() << "s\n\n";
	
	return 0;
	
	
}


