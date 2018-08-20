//g++ -O3 -std=c++0x src/sk_weed.cpp -lz -o bin/sk_weed
#include <unordered_map> //std::unordered_map
#include <iostream> //std::cout
#include <fstream> //std::ofstream
#include <string> //std::string
#include <vector> //std::vector
#include <chrono> //timing
#include <cmath>       /* ceil */
#include "kmers.hpp"
#include "general.hpp"
#include "DNA.hpp"

using namespace std;

//bin/sk_weed outputfile 

//int main(int argc, char *argv[])
int weedKmers(const vector <string> & weedfiles, const string & kmerfile)
{

	const chrono::high_resolution_clock::time_point start = chrono::high_resolution_clock::now();
	// Create the kmer map
	unordered_map < string, char > kmerMap;
	
	ifstream fileStream;

	if (openFileStream(kmerfile, fileStream)){return 1;};

	int kmersize;
	vector < string > weednames;
	
	readKmerHeader(fileStream, kmersize, weednames);
	
	char base;
	char kmerbuffer[kmersize*2/3];
	char asciibuffer[int(ceil(float(weednames.size())/6))];
	while (fileStream.read(asciibuffer, sizeof(asciibuffer))){
		string asciibits (asciibuffer, sizeof(asciibuffer));
		
		while (fileStream.peek()!='\n' && fileStream.get(base)){
			base=toupper(base);
			fileStream.read(kmerbuffer, sizeof(kmerbuffer));
			string kmer (kmerbuffer, sizeof(kmerbuffer));
			kmerMap.insert(make_pair(kmer, base));
		}
		fileStream.ignore(256,'\n');//skip the end ofline character
    }
	fileStream.close();
	
	cout << kmerMap.size() << " unique kmers in map" << endl;
	
	int weeded=0;
	int totalweedkmers=0;
	
	for (vector < string >::const_iterator it = weedfiles.begin(); it != weedfiles.end(); ++it){

		if (openFileStream(*it, fileStream)){return 1;};

		int weedkmersize;
		vector < string > names;
		
		readKmerHeader(fileStream, weedkmersize, names);
		
		if (weedkmersize!=kmersize){
			cout << "Files have different kmer sizes" << endl << endl;
			return 0;
		}
		string filename=string(*it);
		string outputfile=filename.substr(0, filename.find_last_of("."))+".weeded"+filename.substr(filename.find_last_of("."), filename.length());

		cout << "Writing weeded kmers to " << outputfile << endl;
	
		ofstream weededfile(outputfile);
		weededfile << kmersize << "\n";
		for ( vector < string >::iterator it=names.begin(); it!=names.end(); ++it){ //print each sample name to output file stream
			weededfile << *it << " "; 
		}
		weededfile << endl;

		bool firstkmer=true;

		char asciibuffer[int(ceil(float(names.size())/6))];
		while (fileStream.read(asciibuffer, sizeof(asciibuffer))){
			string asciibits (asciibuffer, sizeof(asciibuffer));
			
			while (fileStream.peek()!='\n' && fileStream.get(base)){
				base=toupper(base);
				fileStream.read(kmerbuffer, sizeof(kmerbuffer));
				string kmer (kmerbuffer, kmersize*2/3);
				
				unordered_map < string, char >::iterator it2 = kmerMap.find(kmer);//check if the kmer is in the hash
				if ( it2 == kmerMap.end() ){//if the kmer is not in the hash
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
				weededfile << endl;
			}
			firstkmer=true;
			fileStream.ignore(256,'\n');//skip the end ofline character
    	}
		fileStream.close();
		weededfile.close();
	}
	cout << totalweedkmers << " kmers prior to weeding" << endl;
	cout << weeded << " kmers remaining after weeding" << endl;
	cout << totalweedkmers-weeded << " kmers weeded (Note this may be more than the number of kmers in the weed file due to mulitple variants of a kmer being in the file)" << endl;

	printDuration(start);
	
	return 0;
	
	
}


