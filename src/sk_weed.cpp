//g++ -O3 -std=c++0x src/sk_weed.cpp -lz -o bin/sk_weed
#include <unordered_map> //std::unordered_map
#include <iostream> //std::cout
#include <fstream> //std::ofstream
#include <string> //std::string
#include <vector> //std::vector
#include <chrono> //timing
#include <cmath>       /* ceil */
#include "kmers.hpp"

using namespace std;

//bin/sk_weed outputfile 

//int main(int argc, char *argv[])
int weedKmers(const vector<string> & weedfiles, const string & kmerfile)
{

	auto start = chrono::high_resolution_clock::now();
	// Create the kmer map
	unordered_map<string, string> kmerMap;
	
	cout << "Reading " << kmerfile << endl;
	ifstream fileStream;
	fileStream.open(kmerfile, ios::in);
	if (fileStream.fail()) {
		cout << "Failed to open " << kmerfile << endl << endl;
		return 0;
	}
	int kmersize;
	vector < string > weednames;
	try {
		int returnval = readKmerHeader(fileStream, kmersize, weednames);
	}
	catch (int e){
		cout << "An exception occurred when reading file " << kmerfile << ". Please check the format. Exception Nr. " << e << '\n';
		return 1;
	}
	char basebuffer[1];
	char kmerbuffer[kmersize*2/3];
	char asciibuffer[int(ceil(float(weednames.size())/6))];
	while (fileStream.read(asciibuffer, sizeof(asciibuffer))){
		string asciibits (asciibuffer, sizeof(asciibuffer));
		
		while (fileStream.peek()!='\n' && fileStream.read(basebuffer, sizeof(basebuffer))){
			string base (basebuffer, 1);
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
	
	for (auto it = weedfiles.begin(); it != weedfiles.end(); ++it){
		cout << "Weeding kmers from " << *it << endl;
		
		fileStream.open(*it, ios::in);
		if (fileStream.fail()) {
			cout << "Failed to open " << *it << endl << endl;
			return 0;
		}
		int weedkmersize;
		vector < string > names;
		
		try {
			int returnval = readKmerHeader(fileStream, weedkmersize, names);
		}
		catch (int e){
			cout << "An exception occurred when reading file " << *it << ". Please check the format. Exception Nr. " << e << '\n';
			return 1;
		}
		if (weedkmersize!=kmersize){
			cout << "Files have different kmer sizes" << endl << endl;
			return 0;
		}
		string filename=string(*it);
		string outputfile=filename.substr(0, filename.find_last_of("."))+".weeded"+filename.substr(filename.find_last_of("."), filename.length());

		cout << "Writing weeded kmers to " << outputfile << endl;
	
		ofstream weededfile(outputfile);
		weededfile << kmersize << "\n";
		for ( auto it=names.begin(); it!=names.end(); ++it){ //print each sample name to output file stream
			weededfile << *it << " "; 
		}
		weededfile << endl;

		bool firstkmer=true;

		char asciibuffer[int(ceil(float(names.size())/6))];
		while (fileStream.read(asciibuffer, sizeof(asciibuffer))){
			string asciibits (asciibuffer, sizeof(asciibuffer));
			
			while (fileStream.peek()!='\n' && fileStream.read(basebuffer, sizeof(basebuffer))){
				string base (basebuffer, 1);
				fileStream.read(kmerbuffer, sizeof(kmerbuffer));
				string kmer (kmerbuffer, kmersize*2/3);
				
				auto it2 = kmerMap.find(kmer);//check if the kmer is in the hash
				if ( it2 == kmerMap.end() ){//if the kmer is not in the hash
					if (firstkmer){//if it's the first kmer, print the ascii bitstring
						weededfile << asciibits;
						firstkmer=false;
					}
					weededfile << base[0] << kmer;//print the kmer
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

	auto finish = std::chrono::high_resolution_clock::now();
	chrono::duration<double> elapsed = finish - start;
	cout << "Done\n";
	cout << "Total time required: " << elapsed.count() << "s" << endl << endl;
	
	return 0;
	
	
}


