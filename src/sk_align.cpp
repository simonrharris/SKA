//g++ -O3 -std=c++0x splitkmer.cpp -lz -o splitkmer
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "general.hpp"
#include "kmers.hpp"
#include <chrono> //timing
using namespace std;


//int main(int argc, char *argv[])
int alignKmers(const float & minproportion, const string & outputfile, const vector<string> & kmerfiles)
{

	auto start = chrono::high_resolution_clock::now();
	int numfiles=kmerfiles.size();
	
	float maxmissing=(1.0-minproportion)*numfiles;

	if (maxmissing<1){
		maxmissing=1;
	}

	float minrequired=numfiles-maxmissing;
	
	cout << "Keeping variants for which at least " << int(minrequired) << " samples include kmer matches" << "\n\n";
	
	// Create the kmer map
	unordered_map<string, string> kmerMap;
	string emptySequence (numfiles , '-');
	int oldkmersize=0;

	ifstream fileStream;
	for (int s = 0; s < kmerfiles.size(); ++s){
		cout << "Reading " << kmerfiles[s] << "\n";
		fileStream.open(kmerfiles[s], ios::in);

		if (fileStream.fail()) {
			cout << "Failed to open" << kmerfiles[s] << "\n\n";
			return 0;
		}
		int kmersize=readKmerHeader(fileStream);
		if (s==0){
			oldkmersize=kmersize;
		}

		if (kmersize!=oldkmersize){
			cout << "kmer files have different kmer sizes\n\n";
			return 0;
		}

		char basebuffer[1];
		char kmerbuffer[kmersize*2/3];

		while (fileStream.read(basebuffer, sizeof(basebuffer))){
			string base (basebuffer, 1);
			fileStream.read(kmerbuffer, sizeof(kmerbuffer));
			string kmer (kmerbuffer, kmersize*2/3);
			
			auto it = kmerMap.find(kmer);//check if the kmer is in the hash
			if ( it != kmerMap.end() ){//if the kmer is in the hash
				it->second[s]=base[0];//make the hash a repeat
			}
			else {
				if (s+1<=maxmissing){
					string newsequence = emptySequence;
					newsequence[s] = base[0];
					kmerMap.insert(make_pair(kmer, newsequence));
				}
			}
    	}
		fileStream.close();
	}
	cout << kmerMap.size() << " kmers read from " << numfiles << " files\n";
	
	int consta = 0;
	int constc = 0;
	int constg = 0;
	int constt = 0;
	bool variantonly=true;

	auto it = kmerMap.begin();
	auto endIter = kmerMap.end();

	for (; it!=endIter; ){
		int a = 0;
		int c = 0;
		int g = 0;
		int t = 0;
		int afound = 0;
		int cfound = 0;
		int gfound = 0;
		int tfound = 0;
		for (int i=0; i<numfiles; ++i){
			switch (it->second[i])
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
		if ((a+c+g+t)<minrequired){
			kmerMap.erase(it++);
		}
		else if ((afound+cfound+gfound+tfound)==1 && variantonly){
			if (afound>0){
				consta++;
			}
			if (cfound>0){
				constc++;
			}
			if (gfound>0){
				constg++;
			}
			if (tfound>0){
				constt++;
			}
			kmerMap.erase(it++);
		}
		else {
			++it;
		}
		
		
	}
	
	
	cout << "Printing alignment of "<< kmerMap.size() << " sites\n";
	
	ofstream alignfile(outputfile);
	
	for (int i=0; i<numfiles; ++i){
		string filename =splitFileName(kmerfiles[i]);
		alignfile << ">" << filename << "\n";
		float nonns=0.0;
		for (auto it=kmerMap.begin(); it!=kmerMap.end(); ++it){
			alignfile << it->second[i];
			if (it->second[i]=='A' || it->second[i]=='C' || it->second[i]=='G' || it->second[i]=='T'){
				nonns++;
			}
		}
		alignfile << "\n";
		if ((nonns/kmerMap.size())<0.5){
			cout << "Warning: " << filename << " only matches " << nonns/kmerMap.size()*100 << "% of kmers in the alignment\n";
		}
	}
	alignfile.close();

	cout << "Constant sites matching filters (a c g t):\n";
	cout << consta << " " << constc << " " << constg << " " << constt << "\n";

	auto finish = std::chrono::high_resolution_clock::now();
	chrono::duration<double> elapsed = finish - start;
	cout << "Done\n";
	cout << "Total time required: " << elapsed.count() << "s\n\n";
	
	return 0;
	
	
}


