//g++ -O3 -std=c++0x splitkmer.cpp -lz -o splitkmer
#include <iostream>
#include <fstream>
#include <string>
#include <vector> //std::vector
#include "kmers.hpp"
#include "general.hpp"
using namespace std;


//int main(int argc, char *argv[])
int summariseKmerFiles(const vector<string> & kmerfiles)
{

	int numfiles=0;
	//cout << "Summary of bases in splitkmer files:\n";
	cout << "File\tKmer size\tTotal kmers\tAs\tCs\tGs\tTs\tNs\tOthers\tGC Content\n";
	for (int s = 0; s < kmerfiles.size(); ++s){
		++numfiles;
		
		int a = 0;
		int c = 0;
		int g = 0;
		int t = 0;
		int n = 0;
		int other = 0;
		int kmers=0;
		
		ifstream fileStream;
		fileStream.open(kmerfiles[s], ios::in);
		if (fileStream.fail()) {
			cout << "Failed to open " << kmerfiles[s] << "\n\n";
			continue;
		}
		int kmersize=readKmerHeader(fileStream);
		
		char basebuffer[1];
		char kmerbuffer[kmersize*2/3];
		while (fileStream.read(basebuffer, sizeof(basebuffer))){
			string base (basebuffer, 1);
			base[0]=toupper(base[0]);
			fileStream.read(kmerbuffer, sizeof(kmerbuffer));
			string kmer (kmerbuffer, kmersize*2/3);
			
			kmers+=1;
			
			switch (base[0])
			{
				case 'A':
					a++;
					break;
				case 'C':
					c++;
					break;
				case 'G':
					g++;
					break;
				case 'T':
					t++;
					break;
				case 'N':
					n++;
					break;
				default:
					other++;
			}	
	
    		}
		fileStream.close();
	
		string filename=splitFileName(kmerfiles[s]);
		
		float gccontent=(float(c)+float(g))/(float(a)+float(c)+float(g)+float(t));
		cout << filename << "\t" << kmersize << "\t" << kmers << "\t" << a << "\t" << c  << "\t" << g  << "\t" << t << "\t" << n << "\t" << other << "\t" << gccontent << "\n";
	}

	return 0;
		
}


