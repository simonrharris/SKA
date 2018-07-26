//g++ -O3 -std=c++0x splitkmer.cpp -lz -o splitkmer
#include <iostream>
#include <fstream>
#include <string>
#include <vector> //std::vector
#include <cmath>       /* ceil */
#include "kmers.hpp"
#include "general.hpp"
using namespace std;


//int main(int argc, char *argv[])
int summariseKmerFiles(const vector<string> & kmerfiles)
{
	int kmersize;
	vector < string > fileSampleNames;
	//cout << "Summary of bases in splitkmer files:\n";
	cout << "Sample\tKmer size\tTotal kmers\tAs\tCs\tGs\tTs\tNs\tOthers\tGC Content\n";
	for (int s = 0; s < kmerfiles.size(); ++s){
		ifstream fileStream;
		fileStream.open(kmerfiles[s], ios::in);
		if (fileStream.fail()) {
			cout << "Failed to open " << kmerfiles[s] << "\n\n";
			return 0;
		}
		
		try {
			int returnval = readKmerHeader(fileStream, kmersize, fileSampleNames);
		}
		catch (int e){
			cout << "An exception occurred when reading file " << kmerfiles[s] << ". Please check the format. Exception Nr. " << e << '\n';
			return 1;
		}

		vector < int > a(fileSampleNames.size(), 0);
		vector < int > c(fileSampleNames.size(), 0);
		vector < int > g(fileSampleNames.size(), 0);
		vector < int > t(fileSampleNames.size(), 0);
		vector < int > n(fileSampleNames.size(), 0);
		vector < int > other(fileSampleNames.size(), 0);
		vector < int > kmers(fileSampleNames.size(), 0);
		
		char basebuffer[1];
		char kmerbuffer[kmersize*2/3];
		char asciibuffer[int(ceil(float(fileSampleNames.size())/6))];
		while (fileStream.read(asciibuffer, sizeof(asciibuffer))){//read the ascii representation of the taxa
			string asciibits (asciibuffer, sizeof(asciibuffer));
			vector < bool > mybits;
			vectorbool_from_ascii(asciibits, mybits);//read the ascii representation of the taxa to a verctor of bools

			while (fileStream.peek()!='\n' && fileStream.read(basebuffer, sizeof(basebuffer))){
				string base (basebuffer, 1);
				base[0]=toupper(base[0]);
				fileStream.read(kmerbuffer, sizeof(kmerbuffer));
				string kmer (kmerbuffer, kmersize*2/3);
				
				for (int i=0; i<fileSampleNames.size(); ++i){ //add the kmer count to all samples that are true in the bitset
						if (mybits[i]==1){
							kmers[i]++;
						}
					}
				
				switch (base[0])
				{
					case 'A':
						for (int i=0; i<fileSampleNames.size(); ++i){ //add the base to all samples that are true in the bitset
							if (mybits[i]==1){
								a[i]++;
							}
						}
						break;
					case 'C':
						for (int i=0; i<fileSampleNames.size(); ++i){ //add the base to all samples that are true in the bitset
							if (mybits[i]==1){
								c[i]++;
							}
						}
						break;
					case 'G':
						for (int i=0; i<fileSampleNames.size(); ++i){ //add the base to all samples that are true in the bitset
							if (mybits[i]==1){
								g[i]++;
							}
						}
						break;
					case 'T':
						for (int i=0; i<fileSampleNames.size(); ++i){ //add the base to all samples that are true in the bitset
							if (mybits[i]==1){
								t[i]++;
							}
						}
						break;
					case 'N':
						for (int i=0; i<fileSampleNames.size(); ++i){ //add the base to all samples that are true in the bitset
							if (mybits[i]==1){
								n[i]++;
							}
						}
						break;
					default:
						for (int i=0; i<fileSampleNames.size(); ++i){ //add the base to all samples that are true in the bitset
							if (mybits[i]==1){
								other[i]++;
							}
						}
				}
			}
			fileStream.ignore(256,'\n');//skip the end ofline character
    	}
		fileStream.close();
		
		for (int i=0; i<fileSampleNames.size(); ++i){ //add the base to all samples that are true in the bitset		
			float gccontent=(float(c[i])+float(g[i]))/(float(a[i])+float(c[i])+float(g[i])+float(t[i]));
			cout << fileSampleNames[i] << "\t" << kmersize << "\t" << kmers[i] << "\t" << a[i] << "\t" << c[i] << "\t" << g[i] << "\t" << t[i] << "\t" << n[i] << "\t" << other[i] << "\t" << gccontent << "\n";
		}
    	fileSampleNames.clear();
	}

	return 0;
		
}


