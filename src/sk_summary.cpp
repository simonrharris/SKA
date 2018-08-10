#include <iostream>
#include <fstream>
#include <string>
#include <vector> //std::vector
#include <map>
#include <cmath>       /* ceil */
#include "kmers.hpp"
#include "general.hpp"
using namespace std;


int summariseKmerFiles(const vector<string> & kmerfiles)
{
	int kmersize;
	vector < string > fileSampleNames;
	ifstream fileStream;
	
	cout << "Sample\tKmer size\tTotal kmers\tAs\tCs\tGs\tTs\tNs\tOthers\tGC Content\n";
	for (int s = 0; s < kmerfiles.size(); ++s){
		
		if (openFileStream(kmerfiles[s], fileStream, false)){return 1;};
		
		readKmerHeader(fileStream, kmersize, fileSampleNames);

		map < char, vector < int > > basecounts;

		vector < char > bases {'A', 'C', 'G', 'T', 'N', 'O'};

		for (vector<char>::iterator it=bases.begin(); it!=bases.end(); ++it){
			basecounts.insert(make_pair(*it, vector < int > (fileSampleNames.size(), 0)));
		}

		vector < int > kmers(fileSampleNames.size(), 0);
		
		char basebuffer[1];
		char kmerbuffer[kmersize*2/3];
		char asciibuffer[int(ceil(float(fileSampleNames.size())/6))];
		while (fileStream.read(asciibuffer, sizeof(asciibuffer))){//read the ascii representation of the taxa
			string asciibits (asciibuffer, sizeof(asciibuffer));
			vector < bool > mybits;
			vectorbool_from_ascii(asciibits, mybits);//read the ascii representation of the taxa to a vector of bools

			while (fileStream.peek()!='\n' && fileStream.read(basebuffer, sizeof(basebuffer))){
				string base (basebuffer, 1);
				base[0]=toupper(base[0]);
				fileStream.read(kmerbuffer, sizeof(kmerbuffer));
				string kmer (kmerbuffer, kmersize*2/3);
				
				for (int i=0; i<fileSampleNames.size(); ++i){ //add the kmer count to all samples that are true in the bitset
						if (mybits[i]==1){
							kmers[i]++;
							map < char, vector < int > >::iterator it = basecounts.find(base[0]);//check if the kmer is in the map
							if ( it != basecounts.end() ){//if the kmer is in the map
								it->second[i]++;//increment the count for the base
							}
							else {//if the kmer isn't in the map
								basecounts['O'][i]++;//increment the count for other
							}
						}
					}
			}
			fileStream.ignore(256,'\n');//skip the end ofline character
    	}
		fileStream.close();
		
		for (int i=0; i<fileSampleNames.size(); ++i){ //print the results	
			cout << fileSampleNames[i] << "\t" << kmersize << "\t" << kmers[i] << "\t";
			for (vector<char>::iterator it=bases.begin(); it!=bases.end(); ++it){
				cout << basecounts[*it][i] << "\t";
			}
			float gccontent=(float(basecounts['C'][i])+float(basecounts['G'][i]))/(float(basecounts['A'][i])+float(basecounts['C'][i])+float(basecounts['G'][i])+float(basecounts['T'][i]));
			cout << gccontent << endl;
		}
    	fileSampleNames.clear();
	}

	return 0;
		
}


