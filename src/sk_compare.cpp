//g++ -O3 -std=c++0x src/sk_compare -lz -o bin/sk_compare
#include <unordered_map> //std::unordered_map
#include <iostream> //std::cout
#include <fstream> //std::ofstream
#include <string> //std::string
#include <vector> //std::vector
#include <cmath>       /* ceil */
#include "kmers.hpp"
#include "general.hpp"

using namespace std;


//int main(int argc, char *argv[])
int compareKmerFiles(const string & queryfile, const vector<string> & subjectfiles)
{

	// Create the kmer map
	unordered_map<string, string> kmerMap;


	vector < string > sampleNames;
	if (collectSampleNames(subjectfiles, sampleNames)!=0){return 1;}
	int numSamples=sampleNames.size();
	string emptySequence (numSamples , '-');
	vector < int > queryKmers (numSamples, 0);

	int sampleNum=0;
	int subjectkmersize;
	int oldsubjectkmersize;
	vector < string > subjectSampleNames;

	ifstream fileStream;

	for (int s = 0; s < subjectfiles.size(); ++s){
		
		fileStream.open(subjectfiles[s], ios::in);

		if (fileStream.fail()) {
			cout << "Failed to open " << subjectfiles[s] << "\n\n";
			return 0;
		}
		
		
		try {
			int returnval = readKmerHeader(fileStream, subjectkmersize, subjectSampleNames);
		}
		catch (int e){
			cout << "An exception occurred when reading file " << subjectfiles[s] << ". Please check the format. Exception Nr. " << e << '\n';
			return 1;
		}

		if (s==0){
			oldsubjectkmersize=subjectkmersize;
		}

		if (subjectkmersize!=oldsubjectkmersize){
			cout << "kmer files have different kmer sizes\n\n";
			return 0;
		}

		char basebuffer[1];
		char kmerbuffer[subjectkmersize*2/3];
		char asciibuffer[int(ceil(float(subjectSampleNames.size())/6))];
		while (fileStream.read(asciibuffer, sizeof(asciibuffer))){//read the ascii representation of the taxa
			string asciibits (asciibuffer, sizeof(asciibuffer));
			vector < bool > mybits;
			vectorbool_from_ascii(asciibits, mybits);//read the ascii representation of the taxa to a verctor of bools
			while (fileStream.peek()!='\n' && fileStream.read(basebuffer, sizeof(basebuffer))){
				string base (basebuffer, 1);
				base[0]=toupper(base[0]);
				fileStream.read(kmerbuffer, sizeof(kmerbuffer));
				string kmer (kmerbuffer, subjectkmersize*2/3);
				
				auto it = kmerMap.find(kmer);//check if the kmer is in the hash
				if ( it != kmerMap.end() ){//if the kmer is in the hash
					for (int i=0; i<subjectSampleNames.size(); ++i){ //add the base to all samples that are true in the bitset
						if (mybits[i]){
							it->second[i+sampleNum]=base[0];
							queryKmers[i+sampleNum]++;
						}
					}
				}
				else {
					string newsequence = emptySequence;
					for (int i=0; i<subjectSampleNames.size(); ++i){ //add the base to all samples that are true in the bitset
						if (mybits[i]){
							newsequence[i+sampleNum]=base[0];
							queryKmers[i+sampleNum]++;
						}
					}
					kmerMap.insert(make_pair(kmer, newsequence));
				}
			}
			fileStream.ignore(256,'\n');//skip the end ofline character
		}
		sampleNum+=int(subjectSampleNames.size()); //add the number of samples in the file to the count of total samples
		subjectSampleNames.clear();
		fileStream.close();
	}
	
	
	
	
	fileStream.open(queryfile, ios::in);
	if (fileStream.fail()) {
		cout << "Failed to open " << queryfile << "\n\n";
		return 1;
	}
	int querykmersize;
	vector < string > querySampleNames;
	if (querySampleNames.size()>1){
		cout << "Error: Query can only be a single sample, not a kmerge" << endl;
		return 1;
	}
	
	try {
			int returnval = readKmerHeader(fileStream, querykmersize, querySampleNames);
		}
		catch (int e){
			cout << "An exception occurred when reading file " << queryfile << ". Please check the format. Exception Nr. " << e << '\n';
			return 1;
		}
	if (subjectkmersize!=querykmersize){
		cout << queryfile << " has different kmer size to query\n";
		return 1;
	}

	cout << "Subject\tKmers unique to Query\tKmers unique to Subject\tMatches\t% kmers in Query matching\t% kmers in Subject matching\tSNPs\t%ID of matching kmers\t%ID of Query kmers\t%ID of Subject kmers\tNs in Query\tNs in Subject\tNs in both\n";

	
	vector < int > kmerjustinb (sampleNames.size(), 0);
	vector < int > snps (sampleNames.size(), 0);
	vector < int > nina (sampleNames.size(), 0);
	vector < int > ninb (sampleNames.size(), 0);
	vector < int > ninboth (sampleNames.size(), 0);
	vector < int > matches (sampleNames.size(), 0);

	char basebuffer[1];
	char kmerbuffer[querykmersize*2/3];
	char asciibuffer[int(ceil(float(querySampleNames.size())/6))];

	while (fileStream.read(asciibuffer, sizeof(asciibuffer))){//read the ascii representation of the taxa
		string asciibits (asciibuffer, sizeof(asciibuffer));
		//vector < bool > mybits;
		//vectorbool_from_ascii(asciibits, mybits);//read the ascii representation of the taxa to a verctor of bools
		while (fileStream.peek()!='\n' && fileStream.read(basebuffer, sizeof(basebuffer))){
			string base (basebuffer, 1);
			base[0]=toupper(base[0]);
			fileStream.read(kmerbuffer, sizeof(kmerbuffer));
			string kmer (kmerbuffer, querykmersize*2/3);
			
			//cout << kmer << endl;

			auto it = kmerMap.find(kmer);//check if the kmer is in the hash
			if ( it != kmerMap.end() ){//if the kmer is in the hash
				for (int i=0; i<sampleNames.size(); ++i){ //add the base to all samples that are true in the bitset
					//cout << it->second[i];
					if (it->second[i]=='N' && base[0]=='N'){
						ninboth[i]++;
					}
					else if (it->second[i]=='-'){
						kmerjustinb[i]++;
					}
					else if (it->second[i]=='N'){
						nina[i]++;
					}
					else if (base[0]=='N'){
						ninb[i]++;
					}
					else if (it->second[i]==base[0]){
						matches[i]++;
					}
					else {
						snps[i]++;
					}
				}
				//cout << endl;
			}
			else {
				for (int i=0; i<sampleNames.size(); ++i){ //add the base to all samples that are true in the bitset
					kmerjustinb[i]++;
				}
			}
    	}
    	fileStream.ignore(256,'\n');//skip the end ofline character
    }
	fileStream.close();
	
	for (int i=0; i<sampleNames.size(); ++i){

		int kmerjustina=queryKmers[i]-(matches[i]+snps[i]+nina[i]+ninboth[i]+ninb[i]);
		float percentmatcha=float(matches[i]+snps[i]+nina[i]+ninboth[i])/(kmerjustina+matches[i]+snps[i]+nina[i]+ninboth[i])*100;
		float percentmatchb=float(matches[i]+snps[i]+ninb[i]+ninboth[i])/(kmerjustinb[i]+matches[i]+snps[i]+ninb[i]+ninboth[i])*100;
		float percentidofmatches=float(matches[i])/(matches[i]+snps[i])*100;
		float percentidofquery=(percentidofmatches*percentmatcha)/100;
		float percentidofsubject=(percentidofmatches*percentmatchb)/100;
		
		cout << sampleNames[i] << "\t" << kmerjustina << "\t" << kmerjustinb[i] << "\t" << matches[i] << "\t" << percentmatcha << "\t" << percentmatchb << "\t" << snps[i] << "\t" << percentidofmatches << "\t" << percentidofquery << "\t" << percentidofsubject << "\t" << nina[i] << "\t" << ninb[i] << "\t" << ninboth[i] << "\n";

	}
	return 0;
	
}


