//g++ -O3 -std=c++0x src/sk_compare -lz -o bin/sk_compare
#include <unordered_map> //std::unordered_map
#include <iostream> //std::cout
#include <fstream> //std::ofstream
#include <string> //std::string
#include <vector> //std::vector
#include "kmers.hpp"
#include "general.hpp"

using namespace std;


//int main(int argc, char *argv[])
int compareKmerFiles(const string & queryfile, const vector<string> & subjectfiles)
{

	// Create the kmer map
	unordered_map<string, string> kmerMap;
	int numfiles=0;
	
	ifstream fileStream;
	fileStream.open(queryfile, ios::in);

	if (fileStream.fail()) {
		cout << "Failed to open " << queryfile << "\n\n";
		return 0;
	}

	int querykmersize=readKmerHeader(fileStream);

	char basebuffer[1];
	char kmerbuffer[querykmersize*2/3];
	while (fileStream.read(basebuffer, sizeof(basebuffer))){
		string base (basebuffer, 1);
		fileStream.read(kmerbuffer, sizeof(kmerbuffer));
		string kmer (kmerbuffer, querykmersize*2/3);
		
		auto it = kmerMap.find(kmer);//check if the kmer is in the hash
		if ( it != kmerMap.end() ){//if the kmer is in the hash
			it->second=base;//make the hash a repeat
			cout << "shouldn't be here\n";
		}
		else {
			kmerMap.insert(make_pair(kmer, base));
		}
			
	}
	fileStream.close();
	
	cout << "Subject\tKmers unique to Query\tKmers unique to Subject\tMatches\t% kmers in Query matching\t% kmers in Subject matching\tSNPs\t%ID of matching kmers\t%ID of Query kmers\t%ID of Subject kmers\tNs in Query\tNs in Subject\tNs in both\n";

	for (int i=0; i<subjectfiles.size(); ++i){
		int kmerjustinb=0;
		int snps=0;
		int nina=0;
		int ninb=0;
		int ninboth=0;
		int matches=0;
		
		fileStream.open(subjectfiles[i], ios::in);
		if (fileStream.fail()) {
			cout << "Failed to open " << subjectfiles[i] << "\n\n";
			continue;
		}
		int subjectkmersize=readKmerHeader(fileStream);
		if (subjectkmersize!=querykmersize){
			cout << subjectfiles[i] << " has different kmer size to query\n";
			return 0;
		}
		while (fileStream.read(basebuffer, sizeof(basebuffer))){
			string base (basebuffer, 1);
			fileStream.read(kmerbuffer, sizeof(kmerbuffer));
			string kmer (kmerbuffer, subjectkmersize*2/3);
			
			auto it = kmerMap.find(kmer);//check if the kmer is in the hash
			if ( it != kmerMap.end() ){//if the kmer is in the hash
				if (it->second[0]=='N' && base[0]=='N'){
					ninboth++;
				}
				else if (it->second[0]=='N'){
					nina++;
				}
				else if (base[0]=='N'){
					ninb++;
				}
				else if (it->second[0]==base[0]){
					matches++;
				}
				else {
					snps++;
				}
			}
			else {
				kmerjustinb++;
			}
    	}
	fileStream.close();
	
	int kmerjustina=kmerMap.size()-matches;
	float percentmatcha=float(matches+snps+nina+ninboth)/(kmerjustina+matches+snps+nina+ninboth)*100;
	float percentmatchb=float(matches+snps+ninb+ninboth)/(kmerjustinb+matches+snps+ninb+ninboth)*100;
	float percentidofmatches=float(matches)/(matches+snps)*100;
	float percentidofquery=(percentidofmatches*percentmatcha)/100;
	float percentidofsubject=(percentidofmatches*percentmatchb)/100;

	string filename=splitFileName(subjectfiles[i]);
	
	cout << filename << "\t" << kmerjustina << "\t" << kmerjustinb << "\t" << matches << "\t" << percentmatcha << "\t" << percentmatchb << "\t" << snps << "\t" << percentidofmatches << "\t" << percentidofquery << "\t" << percentidofsubject << "\t" << nina << "\t" << ninb << "\t" << ninboth << "\n";
	}

	return 0;
	
}


