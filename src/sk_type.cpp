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
int typeKmerFiles(const string & queryfile, const vector<string> & subjectfiles)
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
	int uniquecount=0;
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
			if (islower(base[0])){
				uniquecount++;
			}
			kmerMap.insert(make_pair(kmer, base));
		}
			
	}
	fileStream.close();
	
	//cout << "Subject\t# unique typing kmers\t% unique typing kmers matching\tunique kmer SNPs\t# non-unique typing kmers\t% non-unique typing kmers matching\tnon-unique kmer SNPs\t# sample kmers\t% samples kmers matching\t% sample kmers matching\tConfidence\n";
	cout << "Subject\t% unique typing kmers matching\t% non-unique typing kmers matching\t% samples kmers matching\tConfidence\n";

	for (int i=0; i<subjectfiles.size(); ++i){
		int kmerjustinb=0;
		int snps=0;
		int n=0;
		int un=0;
		int matches=0;
		int umatches=0;
		int usnps=0;
		int samplekmers=0;
		
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
			samplekmers++;
			string base (basebuffer, 1);
			base[0]=toupper(base[0]);
			fileStream.read(kmerbuffer, sizeof(kmerbuffer));
			string kmer (kmerbuffer, subjectkmersize*2/3);
			
			auto it = kmerMap.find(kmer);//check if the kmer is in the hash
			if ( it != kmerMap.end() ){//if the kmer is in the hash
				string typebase=it->second;
				bool isunique=islower(typebase[0]);
				typebase[0]=toupper(typebase[0]);
				if (typebase[0]=='N' || base[0]=='N'){
					if (isunique){
						un++;
					}
					else{
						n++;
					}
				}
				else if (typebase[0]==base[0]){
					if (isunique){
						umatches++;
					}
					else{
						matches++;
					}
				}
				else {
					if (isunique){
						usnps++;
					}
					else {
						snps++;
					}
				}
			}
			else {
				kmerjustinb++;
			}
    	}
	fileStream.close();
	
	//int kmerjustina=kmerMap.size()-(matches+snps+umatches+usnps+n+un);

	float percentsamplematch=float(matches+snps+umatches+usnps+n+un)/samplekmers*100;
	float percentuniquematches=float(umatches+un)/uniquecount*100;
	float percentnonuniquematch=float(matches+n)/(kmerMap.size()-(umatches+usnps+n))*100;
	float uniquetononuniqueratio=(percentuniquematches/percentnonuniquematch);
	if (uniquetononuniqueratio>1){
		uniquetononuniqueratio=1;
	}

	string filename=splitFileName(subjectfiles[i]);
	//cout << filename << "\t" << uniquecount << "\t" << percentuniquematches << "\t" << usnps << "\t" << kmerMap.size() << "\t" << percentnonuniquematch << "\t" << snps << "\t" << samplekmers << "\t" << percentsamplematch << "\t" << uniquetononuniqueratio*percentsamplematch << "\n";
	
	cout << filename << "\t" << percentuniquematches << "\t" << percentnonuniquematch << "\t" << percentsamplematch << "\t" << uniquetononuniqueratio*percentsamplematch << "\n";
	}

	return 0;
	
}


