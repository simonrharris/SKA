//g++ -O3 -std=c++0x src/sk_compare -lz -o bin/sk_type
#include <unordered_map> //std::unordered_map
#include <map> //std::unordered_map
#include <iostream> //std::cout
#include <fstream> //std::ofstream
#include <string> //std::string
#include <vector> //std::vector
#include <cmath>       /* ceil */
#include <sstream> //std::stringstream
#include "kmers.hpp"
#include "general.hpp"
#include "gzstream.h"
#include "float.h"
#include <set>
#include <algorithm>

using namespace std;


bool AlmostEqualRelative(float A, float B,
                         float maxRelDiff = FLT_EPSILON)
{
    // Calculate the difference.
    float diff = fabs(A - B);
    A = fabs(A);
    B = fabs(B);
    // Find the largest
    float largest = (B > A) ? B : A;
 
    if (diff <= largest * maxRelDiff){
        return true;
    }
    else { 
    	return false;
    }
}

//int main(int argc, char *argv[])
int typeKmerFile(const string & queryfile, const string & profileFile, const vector<string> & subjectfiles)
{
	
	// Create the kmer map
	unordered_map < string, string > kmerMap;
	
	ifstream fileStream;
	//cout << "Reading " << queryfile << endl << endl;
	fileStream.open(queryfile, ios::in);
	if (fileStream.fail()) {
		cout << "Failed to open " << queryfile << "\n\n";
		return 1;
	}
	int querykmersize;
	vector < string > querySampleNames;
	
	try {
		int returnval = readKmerHeader(fileStream, querykmersize, querySampleNames);
	}
	catch (int e){
		cout << "An exception occurred when reading file " << queryfile << ". Please check the format. Exception Nr. " << e << '\n';
		return 1;
	}

	char basebuffer[1];
	char kmerbuffer[querykmersize*2/3];
	char asciibuffer[int(ceil(float(querySampleNames.size())/6))];
	vector < unordered_map < string, vector < int > > > profileMap (querySampleNames.size());
	vector < set < string > > novel ( querySampleNames.size());
	vector < set < string > > hasNs ( querySampleNames.size());
	vector < set < string > > hasgaps ( querySampleNames.size());
	vector < bool > uncertain ( querySampleNames.size(), false);

	while (fileStream.read(asciibuffer, sizeof(asciibuffer))){//read the ascii representation of the taxa
		string asciibits (asciibuffer, sizeof(asciibuffer));
		vector < bool > mybits;
		vectorbool_from_ascii(asciibits, mybits);//read the ascii representation of the taxa to a verctor of bools
		while (fileStream.peek()!='\n' && fileStream.read(basebuffer, sizeof(basebuffer))){
			string base (basebuffer, 1);
			base[0]=toupper(base[0]);
			fileStream.read(kmerbuffer, sizeof(kmerbuffer));
			string kmer (kmerbuffer, querykmersize*2/3);

			auto it = kmerMap.find(kmer);//check if the kmer is in the hash
			if ( it != kmerMap.end() ){//if the kmer is in the hash
				for (int i=0; i<mybits.size(); ++i){
					if (mybits[i]){
						it->second[i]=base[0];
					}
				}
			}
			else {
				auto ret = kmerMap.insert(make_pair(kmer, string (querySampleNames.size(),'-')));
				for (int i=0; i<mybits.size(); ++i){
					if (mybits[i]){
						ret.first->second[i]=base[0];
					}
				}
			}
    	}
    	fileStream.ignore(256,'\n');//skip the end ofline character
    }
	fileStream.close();


	//cout << "Best matching alleles:"<< endl;

	string header;
	string sequence;
	int substringlength=(querykmersize*2)+1;

	for (int s = 0; s < subjectfiles.size(); ++s){

		int numAlleles=0;
		vector <string> alleleNames;
		vector < float > bestid (querySampleNames.size(), 0);
		vector < int > Ncount (querySampleNames.size(), 0);
		vector < int > gapcount (querySampleNames.size(), 0);
		vector < vector < int > > bestMatches (querySampleNames.size());
		vector < string > bestSequences (querySampleNames.size());

		igzstream gzfileStream;
		gzfileStream.open(subjectfiles[s].c_str());

		if (gzfileStream.get()!='>'){
			cout << endl << "Error: " << subjectfiles[s] << " is not in the correct format. Expecting header to start with >." << endl << endl;
			return 1;
		}

		while (getline(gzfileStream, header)){
			getline(gzfileStream, sequence, '>');
			numAlleles++;
		}
 		gzfileStream.close();//annoyingly have to close and reopen the file as seek/rewind aren't implemented for gzstream

		int alleleNumber=0;

		gzfileStream.clear();
		gzfileStream.open(subjectfiles[s].c_str());

		if (gzfileStream.get()!='>'){
			cout << endl << "Error: " << subjectfiles[s] << " is not in the correct format. Expecting header to start with >." << endl << endl;
			continue;
		}

		while (getline(gzfileStream, header)){

			stringstream headerstream;
			headerstream << header;//convert the sequence to stringstream

			string alleleName;
			getline(headerstream, alleleName, ' ');
			alleleNames.push_back(alleleName);
		
			getline(gzfileStream, sequence, '>');

			sequence.erase(std::remove_if( sequence.begin(), sequence.end(), ::isspace ), sequence.end() );

			vector < string> alleleSequence ( querySampleNames.size() , string ( sequence.length() , '-' ) );
			vector < vector < bool > > alleleCoverage ( querySampleNames.size() , vector < bool > ( sequence.length() , false ) );
				
			if (sequence.length()<substringlength){
				continue;
			}
			
			transform(sequence.begin(), sequence.end(), sequence.begin(), ::toupper);//change all letters in the string to upper case
			
			int i=0;

			for (auto iti = sequence.cbegin(), end = sequence.cend()-(substringlength-1); iti != end; ++iti, ++i){
				
				string kmer=sequence.substr(i,substringlength);
				
				bool isrev=false;
				if (reverse_is_min(kmer, querykmersize+1)){
					reverse(kmer.begin(), kmer.end());
					transform(kmer.begin(),kmer.end(),kmer.begin(),complement);
					isrev=true;
				}
				
				char base=kmer[querykmersize];
				kmer.erase(kmer.begin()+querykmersize);
				ascii_codons(kmer);
				
				auto it = kmerMap.find(kmer);//check if the kmer is in the hash
				
				if ( it != kmerMap.end() ){//if the kmer is in the hash
					for (int k=0; k<querySampleNames.size(); ++k){
						if (it->second[k]!='-'){
							for (int j=i; j<i+substringlength; ++j){
								alleleCoverage[k][j]=true;

								if (j==(i+querykmersize)){
									char querybase;
									if (isrev){
										querybase=complement(it->second[k]);
									}
									else {
										querybase=it->second[k];
									}
									if (alleleSequence[k][j]=='-' or islower(alleleSequence[k][j])){
										alleleSequence[k][j]=querybase;
									}
									else if (querybase==base){
										//cout << "HERE" << endl;
										alleleSequence[k][j]=base;
									}
									else if (alleleSequence[k][j]!=base){
										//cout << querybase << " " << alleleSequence[k][j] << " " << base << endl;
										alleleSequence[k][j]='N';
									}
								}
								else if (alleleSequence[k][j]=='-'){
									alleleSequence[k][j]=tolower(sequence[j]);
								}
							}
						}
					}
				}
		
    		}

    		for (int k=0; k<querySampleNames.size(); ++k){
    			float matches=0.0;
    			float Ns=0.0;
    			float gaps=0.0;
    			char base;
    			for (int i=0; i<sequence.length(); ++i){
    				base=toupper(alleleSequence[k][i]);
    				if (base==sequence[i]){
    					matches++;
    				}
    				else if (base=='N'){
    					Ns++;
    				}
    				else if (base=='-'){
    					gaps++;
    				}
    			}
    			float sampleid=(Ns+matches)/sequence.length();
    			//cout  << alleleNumber << " " << sequence.length() << " " << matches << " " << sampleid << " " << Ns << " " << gaps << endl;
				if (sampleid>bestid[k]){
					bestMatches[k].clear();
					bestid[k]=sampleid;
					Ncount[k]=Ns;
					gapcount[k]=gaps;
				}
				//if (sampleid>0 && AlmostEqualRelative(sampleid, bestid[k])){
				if (sampleid>0 && sampleid==bestid[k]){
					bestMatches[k].push_back(alleleNumber);
					bestSequences[k]=alleleSequence[k];
				}
			}
			
    		alleleNumber++;
		}
 		gzfileStream.close();


 		for (int k=0; k<querySampleNames.size(); ++k){

 			int allele;
			string locus;
			//cout << querySampleNames[k];
	 		for (auto it=bestMatches[k].begin(); it!=bestMatches[k].end(); ++it){
	 			//cout << "\t" << alleleNames[*it];
				string::size_type found=alleleNames[*it].find_last_of("_");
				if (found!=string::npos){
					locus=alleleNames[*it].substr(0,found);
					allele=stoi(alleleNames[*it].substr(found+1));
				}
				auto it2 = profileMap[k].find(locus);//check if the kmer is in the hash
				if ( it2 != profileMap[k].end() ){//if the kmer is not in the hash
					it2->second.push_back(allele);
					uncertain[k]=true;
				}
				else{
					vector < int > newLocus;
					newLocus.push_back(allele);
					profileMap[k].insert(make_pair(locus, newLocus));
				}
			}

			if ((bestid[k]!=1.0 || Ncount[k]>0) && bestid[k]>0){//if the allele sequence of the sample has Ns or isn't a perfect match to an allele in the file, print the allele sequence to a file
				if (bestid[k]!=1.0){
					novel[k].insert(locus);
				}
				if (Ncount[k]>0){
					hasNs[k].insert(locus);
				}
				if (gapcount[k]>0){
					hasgaps[k].insert(locus);
				}
				ofstream alignfile(querySampleNames[k]+"_"+locus+".fa");
				alignfile << ">" << querySampleNames[k] << "_" << locus << endl;
				alignfile << bestSequences[k] << endl;
				alignfile.close();
			}
			/*cout << ">" << querySampleNames[k] << "_" << locus << endl;
			cout << bestSequences[k] << endl;*/
			//cout << "\t" << bestid[k] << "\t" << Ncount[k] << "\t" << gapcount[k] << endl;
		}
	}


	if (profileFile!=""){

		vector < int > certainFiles;
		for (int k=0; k<querySampleNames.size(); ++k){
			if (novel[k].size()==0 && not uncertain[k]){
				certainFiles.push_back(k);
			}
		}


		fileStream.clear();
		fileStream.open(profileFile, ios::in);

		string header;
		getline(fileStream, header);

		stringstream headerstream;
		headerstream << header;//convert the header to a stringstream
		string word;
		int i=0;
		int ST=-1;
		int cc=-1;
		map < string, int > alleles;
		int maxallelecolumn=0;

		while (getline(headerstream, word, '\t')){
			if (word=="ST"){
				ST=i;
				maxallelecolumn=i;
			}
			else if (word=="clonal_complex"){
				cc=i;
			}
			else {
				alleles[word]=i;
				maxallelecolumn=i;
			}
			++i;
		}

		if (alleles.size()!=profileMap[0].size()){
			cout << "the alleles in your profile file don't match those in your input files" << endl << endl;
			return 1;
		}

		for (auto it=alleles.begin(); it!=alleles.end(); ++it){
			auto it2 = profileMap[0].find(it->first);//check if the kmer is in the hash
			if ( it2 == profileMap[0].end() ){//if the kmer is not in the hash
				cout << "the alleles in your profile file don't match those in your input files" << endl << endl;
				return 1;
			}
		}

		cout << "Sample" << "\t";
		for (auto it=alleles.begin(); it!=alleles.end(); ++it){
			cout << it->first << "\t";
		}
		cout << "ST" << endl;

		string line;
		int linenumber=0;
		int matches=0;
		while (getline(fileStream, line)){
			linenumber++;
			stringstream linestream;
			linestream << line;//convert the header to a stringstream
			vector < string > words;

			while (getline(linestream, word, '\t')){
				words.push_back(word);
			}

			if (words.size()<maxallelecolumn){
				cout << "Malformed profile file on line " << linenumber << ". Expecting at least " << maxallelecolumn << " columns and found " << words.size() << endl << endl;
				return 1;
			}

			for (auto k=certainFiles.begin(); k!=certainFiles.end(); ++k){
				bool perfectmatch=true;
				for (auto it=alleles.begin(); it!=alleles.end(); ++it){
					auto it2 = profileMap[*k].find(it->first);//check if the kmer is in the hash
					if ( it2 == profileMap[*k].end() ){//if the kmer is not in the hash
						perfectmatch=false;
						cout << "Shouldn't be here";
					}
					else{
						if (stoi(words[it->second])!=it2->second[0]){
							perfectmatch=false;
							continue;
						}
					}
				}
				if (perfectmatch){
					vector < int > STs;
					STs.push_back(stoi(words[ST].substr(words[ST].find_last_of("_")+1)));
					profileMap[*k]["ST"]=STs;
				}
			}
			if (matches==certainFiles.size()){
				break;
			}
		}

		fileStream.close();
		

		for (int k=0; k<querySampleNames.size(); ++k){
			cout << querySampleNames[k] << "\t";
			for (auto it=alleles.begin(); it!=alleles.end(); ++it){
				string suffix="";
				auto it2 = novel[k].find(it->first);//check if the kmer is in the set
				if ( it2 != novel[k].end() ){
					suffix="*";
				}
				it2 = hasNs[k].find(it->first);//check if the kmer is in the set
				if ( it2 != hasNs[k].end() ){
					suffix+="N";
				}
				it2 = hasgaps[k].find(it->first);//check if the kmer is in the set
				if ( it2 != hasgaps[k].end() ){
					suffix+="-";
				}

				if (profileMap[k][it->first].size()==0){
					cout << "-" << "\t";
				}
				else if (profileMap[k][it->first].size()==1){
					cout << profileMap[k][it->first][0] << suffix << "\t";
				}
				else {
					for (int i=0; i<profileMap[k][it->first].size()-1; ++i){
						cout << profileMap[k][it->first][i] << "/";
					}
					cout << profileMap[k][it->first][profileMap[k][it->first].size()-1] << suffix << "\t";
				}
			}
			if (profileMap[k]["ST"].size()==1){
				cout << profileMap[k]["ST"][0] << endl;
			}
			else {
				cout << "-" << endl;
			}
		}

	}
	else{
		set < string > allLoci;
		for (int k=0; k<querySampleNames.size(); ++k){
			for (auto it=profileMap[k].begin(); it!=profileMap[k].end(); ++it){
				allLoci.insert(it->first);
				//cout << it->first << endl;
			}
		}
		cout << "Sample";
		for (auto it=allLoci.begin(); it!=allLoci.end(); ++it){
			cout << "\t" << *it;
		}
		cout << endl;
		for (int k=0; k<querySampleNames.size(); ++k){
			cout << querySampleNames[k] << "\t";
			for (auto it=allLoci.begin(); it!=allLoci.end(); ++it){
				string suffix="";
				auto it2 = novel[k].find(*it);//check if the kmer is in the set
				if ( it2 != novel[k].end() ){
					suffix="*";
				}
				it2 = hasNs[k].find(*it);//check if the kmer is in the set
				if ( it2 != hasNs[k].end() ){
					suffix+="N";
				}
				it2 = hasgaps[k].find(*it);//check if the kmer is in the set
				if ( it2 != hasgaps[k].end() ){
					suffix+="-";
				}

				if (profileMap[k][*it].size()==0){
					cout << "-" << "\t";
				}
				else if (profileMap[k][*it].size()==1){
					cout << profileMap[k][*it][0] << suffix << "\t";
				}
				else {
					for (int i=0; i<profileMap[k][*it].size()-1; ++i){
						cout << profileMap[k][*it][i] << "/";
					}
					cout << profileMap[k][*it][profileMap[k][*it].size()-1] << suffix << "\t";
				}
			}
			cout << endl;
		}
	}
	

	return 0;
	
}


