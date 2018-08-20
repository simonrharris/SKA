#include <unordered_map> //std::unordered_map
#include <map> //std::unordered_map
#include <iostream> //std::cout std::cerr
#include <fstream> //std::ofstream
#include <string> //std::string
#include <vector> //std::vector
#include <cmath> // std::ceil
#include <sstream> //std::stringstream
#include <set> //std::set
#include <algorithm> //std::transform
#include "kmers.hpp"
#include "general.hpp"
#include "DNA.hpp"
#include "gzstream.h"
#include "float.h"


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
int typeKmerFile(const std::string & queryfile, const std::string & profileFile, const std::vector < std::string > & subjectfiles)
{
	
	// Create the kmer map
	std::unordered_map < std::string, std::string > kmerMap;
	
	std::ifstream fileStream;

	if (openFileStream(queryfile, fileStream, false)){return 1;};

	int querykmersize;
	std::vector < std::string > querySampleNames;
	
	readKmerHeader(fileStream, querykmersize, querySampleNames);
	
	char base;
	char kmerbuffer[querykmersize*2/3];
	char asciibuffer[int(std::ceil(float(querySampleNames.size())/6))];
	std::vector < std::unordered_map < std::string, std::vector < int > > > profileMap (querySampleNames.size());
	std::vector < std::set < std::string > > novel ( querySampleNames.size());
	std::vector < std::set < std::string > > hasNs ( querySampleNames.size());
	std::vector < std::set < std::string > > hasgaps ( querySampleNames.size());
	std::vector < bool > uncertain ( querySampleNames.size(), false);

	while (fileStream.read(asciibuffer, sizeof(asciibuffer))){//read the ascii representation of the taxa
		std::string asciibits (asciibuffer, sizeof(asciibuffer));
		std::vector < bool > mybits;
		vectorbool_from_ascii(asciibits, mybits);//read the ascii representation of the taxa to a verctor of bools
		while (fileStream.peek()!='\n' && fileStream.get(base)){
			base=std::toupper(base);
			fileStream.read(kmerbuffer, sizeof(kmerbuffer));
			std::string kmer (kmerbuffer, querykmersize*2/3);

			std::unordered_map < std::string, std::string >::iterator it = kmerMap.find(kmer);//check if the kmer is in the hash
			if ( it != kmerMap.end() ){//if the kmer is in the hash
				for (int i=0; i<mybits.size(); ++i){
					if (mybits[i]){
						it->second[i]=base;
					}
				}
			}
			else {
				std::pair < std::unordered_map < std::string, std::string >::iterator, bool > ret = kmerMap.insert(std::make_pair(kmer, std::string (querySampleNames.size(),'-')));
				for (int i=0; i<mybits.size(); ++i){
					if (mybits[i]){
						ret.first->second[i]=base;
					}
				}
			}
    	}
    	fileStream.ignore(256,'\n');//skip the end ofline character
    }
	fileStream.close();


	//cout << "Best matching alleles:"<< endl;

	std::string header;
	std::string sequence;
	int substringlength=(querykmersize*2)+1;

	for (int s = 0; s < subjectfiles.size(); ++s){

		int numAlleles=0;
		std::vector < std::string > alleleNames;
		std::vector < float > bestid (querySampleNames.size(), 0);
		std::vector < int > Ncount (querySampleNames.size(), 0);
		std::vector < int > gapcount (querySampleNames.size(), 0);
		std::vector < std::vector < int > > bestMatches (querySampleNames.size());
		std::vector < std::string > bestSequences (querySampleNames.size());

		igzstream gzfileStream;
		gzfileStream.open(subjectfiles[s].c_str());
		if (gzfileStream.fail()) {
			cout << "Failed to open " << subjectfiles[s] << endl << endl;
			return 1;
		}

		if (gzfileStream.get()!='>'){
			std::cerr << endl << "Error: " << subjectfiles[s] << " is not in the correct format. Expecting header to start with >." << endl << endl;
			return 1;
		}

		while (std::getline(gzfileStream, header)){
			std::getline(gzfileStream, sequence, '>');
			numAlleles++;
		}
 		gzfileStream.close();//annoyingly have to close and reopen the file as seek/rewind aren't implemented for gzstream

		int alleleNumber=0;

		gzfileStream.clear();
		gzfileStream.open(subjectfiles[s].c_str());

		if (gzfileStream.get()!='>'){
			std::cerr << endl << "Error: " << subjectfiles[s] << " is not in the correct format. Expecting header to start with >." << endl << endl;
			continue;
		}

		while (std::getline(gzfileStream, header)){

			std::stringstream headerstream;
			headerstream << header;//convert the sequence to stringstream

			std::string alleleName;
			std::getline(headerstream, alleleName, ' ');
			alleleNames.push_back(alleleName);
		
			std::getline(gzfileStream, sequence, '>');

			sequence.erase(std::remove_if( sequence.begin(), sequence.end(), ::isspace ), sequence.end() );

			std::vector < std::string > alleleSequence ( querySampleNames.size() , string ( sequence.length() , '-' ) );
			std::vector < std::vector < bool > > alleleCoverage ( querySampleNames.size() , vector < bool > ( sequence.length() , false ) );
				
			if (sequence.length()<substringlength){
				continue;
			}
			
			std::transform(sequence.begin(), sequence.end(), sequence.begin(), ::toupper);//change all letters in the string to upper case
			
			int i=0;

			for (std::string::iterator iti = sequence.begin(), end = sequence.end()-(substringlength-1); iti != end; ++iti, ++i){
				
				std::string kmer=sequence.substr(i,substringlength);
				
				bool isrev=reverseComplementIfMin(kmer);
				extractMiddleBase(kmer, base);
				ascii_codons(kmer);
				
				std::unordered_map < std::string, std::string >::iterator it = kmerMap.find(kmer);//check if the kmer is in the hash
				
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
    				base=std::toupper(alleleSequence[k][i]);
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
			std::string locus;
			//cout << querySampleNames[k];
	 		for (std::vector < int >::iterator it=bestMatches[k].begin(); it!=bestMatches[k].end(); ++it){
	 			//cout << "\t" << alleleNames[*it];
				string::size_type found=alleleNames[*it].find_last_of("_");
				if (found!=string::npos){
					locus=alleleNames[*it].substr(0,found);
					allele=stoi(alleleNames[*it].substr(found+1));
				}
				std::unordered_map < std::string, std::vector < int > >::iterator it2 = profileMap[k].find(locus);//check if the kmer is in the hash
				if ( it2 != profileMap[k].end() ){//if the kmer is not in the hash
					it2->second.push_back(allele);
					uncertain[k]=true;
				}
				else{
					std::vector < int > newLocus;
					newLocus.push_back(allele);
					profileMap[k].insert(std::make_pair(locus, newLocus));
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
				std::ofstream alignfile(querySampleNames[k]+"_"+locus+".fa");
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

		std::vector < int > certainFiles;
		for (int k=0; k<querySampleNames.size(); ++k){
			if (novel[k].size()==0 && not uncertain[k]){
				certainFiles.push_back(k);
			}
		}

		if(openFileStream(profileFile, fileStream)){return 1;}

		std::string header;
		std::getline(fileStream, header);

		std::stringstream headerstream;
		headerstream << header;//convert the header to a stringstream
		std::string word;
		int i=0;
		int ST=-1;
		int cc=-1;
		std::map < std::string, int > alleles;
		int maxallelecolumn=0;

		while (std::getline(headerstream, word, '\t')){
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
			std::cerr << "the alleles in your profile file don't match those in your input files" << endl << endl;
			return 1;
		}

		for (std::map < std::string, int >::iterator it=alleles.begin(); it!=alleles.end(); ++it){
			std::unordered_map < std::string, std::vector < int > >::iterator it2 = profileMap[0].find(it->first);//check if the kmer is in the hash
			if ( it2 == profileMap[0].end() ){//if the kmer is not in the hash
				std::cerr << "the alleles in your profile file don't match those in your input files" << endl << endl;
				return 1;
			}
		}

		std::cout << "Sample" << "\t";
		for (std::map < std::string, int >::iterator it=alleles.begin(); it!=alleles.end(); ++it){
			std::cout << it->first << "\t";
		}
		std::cout << "ST" << endl;

		std::string line;
		int linenumber=0;
		int matches=0;
		while (std::getline(fileStream, line)){
			linenumber++;
			std::stringstream linestream;
			linestream << line;//convert the header to a stringstream
			std::vector < std::string > words;

			while (std::getline(linestream, word, '\t')){
				words.push_back(word);
			}

			if (words.size()<maxallelecolumn){
				std::cerr << "Malformed profile file on line " << linenumber << ". Expecting at least " << maxallelecolumn << " columns and found " << words.size() << endl << endl;
				return 1;
			}

			for (std::vector < int >::iterator k=certainFiles.begin(); k!=certainFiles.end(); ++k){
				bool perfectmatch=true;
				for (std::map < std::string, int >::iterator it=alleles.begin(); it!=alleles.end(); ++it){
					std::unordered_map < std::string, std::vector < int > >::iterator it2 = profileMap[*k].find(it->first);//check if the kmer is in the hash
					if ( it2 == profileMap[*k].end() ){//if the kmer is not in the hash
						perfectmatch=false;
						std::cerr << "Shouldn't be here";
					}
					else{
						if (std::stoi(words[it->second])!=it2->second[0]){
							perfectmatch=false;
							continue;
						}
					}
				}
				if (perfectmatch){
					std::vector < int > STs;
					STs.push_back(std::stoi(words[ST].substr(words[ST].find_last_of("_")+1)));
					profileMap[*k]["ST"]=STs;
				}
			}
			if (matches==certainFiles.size()){
				break;
			}
		}

		fileStream.close();
		

		for (int k=0; k<querySampleNames.size(); ++k){
			std::cout << querySampleNames[k] << "\t";
			for (std::map < std::string, int >::iterator it=alleles.begin(); it!=alleles.end(); ++it){
				std::string suffix="";
				std::set < std::string >::iterator it2 = novel[k].find(it->first);//check if the kmer is in the set
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
					std::cout << "-" << "\t";
				}
				else if (profileMap[k][it->first].size()==1){
					std::cout << profileMap[k][it->first][0] << suffix << "\t";
				}
				else {
					for (int i=0; i<profileMap[k][it->first].size()-1; ++i){
						std::cout << profileMap[k][it->first][i] << "/";
					}
					std::cout << profileMap[k][it->first][profileMap[k][it->first].size()-1] << suffix << "\t";
				}
			}
			if (profileMap[k]["ST"].size()==1){
				std::cout << profileMap[k]["ST"][0] << endl;
			}
			else {
				std::cout << "-" << endl;
			}
		}

	}
	else{
		std::set < std::string > allLoci;
		for (int k=0; k<querySampleNames.size(); ++k){
			for (std::unordered_map < std::string, std::vector < int > >::iterator it=profileMap[k].begin(); it!=profileMap[k].end(); ++it){
				allLoci.insert(it->first);
				//cout << it->first << endl;
			}
		}
		std::cout << "Sample";
		for (std::set < std::string >::iterator it=allLoci.begin(); it!=allLoci.end(); ++it){
			std::cout << "\t" << *it;
		}
		std::cout << endl;
		for (int k=0; k<querySampleNames.size(); ++k){
			std::cout << querySampleNames[k] << "\t";
			for (std::set < std::string >::iterator it=allLoci.begin(); it!=allLoci.end(); ++it){
				std::string suffix="";
				std::set < std::string >::iterator it2 = novel[k].find(*it);//check if the kmer is in the set
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
					std::cout << "-" << "\t";
				}
				else if (profileMap[k][*it].size()==1){
					std::cout << profileMap[k][*it][0] << suffix << "\t";
				}
				else {
					for (int i=0; i<profileMap[k][*it].size()-1; ++i){
						std::cout << profileMap[k][*it][i] << "/";
					}
					std::cout << profileMap[k][*it][profileMap[k][*it].size()-1] << suffix << "\t";
				}
			}
			std::cout << endl;
		}
	}

	return 0;
	
}


