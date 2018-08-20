#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <chrono> //timing
using namespace std;


string splitFileName(const string & str){
	string filename=str;
	string path="";
	string::size_type found=filename.find_last_of("/\\");
	if (found!=string::npos){
		path=str.substr(0,found);
		filename=str.substr(found+1);
	}
	return filename;
}

int openFileStream(const string & fileName, ifstream & fileStream, bool verbose){

	if (fileStream.is_open()){
		fileStream.close();
		fileStream.clear();
	}

	if (verbose){
		cout << "Reading " << fileName << endl;
	}

	fileStream.open(fileName, ios::in);
	if (fileStream.fail()) {
		cerr << "Failed to open " << fileName << endl << endl;
		return 1;
	}
	else{
		return 0;
	}
}

int fileToVector(const string & filename, vector <string> & fileargs){
	ifstream fileStream;
	fileStream.open(filename, ios::in);
	if (fileStream.fail()) {
		cerr << "Failed to open " << filename << endl << endl;
		return 1;
	}
	string word;
	while (fileStream >> word){
		if (word.length()>500){
			cerr << "Names > 500 characters are not allowed" << endl;
			return 1;
		}
		else if (word.length()>0){
			fileargs.push_back(word);
		}
	}
	return 0;
}

void printDuration(const chrono::high_resolution_clock::time_point start){
	chrono::high_resolution_clock::time_point finish = std::chrono::high_resolution_clock::now();
	chrono::duration<double> elapsed = finish - start;
	cout << "Done" << endl;
	cout << "Total time required: " << elapsed.count() << "s" << endl;
}