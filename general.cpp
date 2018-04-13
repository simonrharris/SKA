#include <string>
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