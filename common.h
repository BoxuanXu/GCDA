/*
 * $Id$
 */
#ifndef COMMON_H
#  define COMMON_H


#include <map>
#include <fstream>
//#include "armadillo"
#include "xg_str.h"
using namespace std;
//using namespace arma;

void cmd(string s, bool show = true){
	FILE *in;
	char buff[512];

	if(!(in = popen(s.c_str(), "r")))
        exit(1);

	while(fgets(buff, sizeof(buff), in)!=NULL){
        if(show) cout << buff;
	}
	pclose(in);
}

/*string cmd(string s,bool show =true){
	FILE *in;
	char buff[512];
	if(!(in = popen(s.cstr(),"r")))
	return "0";
	
 	while(fgets(buff,sizeof(buff),in)!=NULL){
	if(show) return buff;
	}
	
 	return "0";
}*/

void WriteConfig(const char * conf_file, const map<string, string>& config){
    string buf;
    ofstream f(conf_file);
    if(f.is_open()){
        for(map<string, string>::const_iterator it = config.begin(); it != config.end(); ++it){
            f << it->first << "=" << it->second << endl;
        }
        f.close();
    }
    else{
        //ERROR
    }
}

void ReadConfig(const char * conf_file, map<string, string>& config){
    string buf;
    ifstream f(conf_file);
    if(f.is_open()){
        while(getline(f, buf)){
            size_t p = buf.find_first_of("=");
            if(p != string::npos){
                config[trim(buf.substr(0, p))] = trim(buf.substr(p + 1));
            }
            else{
                // ERROR LINE
            }
        }
        f.close();
    }
    else{
        //ERROR
    }
}

void DebugConfig(const map<string,string> & config){
    for(map<string, string> :: const_iterator it = config.begin(); it != config.end(); ++it){
        clog << it->first << " : " << it->second << endl;
    }
}



#endif /* ifndef COMMON_H */

