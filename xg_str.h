/*
 * $Id$
 */
#ifndef XG_STR_H
#  define XG_STR_H


namespace xg{
    string tolower(string x){
        transform(x.begin(), x.end(), x.begin(), ::tolower);
        return x;
    }
    string toupper(string x){
        transform(x.begin(), x.end(), x.begin(), ::toupper);
        return x;
    }
    string trim(const string& str){
       string::size_type pos = str.find_first_not_of(' ');
  	if(pos == string::npos)
	{
		return str;
	}	
	string::size_type pos2 =str.find_last_not_of(' ');
	if(pos2 != string::npos)
	{
		return str.substr(pos,pos2-pos+1);

	}
	return str.substr(pos);
    }
}

#endif /* ifndef XG_STR_H */

