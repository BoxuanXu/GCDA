/************************
 * Author:xubx1992
 * Created Time:
 ************************/
#include <iostream>
#include <cstdio>
#include <set>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include "bpch.h"
#include <map>
using namespace std;

int main(int argc,char *argv[])
{
	//int x = lexical_cast<int>(argv[2]);
	//int y = lexical_cast<int>(argv[3]);
	//int z = lexical_cast<int>(argv[4]);
	//datetime s_date (argv[2], "YYYY-MM-DD");
	//datetime e_date (argv[3], "YYYY-MM-DD");

	bpch bpchout(argv[1]);
	vector<double> data;
	for(size_t m=0;m<bpchout.size();++m){
                if(kmp_match(bpchout[m].category,"SCL-FCT")>=0){
                //if(kmp_match(bpchout[m].category,"IJ-AVG-$")>=0){
			for(size_t x = 0;x < bpchout[m].dim[0];++x)
			for(size_t y = 0;y < bpchout[m].dim[1];++y)
			for(size_t z = 0;z < bpchout[m].dim[2];++z)
				printf("%d %d %d %d %lf %.9lf\n",x,y,bpchout[m].dim[3],bpchout[m].dim[4],bpchout[m].tau0,bpchout[m].data[x][y][z]);
		
				}
			}
			return 0;
}
