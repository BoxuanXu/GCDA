/* *************************
 * Author: xg1990
 * Created Time:  
 * LastModified:  Thu 05 Dec 2013 10:49:08 PM CST
 * C File Name: 
 * ************************/
#include <iostream>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include "bpch.h"
#include "xg_str.h"
#include "xg_datetime.h"
#include <map>
#include <boost/algorithm/string.hpp>
using namespace std;
#define out(v) cerr << #v << ": " << (v) << endl

const double R = 6378100; // [meter]
//const double R = 6.371e6; 
const double pi = 3.141592653; // [meter]
const double deg2rad = pi/180;
const double dxx=360.0/360*deg2rad;
const double dyy=180.0/180*deg2rad;
int main(int argc, char *argv[])
{
    bpch bpchmean(argv[1]);
    datetime start_date (argv[2], "YYYY-MM-DD");
    datetime end_date (argv[3], "YYYY-MM-DD");
    map<int,map<int, vector<double> > > x_y_v;
    double global_co2;
    datablock bpch_zmean(144,91,1,1,1,1);
    bpch_zmean.modelres_lat = 2;
    bpch_zmean.modelres_lon = 2.5;
    bpch_zmean.tracer = 1;
    str_cpy(bpch_zmean.category, "CO2-SRCE");
   
	for(size_t m = 0; m < bpchmean.size(); ++m){
                //if(kmp_match(bpchmean[m].category,"IJ-AVG-$")>=0){
			if(bpchmean[m].tau0 >= start_date.tau() && bpchmean[m].tau1 <= end_date.tau())
                	if(kmp_match(bpchmean[m].category,"CO2-SRCE")>=0){
	               	//if(bpchmean[m].tracer == 6)
			{
				/*if(bpchmean[m].tracer==8)
				{
					global_co2 = global_co2;
					cout<<global_co2<<endl;
					global_co2=0;
				}
				else{*/
				for(size_t x = 0;x < bpchmean[m].dim[0]; ++x)
					for(size_t y = 0;y < bpchmean[m].dim[1]; ++y)
					{
						double xpos = 2.5 * (x + 1)  - 181.25;
						double ypos = 2 * (y + 1) - 91;
						double area = abs( dxx*(sin(ypos*deg2rad)-sin(ypos*deg2rad-dyy*abs(2)))*R * R * 2.5);
						global_co2 = global_co2 + bpchmean[m].data[x][y][0] * area * 1e4 * 60 * 60 * 24 * 365 * 12e-15 / (6.022e23);
				
					}		//for(size_t z = 0;z < 47;++z)
				//}
	               		if(bpchmean[m].tracer == 10)
					cout<<global_co2<<endl;
			}
		}
    	}
					global_co2=0;
    return 0;
}


