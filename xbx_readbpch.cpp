/************************
 * Author:xubx1992
 * Created Time:
 ************************/
#include <iostream>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include "bpch.h"
#include "xg_datetime.h"
#include "xg_str.h"
#include <boost/lexical_cast.hpp>
using boost::lexical_cast;
#define out(v) cerr << #v << ":" << (v) <<endl

int main(int argc,char *argv[])
{
	bpch bpchout(argv[1]);
	/*int x = lexical_cast<int>(argv[2]);
	int y = lexical_cast<int>(argv[3]);
	int z = lexical_cast<int>(argv[4]);*/
	datetime start_date (argv[2], "YYYY-MM-DD");
	datetime end_date (argv[3], "YYYY-MM-DD");
	vector<double> data;
	//cout<<bpchout.size()<<endl;
	for(size_t m=0;m<bpchout.size();++m){
                if(kmp_match(bpchout[m].category,"IJ-AVG-$")>=0){
                //if(kmp_match(bpchmean[m].category,"CO2-SRCE")>=0){
                //if(kmp_match(bpchout[m].category,"SCL-FCT")>=0){
                        //for(size_t i=0;i<bpchout[m].dim[0];i=i+1)
                        //{
                          //      for(size_t j=0;j<bpchout[m].dim[1];j=j+1){
					//for(size_t z=0;z<bpchout[m].dim[2];z=z+1){
                                        	//printf("%d %d %d %lf %lf %.10lf %s\n",bpchout[m].dim[3],bpchout[m].dim[4],bpchout[m].dim[2],bpchout[m].tau0,bpchout[m].tau1,bpchout[m].data[i][j][bpchout[m].dim[2]-1],bpchout[m].unit);
                                        	//printf("%d %d %d %lf %lf %.10lf %s\n",bpchout[m].dim[3],bpchout[m].dim[4],bpchout[m].dim[2],bpchout[m].tau0,bpchout[m].tau1,bpchout[m].data[x][y][z-1],bpchout[m].unit);
					//}
			//	}
                        //}
			if(bpchout[m].tau0 >= start_date.tau() && bpchout[m].tau1 <= end_date.tau())
			{
                                for(size_t x=0;x<bpchout[m].dim[0];x=x+1)
                                	for(size_t y=0;y<bpchout[m].dim[1];y=y+1)
						data.push_back(bpchout[m].data[x][y][0]);
					//for(size_t z=0;z<bpchout[m].dim[2];z=z+1)
				//cout<<bpchout[m].data[x][y][z]<<endl;
				printf("%.3f\n",mean(data) * 1e6);
				/*double max,min;
				min = data[0];
				for(int i = 0;i < data.size();++i)
				{
					
					if(max<data[i])
							max = data[i];
					if(min>data[i])
							min = data[i];

				}
				printf("min is:%.3f,max is:%.3f",min * 1e6,max * 1e6);*/
			}
		//cout<<bpchout.size()<<endl;
        }
	}
	return 0;
}
