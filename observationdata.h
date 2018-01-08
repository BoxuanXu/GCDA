/* *************************
 * Author: xbx1992
 * Created Time:  Tue 11  July 2015
 * LastModified:  Tue 11July 2015 04:00:27 PM CST
 * C File Name:  observationdata.h
 * ************************/
#ifndef OBERVATIONDATA_H
#  define OBERVATIONDATA_H

#include <iostream>
#include <cstring>
#include <map>
#include <vector>
#include "xg_datetime.h"
#include "configclass.h"
using namespace std;
class obs_c
{
	public:
	obs_c(){}
	OBS_DATA Initialize(const datetime & DA_start, const datetime & DA_end_lag,daconfig dc, fileconfig fc)
	{
		OBS_DATA obs;
        	FILE *fp = fopen(fc.obs_f,"r");
		observation buf_o;
        	char buf_s[128];
        	double dx = dc.grid_xres;
        	double dy = dc.grid_yres;
		while(fscanf(fp, "%s %lf %lf %d %lf %lf", buf_s, &(buf_o.lon), &(buf_o.lat), &(buf_o.alt), &(buf_o.value), &(buf_o.mdm)) != EOF){
		    //printf("bufs = %s\n", buf_s);
		    buf_o.t = datetime(buf_s, "YYYYMMDDhhmmss");
		    double tau = buf_o.t.tau();
		    if(tau >= DA_start.tau() && tau < DA_end_lag.tau()){
			buf_o.geosI = 1 + (buf_o.lon + 180) / dx;
			//buf_o.geosI = buf_o.lon;
			buf_o.geosJ = 1 + (buf_o.lat +  90) / dy;
			//buf_o.geosJ = buf_o.lat;
			obs.push_back(buf_o);
		    }
		}
        	fclose(fp);
		return obs;
	}
	
	MatrixXf Get_R(OBS_DATA& obs)
	{
		   /********************
 		    * error covariance
 		    * *****************/ 
		    MatrixXf R_test(obs.size(),obs.size());
		    for(size_t i = 0; i < obs.size(); ++i){
                        R_test(i,i) = pow(obs[i].mdm * 1e-6 ,2);

                    }
		    return R_test;	
	}
};


OBS_INDEX build_index(const vector<observation>& obs){
    OBS_INDEX ans;
    for(size_t i = 0; i < obs.size(); ++i){
        ans[pair<int,int>(obs[i].geosI, obs[i].geosJ)].push_back(pair<double, int>(obs[i].t.tau(), i));
    }
    for(OBS_INDEX::iterator it = ans.begin(); it!=ans.end(); ++it){
        sort(it->second.begin(), it->second.end(), cmp_1st<double, int>);
    }
    return ans;
}
#endif
