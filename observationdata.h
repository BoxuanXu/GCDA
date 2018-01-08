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
		observation buf_o;
        	char buf_s[128];
        	double dx = dc.grid_xres;
        	double dy = dc.grid_yres;
        	FILE *fp = fopen("/data1/xubx/GEOS-Chem/tools/main_experience_P_a_test/zhf_site_geos","r");
		vector<pair<pair<int, int>, int > > site;
		int lon,lat,alt;
		while(fscanf(fp, "%d %d %d", &(lon), &(lat), &(alt)) != EOF){	
       			pair<int,int> P = pair<int,int>(lon,lat);
			pair<pair<int,int>, int> PP;
			PP.first = P,PP.second = alt;
			site.push_back(PP);
		}
			
                char cmd_str[128];//string buffer of command 
		sprintf(cmd_str, "/data1/xubx/GEOS-Chem/tools/generate_prior_co2_dir/apri_justldhy/%s0000.bpch", DA_start.str("YYYYMMDD").c_str());
		printf("reading stations output :%s\n", cmd_str);
		bpch station(cmd_str);
	
                for(size_t j = 0; j < station.size(); ++j){
			int ii = 0;
                        //vector<double> rands(551);
                        vector<double> rands(116);
			//vector_randn_boost(rands, 551, 0, 0.04, -0.4, 0.4);
			//vector_randn_boost(rands, 551, 0, 0.01, -0.1, 0.1);
			vector_randn_boost(rands, 116, 0, 0.04, 0.01, -0.01);
			int days = (station[j].tau0 - dc.start_date.tau()) / 24; 
			if(kmp_match(station[j].category, "IJ-AVG-$") >= 0)
			    //for(int i = 0;i<station[j].dim[0];i=i+10)
			    	//for(int m = 0;m<station[j].dim[1];m=m+10){
				for(size_t i = 0; i < site.size();++i){	
					buf_o.geosI = site[i].first.first;
					buf_o.geosJ = site[i].first.second;
					buf_o.alt   = site[i].second;
					buf_o.tau = station[j].tau0;
					buf_o.mdm = rands[ii]*10;
					buf_o.value = station[j].data[site[i].first.first-1][site[i].first.second-1][site[i].second-1] * 1e6 + rands[ii]*10;
					obs.push_back(buf_o);
					++ii;
				}
                	}
			return obs;
	}
	
	map<int, double> Get_R(OBS_DATA& obs)
	{
		   /********************
 		    * error covariance
 		    * *****************/ 
		    map<int, double> R_test;
		    for(size_t i = 0; i < obs.size(); ++i){
			    R_test[i] = pow(obs[i].mdm ,2);

                    }
		    return R_test;	
	}
};


OBS_INDEX build_index(const vector<observation>& obs){
    OBS_INDEX ans;
    for(size_t i = 0; i < obs.size(); ++i){
        ans[pair<int,int>(obs[i].geosI, obs[i].geosJ)].push_back(pair<double, int>(obs[i].tau, i));
    }
    for(OBS_INDEX::iterator it = ans.begin(); it!=ans.end(); ++it){
        sort(it->second.begin(), it->second.end(), cmp_1st<double, int>);
    }
    return ans;
}
#endif
