/* *************************
 * Author: xbx1992
 * Created Time:  Tue 11  July 2015
 * LastModified:  Tue 11July 2015 04:00:27 PM CST
 * C File Name:  sample.h
 * ************************/
#ifndef SAMPLE_H
#  define SAMPLE_H

#include "bpch.h"
#include "observationdata.h"
#include "xg_datetime.h"
#include "configclass.h"

int lookup_obs(const OBS_INDEX& index, int x, int y, double tau, double max_diff,int i);

class Sample{
	public:
	Sample(){}
	void Initialize(OBS_DATA& obs, datetime DA_start,datetime DA_end, daconfig dc, fileconfig fc)
	{
	    	vector<mod_data> mod;
	    	mod_data buf_mod;
	    	vector<size_t> index_v(obs.size(),-1); 
    		
		clock_t time_f = clock(),time_s;
            	
        	FILE *fp = fopen("/data1/xubx/GEOS-Chem/tools/main_experience_P_a_test/zhf_site_geos","r");
        	//FILE *fp = fopen("/data1/xubx/GEOS-Chem/tools/main_experience_ocean0.5_testregion/zhf_site_geos","r");
		vector<pair<pair<int, int>, int > > site;
		int lon,lat,alt;
		while(fscanf(fp, "%d %d %d", &(lon), &(lat), &(alt)) != EOF){	
       			pair<int,int> P = pair<int,int>(lon,lat);
			pair<pair<int,int>, int> PP;
			PP.first = P,PP.second = alt;
			site.push_back(PP);
		}
		
		for(size_t i = 0; i < dc.K; ++i) {
                	char cmd_str[128];//string buffer of command 
                	/*
                 	 * read modeled observation 
                 	 */
			
			//for(datetime start=DA_start;start < DA_end; start += timespan(0,0,0,0,dc.DA_step))
			//{
			sprintf(cmd_str, "%s/ensemble-%.3lu/%s0000.bpch", fc.DA_dir, i, DA_start.str("YYYYMMDD").c_str());
                    	printf("reading stations output :%s\n", cmd_str);
                    	bpch station(cmd_str);
                    	for(size_t j = 0; j < station.size(); ++j){
                        	if(kmp_match(station[j].category, "IJ-AVG-$") >= 0){
			    	//for(int i = 0;i<station[j].dim[0];i=i+10)
			    		//for(int m = 0;m<station[j].dim[1];m=m+10){
				for(size_t i = 0; i < site.size();++i){	
					buf_mod.geosI = site[i].first.first;
					buf_mod.geosJ = site[i].first.second;
					buf_mod.tau = station[j].tau0;
					buf_mod.value = station[j].data[site[i].first.first-1][site[i].first.second-1][site[i].second-1] * 1e6;
					mod.push_back(buf_mod);
					}
				}
                	}
		
			printf("finish read!"); 
    		
			time_s = Cal_time(time_f); 
			
			OBS_INDEX ans;
			for(size_t j = 0;j < mod.size();++j){
				ans[pair<int,int>(mod[j].geosI,mod[j].geosJ)].push_back(pair<double,int>(mod[j].tau,j));
			}
			

			for(OBS_INDEX::iterator it = ans.begin();it!=ans.end();++it){
				sort(it->second.begin(), it->second.end(), cmp_1st<double, int>);	
			}
			
			printf("finish short!"); 
    			time_f = Cal_time(time_s);
			
			for(size_t j =0; j < obs.size(); ++j){
				if(obs[j].tau <= (DA_start + timespan(0,0,0,0,dc.DA_length)).tau())
				{
				if(i == 0)
					index_v[j] = lookup_obs(ans,obs[j].geosI,obs[j].geosJ,obs[j].tau, 0.25, 0);
				if(index_v[j] >= 0)
					obs[j].modeled_v.push_back(mod[index_v[j]].value);
				}
				else
					break;	
			}
			mod.clear();
			
			vector<mod_data>(mod).swap(mod);
			
			printf("finish sample!"); 
    			time_s = Cal_time(time_f);
			//}
            	}
            	/*
             	 * DEBUG OUTPUT
             	 *
             	 * content of obs_index
             	 */
            	/*for(OBS_INDEX::iterator it = obs_index.begin(); it!=obs_index.end(); ++it){
                	cout << it->first.first<<","<<it->first.second<< "::"; 
                	for(size_t i = 0; i < it->second.size(); ++i)
                    	cout << it->second[i].first << ","<<it->second[i].second;
                	cout << endl;
            	}*/

            	/*
             	 * DEBUG OUTPUT
             	 *
             	 * content of obs
             	 */
            	for(size_t i = 0; i < obs.size(); ++i){
                	printf("obs-%lu:\n(%d, %d, %d, [%lf])-> %.10lf\n", i, obs[i].geosI, obs[i].geosJ, obs[i].alt,obs[i].tau, obs[i].value);
                	cout << obs[i].modeled_v << endl;
            	}
	
	}	
};

int lookup_obs(const OBS_INDEX& index, int x, int y, double tau, double max_diff = 0.25,int i=0){
    
    printf("look for tau = %lf...", tau);
    OBS_INDEX::const_iterator it = index.find(pair<int,int>(x,y));
    if(it != index.end()){
        int a = 0, b = it->second.size() - 1;
        while(a < b - 1){
            int mid = (a + b) / 2;
            if( it->second[mid].first > tau){
                b = mid;
            }
            else 
                a = mid;
        }
        double da = fabs(tau - it->second[a].first);
        double db = fabs(tau - it->second[b].first);
	//if(i==0)
	//	printf("x:%d,y:%d,a:%d,b:%d,da:%lf,db:%lf,tau:%lf\n",x,y,a,b,da,db,tau);
        if(da <= db && da <= max_diff){
            printf(" index = %d\n",it->second[a].second );
            return it->second[a].second;
        }
        else if(da > db && db <= max_diff){
            printf(" index = %d\n",it->second[b].second );
            return it->second[b].second;
        }
	else {
            printf(" index = -1\n");
            return -1;
        }
    }
    else {
        printf(" index = -1\n");
        return -1;
    }
}


#endif
