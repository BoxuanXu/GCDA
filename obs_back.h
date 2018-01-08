/* *************************
 * Author: xbx1992
 * Created Time:  Tue 11  July 2015
 * LastModified:  Tue 11July 2015 04:00:27 PM CST
 * C File Name:  observeoperation.h
 * ************************/
#ifndef OBERVEOPERATION_H
#  define OBERVEOPERATION_H

#include <iostream>
#include <cstring>
#include <map>
#include <vector>
#include "common.h"
#include "configclass.h"
using namespace std;

void write_scaling_factor(const string& f, const vector<double>& v, daconfig dc, fileconfig fc, double tau_start, double tau_end);

void write_input_geos(const char* geos_output, const vector<observation>& obs, const datetime & start_date, const datetime& end_date_inter,const datetime& end_date);

class obs_operation
{
	public:
	obs_operation(){}
	void Initialize(daconfig dc,fileconfig fc, mpi_config mpic, int mpi_rank,statevector<double> scal_x_a,const vector<observation>& obs, const datetime & DA_start, const datetime & DA_end, const datetime & DA_end_lag)
	{
        	for(size_t i = 0; i < dc.K; ++i)
		{
	   	if(size_t(mpi_rank % mpic.size) == (i % mpic.size)){
			//printf("PROCESS %d start handling ensemble-%.3lu on %s...\n", mpic.rank, i, mpic.host_name);
			printf("PROCESS %d start handling ensemble-%.3lu...\n", mpi_rank, i);

               		char cmd_str[128];//string buffer of command 
      
                	/*
                 	 * go into the running directory of ensemble-i
                 	 */
		
              		sprintf(cmd_str, "mkdir %s/ensemble-%.3lu", fc.DA_dir, i );
                	cmd(cmd_str);
                
                	sprintf(cmd_str, "%s/ensemble-%.3lu/scaling_factor.geos.2x25", fc.DA_dir, i);
                	write_scaling_factor(cmd_str, scal_x_a.get_col(i), dc, fc, datetime(DA_start.str("YYYY-MM-DD 00:00:00")).tau(), DA_end_lag.tau());
        

                	/*
                 	 * prepare restart file
                 	 */
              		if(DA_start == dc.start_date)
	      		{
				sprintf(cmd_str, "cp %s %s/ensemble-%.3lu/", fc.restart_f, fc.DA_dir, i );
              			cmd(cmd_str);
	      		}
                	/* 
                 	 * write input.geos file
                 	 */
               		sprintf(cmd_str, "%s/ensemble-%.3lu/input.geos", fc.DA_dir, i);
                 	write_input_geos(cmd_str, obs, DA_start, DA_end, DA_end_lag);

                	/* 
                 	 * run geos model 
                 	 */
       	       		cout<<"start begining run ensemble-"<<i<<endl;
               		sprintf(cmd_str, "(cd %s/ensemble-%.3lu; %s ;cd ..;)",fc.DA_dir, i, fc.geos_f);
               		if(mpi_rank == 0)
               			cmd(cmd_str, true);
               		else
               			cmd(cmd_str, false);
			}
        	}// END OF GENERATING y_b for each ensemble member
		printf("PROCESS %d on %s FINISHED...\n", mpi_rank, mpic.host_name);
	}

	void Run_posterior(daconfig dc, fileconfig fc, vector<double> scal_posterior_x_a,OBS_DATA obs, datetime DA_start, datetime DA_end)
	{
			
		    /**********************
		     *  SAVE scal_x_a into file
		     **********************/
		    char buf_str[128];//string buffer of command 
		    sprintf(buf_str, "%s/posterior_scaling_factor.%s.geos.2x25", fc.DA_dir, DA_start.str("YYYYMMDDhhmmss").c_str());
		    write_scaling_factor(buf_str, scal_posterior_x_a, dc, fc, datetime(DA_start.str("YYYY-MM-DD 00:00:00")).tau(), DA_end.tau());
		     /*************************************************
		     *
		     * START POSTERIOR RUNNING
		     *
		     * ***********************************************/
		    printf("START POSTERIOR RUNNING...\n");
		    char cmd_str[128];//string buffer of command 
		    
		    /*
		     * go into the running directory of ensemble-i
		     */
		    sprintf(cmd_str, "mkdir %s/posterior", fc.DA_dir );
		    cmd(cmd_str);
		    sprintf(cmd_str, "%s/posterior/scaling_factor.geos.2x25", fc.DA_dir);
		    write_scaling_factor(cmd_str, scal_posterior_x_a, dc, fc, datetime(DA_start.str("YYYY-MM-DD 00:00:00")).tau(), DA_end.tau());
			
		    /*
		     * prepare restart file
		     */
		    sprintf(cmd_str, "cp %s %s/posterior/", fc.restart_f, fc.DA_dir);
		    cmd(cmd_str);
		    /* 
		     * write input.geos file
		     */
		    sprintf(cmd_str, "%s/posterior/input.geos",fc. DA_dir);
		    write_input_geos(cmd_str, obs, DA_start, DA_end, DA_end);
		    /* 
		     * run geos model 
		     */
		    sprintf(cmd_str, "(cd %s/posterior; %s )",fc.DA_dir, fc.geos_f);
		    cmd(cmd_str, true);
	}
};

void write_input_geos(const char* geos_output, const vector<observation>& obs, const datetime & start_date, const datetime& end_date_inter,const datetime& end_date){
    ifstream geos_conf;
    geos_conf.open("/data1/xubx/GEOS-Chem/tools/xbx_test_dir/input.geos", ifstream::in);
    ofstream out_conf;
    out_conf.open(geos_output, ifstream::out);
    if(!(geos_conf.is_open() && out_conf.is_open())){
        printf(" file open failed !\n");
        exit(1);
    }
    /*
     * get the lat/lon index of observation station
     */
    set<pair<pair<int, int>, int > > STATION_IJ;
    for(size_t i = 0; i < obs.size(); ++i){
    	pair<int,int> P = pair<int,int>(obs[i].geosI, obs[i].geosJ);
    	pair<pair<int,int>, int> PP;
    	PP.first = P,PP.second = obs[i].alt;
        STATION_IJ. insert(PP);
    }
    int station_number = STATION_IJ.size();
            
    char sbuf[1024];
    while(!geos_conf.eof()){
        geos_conf.getline(sbuf, 1024);
        /*
         * set up Start and End time
         */
        if(kmp_match(sbuf, "Start YYYYMMDD, HHMMSS") >= 0){
            sprintf(sbuf, "Start YYYYMMDD, HHMMSS  : %s",start_date.str("YYYYMMDD hhmmss").c_str()); 
        }
        else if(kmp_match(sbuf, "End   YYYYMMDD, HHMMSS") >= 0){
            sprintf(sbuf, "End   YYYYMMDD, HHMMSS  : %s",end_date.str("YYYYMMDD hhmmss").c_str()); 
        }
        /*
         * setup output restart file
         */
        /*else if(kmp_match(sbuf, "Schedule output for ") >= 0){
            char restart_buf[32] = {0};
            if(kmp_match(sbuf, monStr[end_date.month - 1]) >= 0){
                int m = end_date.monthmax();
                for(int i = 0; i < m; ++i) restart_buf[i] = '0';
                restart_buf[end_date.day - 1] = '3';
                if(end_date.month == end_date_inter.month)
		{
			restart_buf[end_date_inter.day - 1] = '3';
		}
		sprintf(sbuf, "Schedule output for %s : %s",monStr[end_date.month - 1], restart_buf); 
            }
	    if((end_date.month - 1) == end_date_inter.month)
            	if(kmp_match(sbuf, monStr[end_date_inter.month - 1]) >= 0)
	    	{
                	int m = end_date_inter.monthmax();
                	for(int i = 0; i < m; ++i) restart_buf[i] = '0';
                	restart_buf[end_date_inter.day - 1] = '3';
			sprintf(sbuf, "Schedule output for %s : %s",monStr[end_date_inter.month - 1], restart_buf); 
			
		}
        }*/
        /* 
         * Write ND48 menu
         */
        else if(kmp_match(sbuf, "ND48 MENU") >= 0){
            while(true){
                //printf("whiling...\n\t\t[%s]\n", sbuf);

                /* turn on ND48 */
                if(kmp_match(sbuf, "Turn on ND48") >= 0){
		    if(station_number > 0 )
                    	out_conf << "Turn on ND48 stations   : T" << endl;
		    else
		    	out_conf << "Turn on ND48 stations   : F" << endl;
                }

                /* set number of stations */
                else if(kmp_match(sbuf, "Number of stations") >= 0 ){
                    out_conf << "Number of stations      :   " << station_number << endl;
                    //out_conf << "Number of stations      :  0" << endl;
                    int i = 0;
	       	    for(set<pair<pair<int,int>, int> >::iterator it = STATION_IJ.begin(); it != STATION_IJ.end(); ++it){
                         if(i < 9)
			 	out_conf << "Station #"<<++i<<" (I,J,Lmax,N) : "<<it->first.first<<" "<<it->first.second<<" "<<it->second<<" 1" << endl;
			 else if(i < 99)
			 	out_conf << "Station #"<<++i<<"(I,J,Lmax,N) : "<<it->first.first<<" "<<it->first.second<<" "<<it->second<<" 1" << endl;
			 else 
	 						 
			 	out_conf << "Station#"<<++i<<"(I,J,Lmax,N) : "<<it->first.first<<" "<<it->first.second<<" "<<it->second<<" 1" << endl;
                    }
                }

                /* skip Station info */
                else if(kmp_match(sbuf, "Station #") >= 0){
                }

                /* write other lines */
                else{
                    out_conf << sbuf << endl;
                }

                geos_conf.getline(sbuf, 1024);

                if(kmp_match(sbuf, "-------------+------------") >= 0) break;
            }
        }
        out_conf << sbuf << endl;
    }

    out_conf.close();
    geos_conf.close();
}

void write_scaling_factor(const string& f, const vector<double>& v, daconfig dc, fileconfig fc, double tau_start, double tau_end){
    //debug("write");
    datablock db(dc.XSIZE,dc.YSIZE,1,1,1,1);
    db.modelres_lat = 180.0 / (dc.YSIZE - 1);
    db.modelres_lon = 360.0 / dc.XSIZE;
    db.tracer = 1;
    str_cpy(db.category, "SCL-FCT");
    double step = dc.DA_step / 60.0;
    double length = dc.DA_length / 60.0;
    bpch sf("CTM bin 02","scaling factor for assimilation");
    int k = 0;
    while(tau_start + step * k < tau_end){
        db.tau0 = db.tau1 = tau_start + step * k;
        for(int x = 0; x < dc.XSIZE; ++x){
            for(int y = 0; y < dc.YSIZE; ++y){
                //db[x][y][0] = v[y * dc.XSIZE + x + dc.XSIZE * dc.YSIZE * int((step * k) / length)];
                db[x][y][0] = v[y * dc.XSIZE + x + dc.XSIZE * dc.YSIZE * k];
            	
	    }
        }
        sf.push(db);
        ++ k;
    }
	sf.writeF(f.c_str());
}
#endif
