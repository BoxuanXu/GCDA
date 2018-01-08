/* *************************
 * Author: xbx1992
 * Created Time:  
 * LastModified:  Tue 11 Dec 2014 04:00:27 PM CST
 * C File Name: 
 * ************************/
#include <iostream>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <set>
#include "bpch.h"
#include "common.h"
#include "xg_datetime.h"
#include "state.h"
#include "init.h"
#include <Eigen/Dense>
#include "mpi.h"
#include "xg_math_vector.h"
//#ifdef _OPENMP
//#include "omp.h"
//#endif
using namespace Eigen;
using namespace std;

map<string, string> DA_config;


typedef map<pair<int, int>, vector<pair<double, int> > > OBS_INDEX;

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
        else if(kmp_match(sbuf, "Schedule output for ") >= 0){
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
        }
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

void write_scaling_factor(const string& f, const vector<double>& v, int XSIZE, int YSIZE, double tau_start, double tau_end, double tau_step, double tau_length){
    //debug("write");
    datablock db(XSIZE,YSIZE,1,1,1,1);
    db.modelres_lat = 180.0 / (YSIZE - 1);
    db.modelres_lon = 360.0 / XSIZE;
    db.tracer = 1;
    str_cpy(db.category, "SCL-FCT");
    
    bpch sf("CTM bin 02","scaling factor for assimilation");

    int k = 0;
    while(tau_start + tau_step * k < tau_end){
        db.tau0 = db.tau1 = tau_start + tau_step * k;
        for(int x = 0; x < XSIZE; ++x){
            for(int y = 0; y < YSIZE; ++y){
                db[x][y][0] = v[y * XSIZE + x + XSIZE * YSIZE * int((tau_step * k) / tau_length)];
            	
	    }
        }
        sf.push(db);
        ++ k;
    }
	sf.writeF(f.c_str());
}

typedef map<int, map<int, int > > REGIONS_MAP;
REGIONS_MAP regions,regions_index;

void read_Regions(const char *f)
{
	FILE *fp = fopen(f,"r");
	int c,line=90;
        
	while(fscanf(fp,"%d",&c) != EOF){
		regions[line][0] = c;
		for(size_t i=1;i<144;++i)
		{
			fscanf(fp,"%d",&c);
			regions[line][i] = c;
                        //debug(c);
                        //debug(line);
                        //debug(i);
			regions_index[c][line] = i;	
		}
		--line;
	}
	
}

template<typename mT>
MatrixXf m2M(matrix<mT>& m)
{
	MatrixXf M(m.nrow(),m.ncol());
	for(size_t i = 0;i < m.nrow();++i)
		for(size_t j = 0;j < m.ncol();++j)
			M(i,j) = m[i][j];
	return M;
}
template<typename mT>
matrix<mT> M2m(Eigen::MatrixXf& M)
{
	matrix<mT> m(M.rows(),M.cols());
	for(size_t i = 0;i < m.rows();++i)
		for(size_t j = 0;j < m.cols();++j)
			m[i][j] = M(i,j);
	
	return m;
}

map<int,MatrixXf> x_b_lagMembers;

template<typename mT>
void get_covariance(vector< vector<mT> >& COV,const int YSIZE)
{
	matrix<double> cov = COV;
	int a_x,a_y,b_x,b_y;
	for(size_t i = 0;i < cov.nrow();++i)
	{
		if(i  < 200)
			cov[i][i] = 0.64;
		else
			cov[i][i] = 0.16;
	}
	for(size_t i = 0;i < cov.nrow();++i)
		for(size_t j = 0;j < i;++j)
	{
		for(size_t x = 0;x < YSIZE;++x){
			if(regions_index[i][x] != 0){
				a_x =  x; 
				a_y = regions_index[i][x];
			}
		}	
		for(size_t x = 0;x < YSIZE;++x){
			if(regions_index[j][x] != 0){
				b_x =  x; 
				b_y = regions_index[j][x];
			}
		}	
		int L = sqrt(pow(a_x - b_x,2) +  pow(a_y - b_y,2));
		(i < 200)?cov[i][j] = cov[j][i] = 0.64 * exp(-L/7):cov[i][j] = cov[j][i] = 0.16 * exp(-L/10);

	}
	//debug(cov);
	COV = cov;
}

int main(int argc, char *argv[])
{
    /**********************************
     * Init MPI environment
     * *******************************/
    int mpi_size, mpi_rank, mpi_name_len;
    char mpi_host_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Init(&argc, &argv); // initialize MPI
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Get_processor_name(mpi_host_name, &mpi_name_len);
    /************************************
     * read configuration file 
     ************************************/
    ReadConfig(argv[1], DA_config);

    /* DA run dir */
    const char * DA_dir = DA_config["DA_dir"].c_str();

    /* load ensemble size   */
    size_t K = atoi(DA_config["ensemble"].c_str());// ensemble size

    /* load grid information*/
    double grid_xres = atof(DA_config["res_lon"].c_str());
    double grid_yres = atof(DA_config["res_lat"].c_str());
    int    XSIZE     = int(360 / grid_xres + 1e-5);
    int    YSIZE     = int(180 / grid_yres + 1e-5) + 1;

    /* load start/end time */
    datetime start_date (DA_config["start_time"], DA_config["date_format"]);
    datetime end_date   (DA_config["end_time"],   DA_config["date_format"]);

    /* load time step   */
    int    DA_step   = atoi(DA_config["DA_step"].c_str()); // assimilation time step     [minute]
    int    DA_length = atoi(DA_config["DA_length"].c_str()); // assimilation step length   [minute]
    size_t F         = DA_length / DA_step;        // number of scaling factor for each assimilation step
    /*nlags*/
    int    nlag = 5; // number of lag	
    //nlag = atoi(DA_config["nlag"].c_str()); // assimilation step length   [minute]
    debug(nlag);
    /* load inflation factor*/
    double rho = atof(DA_config["inflation"].c_str()); // inflation factor.
    /* load  REGIONS*/
    read_Regions(argv[2]);
    /* load Region 22's center*/
    vector<pair<int,int> > R_c;
    int a,b,c;
    FILE *fp2 = fopen(argv[3],"r");
    while(fscanf(fp2, "%d %d %d", &(a), &(b), &(c)) != EOF){
	    R_c.push_back(make_pair<int,int>(b,c));
    }
    fclose(fp2);
    /*****************************************************
     * Set up the initial state vector x (SCALING FACTOR)
     *****************************************************/
    int Class = 247;
    size_t M = Class * nlag,M1= XSIZE * YSIZE * nlag; // size of state vector;
    vector<double> x_init(M, 1.0);
    string restart_f = DA_config["restart_file"];
    /*
     * generate the first ensemble of state vector x 
     *
     *  DONE by PROCESS-0, then broadcast to MPI_COMM_WORLD
     */
    statevector<double> x_b_lagMembers(nlag , Class, K),scal_x_a(nlag, XSIZE * YSIZE, K);
    MatrixXf x_b(M, K);
    matrix<double> cov(Class, Class);
    vector<double> rand_t(M),x_b_p_tre(M1 * K);
    vector<double> rand_s(M);
    //double x_b_p[M * K];
    double *x_b_p = (double*)malloc(M*K*sizeof(double));
    int P;
    if(mpi_rank == 0) {
    	    get_covariance(cov,YSIZE);
	    x_b_lagMembers.Initialization();
	    x_b_p = x_b_lagMembers.expand();
	    MPI_Bcast(x_b_p, M * K, MPI_DOUBLE, 0, MPI_COMM_WORLD);
   }
   else{
    
    	MPI_Bcast(x_b_p, M * K, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	for(size_t j = 0;j< M * K;++j)
    	       x_b_lagMembers((j % (Class * nlag)) / Class, (j % (Class * nlag)) % Class, j / (Class * nlag)) = x_b_p[j];
   }
    //free(x_b_p);
    /************************************************
     *
     * START first year
     *
     ************************************************/
    /************************************************
     *
     * START ASSIMILATION CYCLE
     *
     ************************************************/
    for(datetime DA_i = start_date; DA_i < end_date; DA_i += timespan(0,0,0,0,DA_length)){

        /* set up time */
 	datetime DA_start = DA_i;
        datetime DA_end   = DA_i + timespan(0,0,0,0, DA_length);
        datetime DA_end_lag   = DA_i + timespan(0,0,0,0, DA_length*nlag);
        double   tau_start= DA_start.tau();
        double   tau_end  = DA_end  .tau();
	double   tau_end_lag =DA_end_lag.tau();
        printf("DA_length = %d, START ASSIMILATION CYCLE %s(%lf) - %s(%lf)\n", DA_length, DA_start.str().c_str(), tau_start, DA_end.str().c_str(), tau_end);
        /************************
         * read observation data
         *
         * which is y_o
         ************************/
	vector<observation> obs;
	char file_name[128];
	
        FILE *fp = fopen(DA_config["obs"].c_str(),"r");
	observation buf_o;
        char buf_s[128];
        double dx = atof(DA_config["res_lon"].c_str());
        double dy = atof(DA_config["res_lat"].c_str());
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
        /* build up index to speed up searching */
        OBS_INDEX obs_index = build_index(obs);
	/****************************
	 * statevector to grid
	 *
	 * *****************************/
	    char mk_str[128];
            
    	if(mpi_rank == 0) {
	    if(DA_i == start_date)
	    {
       	    	sprintf(mk_str, "mkdir %s/diagnose", DA_dir);
            	cmd(mk_str);
	    }
            sprintf(mk_str, "statevector.%s.nc", DA_i.str("YYYYMMDDhhmmss").c_str());
    	    bool flag = x_b_lagMembers.save_apri_vector(mk_str);	
	    }   
	    x_b = x_b_lagMembers.transfer(); 
	    for(size_t lag = 0; lag < nlag; ++lag) 
            {    
    		for(size_t i = 0; i < K; ++i)
		{	
			for(size_t x = 0;x < XSIZE;++x){
				for(size_t y = 0;y < YSIZE;++y){
					int index = XSIZE * y + x,index2 = regions[y][x];       
					scal_x_a(lag,index,i) = x_b_lagMembers(lag,index2,i);
		    		}
			}
		}
	    }
	 
	/***************************
	 * run ensemble
	 *
	 *************************/
        for(size_t i = 0; i < K; ++i)
	{ 
            /*
             * run geos on each process
             */
		   if(size_t(mpi_rank % mpi_size) == (i % mpi_size)){
			printf("PROCESS %d start handling ensemble-%.3lu on %s...\n", mpi_rank, i, mpi_host_name);

               char cmd_str[128];//string buffer of command 
      
                /*
                 * go into the running directory of ensemble-i
                 */
		
              sprintf(cmd_str, "mkdir %s/ensemble-%.3lu", DA_dir, i );
                cmd(cmd_str);
                
		//debug(scal_x_a.get_col(i));
                sprintf(cmd_str, "%s/ensemble-%.3lu/scaling_factor.geos.2x25", DA_dir, i);
                write_scaling_factor(cmd_str, scal_x_a.get_col(i), XSIZE, YSIZE, datetime(DA_start.str("YYYY-MM-DD 00:00:00")).tau(), tau_end_lag, DA_step / 60.0, DA_length / 60.0);
        

                /*
                 * prepare restart file
                 */
              if(DA_i == start_date)
	      {
			sprintf(cmd_str, "cp %s %s/ensemble-%.3lu/", restart_f.c_str(), DA_dir, i );
              		cmd(cmd_str);
	      }
		//else
	      //		sprintf(cmd_str, "cp -r %s/posterior/restart.%s %s/ensemble-%.3lu/", DA_dir, DA_start.str("YYYYMMDDhhmm").c_str(), DA_dir ,i);
                /* 
                 * write input.geos file
                 */
               sprintf(cmd_str, "%s/ensemble-%.3lu/input.geos", DA_dir, i);
                 write_input_geos(cmd_str, obs, DA_start, DA_end, DA_end_lag);

                /* 
                 * run geos model 
                 */
       	       cout<<"start begining run ensemble-"<<i<<endl;
               sprintf(cmd_str, "(cd %s/ensemble-%.3lu; %s ;cd ..;)",DA_dir, i, DA_config["geos"].c_str());
               if(mpi_rank == 0)
               	cmd(cmd_str, true);
               else
               	cmd(cmd_str, false);
		}

             }// END OF GENERATING y_b for each ensemble member
	printf("PROCESS %d on %s FINISHED...\n", mpi_rank, mpi_host_name);
	
   	 MPI_Barrier(MPI_COMM_WORLD);
	 
	
	if(mpi_rank == 0){
	  vector<mod_data> mod;
	    mod_data buf_mod;
	    	
            for(size_t i = 0; i < K; ++i) {
                char cmd_str[128];//string buffer of command 
                /*
                 * read modeled observation 
                 */
		    sprintf(cmd_str, "%s/ensemble-%.3lu/stations.%s", DA_dir, i,DA_start.str("YYYYMMDD").c_str());
                    printf("reading stations output :%s\n", cmd_str);
                    bpch station(cmd_str);
                    for(size_t j = 0; j < station.size(); ++j){
                        if(kmp_match(station[j].category, "IJ-AVG-$") >= 0){
				buf_mod.geosI = station[j].dim[3];
				buf_mod.geosJ = station[j].dim[4];
				buf_mod.tau = station[j].tau0;
				buf_mod.value = station[j].data[0][0][station[j].dim[2] - 1];
				mod.push_back(buf_mod);
				}
                	}
		OBS_INDEX ans;
		for(size_t j = 0;j < mod.size();++j){
			ans[pair<int,int>(mod[j].geosI,mod[j].geosJ)].push_back(pair<double,int>(mod[j].tau,j));
			
			//cout<<mod[j].geosI<<" "<<mod[j].geosJ<<" "<<mod[j].tau<<" "<<mod[j].value<<endl;
    			//OBS_INDEX::const_iterator it = ans.find(pair<int,int>(mod[j].geosI,mod[j].geosJ));
		        //cout<<"size:"<<it->second.size()<<"second:"<<it->second[0].second<<endl;
		}
		mod.clear();
		for(OBS_INDEX::iterator it = ans.begin();it!=ans.end();++it){
			sort(it->second.begin(), it->second.end(), cmp_1st<double, int>);	
		}
		for(size_t j =0; j < obs.size(); ++j){
			int idx_mod = lookup_obs(ans,obs[j].geosI,obs[j].geosJ,obs[j].t.tau());
			if(idx_mod >= 0){
				obs[j].modeled_v.push_back(mod[idx_mod].value);
			}	
		}
		/*
                 * remove the run dir of ensemble-i
                 */
                sprintf(cmd_str, "rm -rf %s/ensemble-%.3lu", DA_dir, i );
                //cmd(cmd_str);
            }
            /*
             * DEBUG OUTPUT
             *
             * content of obs_index
             */
            for(OBS_INDEX::iterator it = obs_index.begin(); it!=obs_index.end(); ++it){
                cout << it->first.first<<","<<it->first.second<< "::"; 
                for(size_t i = 0; i < it->second.size(); ++i)
                    cout << it->second[i].first << ","<<it->second[i].second;
                cout << endl;
            }

            /*
             * DEBUG OUTPUT
             *
             * content of obs
             */
            for(size_t i = 0; i < obs.size(); ++i){
                printf("obs-%lu:\n(%d, %d, %d, %s[%lf])-> %.10lf\n", i, obs[i].geosI, obs[i].geosJ, obs[i].alt, obs[i].t.str().c_str(),obs[i].t.tau(), obs[i].value);
                cout << obs[i].modeled_v << endl;
            }
           /*
             * CALCULATE X_b
             */
            MatrixXf X_b(M,K);
	    X_b = x_b;
	    debug(x_b);
	    MatrixXf x_b_bar(M,1);
	    x_b_bar = x_b.rowwise().mean();			    
	    for(size_t i = 0; i < X_b.cols(); ++i)
            {
				
				for(size_t j = 0; j < x_b_bar.rows();++j){
					X_b(j,i) = X_b(j,i) - x_b_bar(j,0); 
				} 
	    }  
	    debug(X_b);
	    MatrixXf P_b =  (X_b * X_b.transpose()) / (K - 1);
            /*
             * BUILD INDEX of obs with sufficient y_b
             */
            MatrixXf  y_o(obs.size(), 1),y_b(obs.size(),K),y_b_bar(obs.size(),1),y_b_bar_met(obs.size(),1);
	    for(size_t i = 0; i < obs.size(); ++i){
                	
		y_o(i,0) = obs[i].value;
		if(obs[i].modeled_v.size() == K){
		     for(size_t j=0;j < K;++j)
			y_b(i,j) = obs[i].modeled_v[j];
               }
	    }
	    y_b_bar_met = y_b_bar = y_b.rowwise().mean();
	    debug(y_b_bar);
	    MatrixXf Y_b(y_b);
	    
		for(size_t i = 0;i < Y_b.rows(); ++i)
            		for(size_t j = 0; j < Y_b.cols(); ++j)
			Y_b(i,j) = Y_b(i,j) - y_b_bar(i,0); 
            debug(Y_b);
	    
            /*
             *start assimilate GRID POINT by GRID POINT 
             */
            debug("START ASSIMILATING REGION");
		   /********************
 		    * error covariance
 		    * *****************/ 
		    MatrixXf R_test(obs.size(),obs.size());
		    for(size_t i = 0; i < obs.size(); ++i){
                        R_test(i,i) = pow(obs[i].mdm * 1e-6 ,2);

                    }
		     /******************
 		     * Pre for obs
 		     * *****************/ 		
		    MatrixXf x_b_inflation(M,K),PHt(M,1),HPHt(1,1),K_gain(M,obs.size()),alpha(M,M);
		    map<int,map<int,double> > HPHR;
		    int index_x;
		    K_gain = MatrixXf::Zero(M,obs.size());
		    alpha = MatrixXf::Zero(M,M);
		    for(size_t n = 0;n < obs.size(); ++n){
		     	   /******************
                     	    * STEP 0 get index_x
                            * ***************/
                          	index_x =int(int(obs[n].t.tau() - tau_start) / (24 * 7));
		     	   /******************
                     	    * STEP 1 check Y_b
                            * ***************/
				if(Y_b(n,0) == 0 || (y_b_bar(n,0) - obs[n].value) >= (2 * obs[n].mdm) )
					continue;
		     	   /******************
                     	    * STEP 2 calculate PHt
                            * ***************/
			    PHt = X_b * (Y_b.row(n).transpose());
			   for(size_t i = 0;i < PHt.cols();++i)
				for(size_t j = 0; j < PHt.rows();++j)
					PHt(j,i) = PHt(j,i) / ( K - 1);
                           debug(PHt);
    		   	   /******************************
     		            * generate Coef
     		            ****************************/
		            debug("before PHt");
			    int index_r = regions[obs[n].geosJ][obs[n].geosI];
			    int a_x,a_y,b_x,b_y;
			    for(size_t x = 0;x < YSIZE;++x){
				if(regions_index[index_r][x] != 0){
					a_x =  x; 
					a_y = regions_index[index_r][x];
				}
			    }	
		            double L;
			    MatrixXf Coef(M,1);
			    for(size_t i = 0; i < M; ++i) 
			    {
				    Coef(i,0) = 0;
				    if((i / Class) == index_x)
				    {
			    		for(size_t x = 0;x < YSIZE;++x){
						if(regions_index[i % Class][x] != 0){
							b_x =  x; 
							b_y = regions_index[i % Class][x];
						}
			    		}	
					L = sqrt(pow(a_x - b_x,2) +  pow(a_y - b_y,2));
				 	((i % Class) < 200)?Coef(i,0) = exp(-L/7):Coef(i,0) = exp(-L/10);
				    }
			   }
			   Coef(index_r+index_x * Class) = 1;
			   debug(Coef);
			    
		     	   /******************
                            * STEP 3 calculate HPHR
                            * ***************/
			    MatrixXf a(Y_b.row(n));
			    HPHt = a * a.transpose();
                            HPHR[n][n] = HPHt(0,0) / (K - 1) + R_test(n,n);
			    debug(HPHt);
			    debug(HPHR[n][n]);
		     	   /******************
                            * STEP 4 calculate K_gain and post_x_b
                            * ***************/
			   for(size_t i = 0;i < K_gain.rows();++i)
			   {
				K_gain(i,n) = PHt(i,0) / HPHR[n][n];
				x_b_bar(i,0) = x_b_bar(i,0) + Coef(i,0) * K_gain(i,n) * (y_o(n,0) - y_b_bar(n,0));
			   }	
		           //debug(K_gain);
		           //debug(x_b_bar);
			  alpha(index_r,index_r) = 1.0 / (1.0 + sqrt(R_test(n,n) / HPHR[n][n])); 
			   //debug(alpha);
			   /*****************
 			   * update the deviations from the mean state vector
 			   * ***************/ 
			   for(size_t k = 0;k < K;++k)
                           {
				for(size_t i = 0;i < M;++i)
					X_b(i,k) = X_b(i,k) - Coef(i,0) * alpha(index_r, index_r) * K_gain(i,n) *(Y_b(n,k));
                           } 
			   //debug(X_b);
		    	   /*****************
		      	    * STEP 6 inflation x_b
		     	    ****************/
		    		for(size_t i = 0;i < M; ++i)
		    	    		for(size_t j = 0;j < K; ++j)
					x_b_inflation(i,j) = x_b_bar(i,0) + (1 + Coef(i,0) * (rho - 1)) * X_b(i,j); 
			   /***************
 			   * updat the ensemble of sampled CO2 concentrations
 			   * **************/  
			   int length_L = 0;
			   for(size_t m = n+1;m < obs.size();++m)
			   {
				L = sqrt(pow(obs[n].geosJ - obs[m].geosJ,2) +  pow(obs[n].geosI - obs[m].geosI,2));
			    	(regions[obs[m].geosJ][obs[m].geosI] < 200)?:length_L=7;length_L=10;
				MatrixXf a2(Y_b.row(n)),b2(Y_b.row(m)),Y_b_term(1,1);
			    	Y_b_term = a2 * b2.transpose();
				double fac =(Y_b_term(0,0) / HPHR[n][n]) / ( K - 1);
				y_b_bar(m,0) = y_b_bar(m,0) + exp(-L/length_L) * fac * (y_o(n,0) - y_b_bar(n,0));
				
				for(size_t k = 0;k < K;++k)
				Y_b(m,k) = Y_b(m,k) - exp(-L/length_L) * alpha(index_r,index_r) * fac *(Y_b(n,k));
			   }
			   //debug(y_b_bar);
			   //debug(Y_b);
			   for(size_t m = 0;m < n + 1;++m)
			   {
				L = sqrt(pow(obs[n].geosJ - obs[m].geosJ,2) +  pow(obs[n].geosI - obs[m].geosI,2));
			    	(regions[obs[m].geosJ][obs[m].geosI] < 200)?:length_L=7;length_L=10;
			    	MatrixXf a2(Y_b.row(n)),b2(Y_b.row(m)),Y_b_term(1,1);
			    	 Y_b_term = a2 * b2.transpose();
				double fac = (Y_b_term(0,0) / HPHR[n][n]) / (K - 1);
				y_b_bar(m,0) = y_b_bar(m,0) + exp(-L/length_L) * fac * (y_o(n,0) - y_b_bar(n,0));
				
				for(size_t k = 0;k < K;++k)
				Y_b(m,k) = Y_b(m,k) - exp(-L/length_L) * alpha(index_r,index_r) * fac *(Y_b(n,k));
			   }
			   debug(n);
		    }
		    /*****************
		     * STEP 4 x_b
		     ****************/
		    for(size_t i = 0;i < M; ++i)
		    	    for(size_t j = 0;j < K; ++j)
				x_b(i,j) = x_b_bar(i,0) + X_b(i,j); 
				//debug(x_b);
				x_b_bar = x_b.rowwise().mean();
				debug(x_b_bar);
		    		x_b = x_b_inflation;
				debug(x_b);
				
		    /*****************
 		     *
                     * STEP 5 print Y_b_bar and obs
                     * ***************/
		    for(size_t i = 0;i < y_b_bar.rows(); ++i)
		    	for(size_t j = 0;j < K; ++j)
			    Y_b(i,j) = y_b_bar(i,0) + Y_b(i,j); 
		    y_b_bar = Y_b.rowwise().mean();	    
		
		    for(size_t i = 0;i < y_b_bar.rows(); ++i)
			    printf("%d %d %s %.10lf %.10lf %.10lf\n",obs[i].geosI,obs[i].geosJ,obs[i].t.str().c_str(),obs[i].value,y_b_bar(i,0),y_b_bar_met(i,0));
		    
		     /*****************:
 		     *
                     * STEP 5 generate all_scaling
                     * ***************/
		    debug("generate all scaling");
            	    vector<double> scal_posterior_x_a(XSIZE * YSIZE);
		    for(size_t x = 0;x < XSIZE;++x){
		    		for(size_t y = 0;y < YSIZE;++y){
				int index = XSIZE * y + x,index2 = regions[y][x];       
	                   		scal_posterior_x_a[index] = x_b_bar(index2,0);   
				}
            	    	}
		    debug(scal_posterior_x_a);
		    /******************
 		     *STEP 6 Propagate
 		     *****************/
		    x_b_lagMembers = reverse_transfer(x_b,nlag,Class);	    	    
                    debug("save post");
		    sprintf(mk_str, "statevector.%s.nc", DA_i.str("YYYYMMDDhhmmss").c_str());
    	    	    bool flag = x_b_lagMembers.save_post_vector(mk_str);	
            	    sprintf(mk_str, "mv statevector.%s.nc %s/diagnose/",DA_i.str("YYYYMMDDhhmmss").c_str(), DA_dir);
            	    cmd(mk_str);
		    
		    debug("propagate");

		    bool rt = x_b_lagMembers.Propagate();
    		    
		    x_b_p = x_b_lagMembers.expand();
    		    //debug(x_b_p2); 
		    MPI_Bcast(x_b_p, M * K, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	            /****************
 		    *
 		    * print P_a
 		    * ***************/ 
		    debug("P_a");
		   // MatrixXf tmp_P_a = (MatrixXf::Identity(M,M) + X_b * Y_b.transpose() * R_test.inverse() * Y_b * X_b.transpose());
		   //MatrixXf P_a = MatrixXf::Identity(M,M) + alpha * X_b * Y_b.transpose * R_test.inverse() * Y_b;
 	    	   /*MatrixXf Y_b_tmp = Y_b.tranpose() * Y_b;
		   for(size_t i = 0;i < M;++i)
	    	    	for(size_t j = 0;j < M;++j)
		    	{
				if(i == j)
					P_a = sqrt(1 + Coef(i,0) *  Y_b_tmp(i,j) / R_test(i,i)) * P_b[i][j];
				else
					P_a[i][j] = P_b[i][j];
			}
		    debug(P_a);*/
		    
		    /**********************
		     *  SAVE scal_x_a into file
		     **********************/
		    char buf_str[128];//string buffer of command 
		    sprintf(buf_str, "%s/posterior_scaling_factor.%s.geos.2x25", DA_dir, DA_i.str("YYYYMMDDhhmmss").c_str());
		    write_scaling_factor(buf_str, scal_posterior_x_a, XSIZE, YSIZE, datetime(DA_start.str("YYYY-MM-DD 00:00:00")).tau(), tau_end, DA_step / 60.0, DA_length / 60.0);
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
		    sprintf(cmd_str, "mkdir %s/posterior", DA_dir );
		    cmd(cmd_str);
		    sprintf(cmd_str, "%s/posterior/scaling_factor.geos.2x25", DA_dir);
		    write_scaling_factor(cmd_str, scal_posterior_x_a, XSIZE, YSIZE, datetime(DA_start.str("YYYY-MM-DD 00:00:00")).tau(), tau_end, DA_step / 60.0, DA_length / 60.0);
			
		    /*
		     * prepare restart file
		     */
		    sprintf(cmd_str, "cp %s %s/posterior/", restart_f.c_str(), DA_dir);
		    cmd(cmd_str);
		    /* 
		     * write input.geos file
		     */
		    sprintf(cmd_str, "%s/posterior/input.geos", DA_dir);
		    write_input_geos(cmd_str, obs, DA_start, DA_end, DA_end);
		    /* 
		     * run geos model 
		     */
		    sprintf(cmd_str, "(cd %s/posterior; %s )",DA_dir, DA_config["geos"].c_str());
		    cmd(cmd_str, true);
		}	
   		else{
    
		MPI_Bcast(x_b_p, M * K, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		for(size_t j = 0;j< M * K;++j)
    	       		x_b_lagMembers((j % (Class * nlag)) / Class, (j % (Class * nlag)) % Class, j / (Class * nlag)) = x_b_p[j];
   		}
   	 	MPI_Barrier(MPI_COMM_WORLD);
	 


	    }// END OF DA CYCLE
		
	    MPI_Finalize();
	    return 0;
	}

			    


