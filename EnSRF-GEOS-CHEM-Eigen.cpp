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
#include <Eigen/Dense>
#include "mpi.h"
#include <netcdfcpp.h>
#include "xg_math_vector.h"
//#ifdef _OPENMP
//#include "omp.h"
//#endif
using namespace Eigen;
using namespace std;

map<string, string> DA_config;

struct observation{
    double lat, lon;
    int geosI,geosJ,alt;
    datetime t;
    double tau;
    double value,mdm;
    vector<double>modeled_v;
};
struct mod_data{
    double geosI,geosJ;
    double tau;
    double value,mdm;
};

		bool save_apri_vector(const char * file_s,MatrixXf x_b,int K)
		{

			NcFile s_f(file_s, NcFile::Replace);
			if(!s_f.is_valid())
			{
				cout<<"couldn't open file"<<endl;
				return 0;
			}
			s_f.add_att("disclasimer","This data is The GCDA's statevector");
			s_f.add_att("Auther","BOXUAN XU");
			s_f.add_att("Email","xuboxuan6@163.com");
			s_f.add_att("Source","GCDA release 1.0");



			NcDim* nmembers = s_f.add_dim("nmembers",K); 
			NcDim* nstate = s_f.add_dim("nstate",23);
			//NcVar* meanstate = s_f.add_var("meanstate",ncDouble,nstate,nmembers);
			NcVar* state = s_f.add_var("apri_state",ncDouble,nstate,nmembers);
			s_f.get_var("apri_state")->add_att("dims","3");

			double temp[23 * K]; 
				for(size_t i =0;i<23;++i)
					for(size_t j =0;j<K;++j)
					temp[i*K+j] = x_b(i,j);
			s_f.get_var("apri_state")->put(temp,23,K);	
		 	//free(temp);
			return 1;
		}		
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


void write_input_geos(const char* geos_output, const vector<observation>& obs, const datetime & start_date, const datetime& end_date){
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
            if(kmp_match(sbuf, monStr[end_date.month - 1]) >= 0){
                char restart_buf[32] = {0};
                int m = end_date.monthmax();
                for(int i = 0; i < m; ++i) restart_buf[i] = '0';
                restart_buf[end_date.day - 1] = '3';
                sprintf(sbuf, "Schedule output for %s : %s",monStr[end_date.month - 1], restart_buf); 
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

void write_scaling_factor(const string& f, const vector<double>& v, int XSIZE, int YSIZE, double tau_start, double tau_end, double tau_step){
    //debug("write");
    datablock db(XSIZE,YSIZE,1,1,1,1);
    db.modelres_lat = 180.0 / (YSIZE - 1);
    db.modelres_lon = 360.0 / XSIZE;
    db.tracer = 1;
    str_cpy(db.category, "SCL-FCT");
    
    bpch sf("CTM bin 02","scaling factor for assimilation");

    int k = 0;
    //debug(v);
    while(tau_start + tau_step * k < tau_end){
        db.tau0 = db.tau1 = tau_start + tau_step * k;
        for(int x = 0; x < XSIZE; ++x){
            for(int y = 0; y < YSIZE; ++y){
                db[x][y][0] = v[y * XSIZE + x ];
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
matrix<mT> Mtom(Eigen::MatrixXf& M)
{
	matrix<mT> m(M.rows(),M.cols());
	for(size_t i = 0;i < m.rows();++i)
		for(size_t j = 0;j < m.cols();++j)
			m[i][j] = M(i,j);
	return m;
}

map<int,MatrixXf> EnsembleMembers;
template<typename mT>
MatrixXf MakeNewmatrix(size_t lag, int nlog,const vector< vector<mT> >& COV,const int K)
{
	matrix<mT> cov = (matrix<mT>)COV;	
	matrix<mT> C = chol(cov);
	VectorXf oneMean = VectorXf::Ones(cov.nrow());
	MatrixXf r_M(cov.nrow(), K);
	if(lag == nlog && nlog >= 3)
	{
		//oneMean += EnsembleMembers[lag-2].rowwise().mean() + EnsembleMembers[lag-3].rowwise().mean();
		//oneMean = oneMean / 3.0;
	}
	for(size_t k = 0;k < K;++k)
	{
		vector<double> rands1(11),rands2(12);
		//vector_randn_boost(rands, cov.nrow(), 0, 1, -1, 1);
		//matrix<double> dot_m = dot(C,rands);
		vector_randn_boost(rands1, 11, 1, 0.64, 0.2, 1.8);
		vector_randn_boost(rands2, 12, 1, 0.16, 0.6, 1.4);
		for(size_t j = 0;j < 11;++j)
			r_M.col(k)(j) = rands1[j]; 
		for(size_t j = 0;j < 12;++j)
			r_M.col(k)(11 + j) = rands2[j]; 
	}
	return r_M;
}

template<typename mT>
void get_covariance(vector< vector<mT> >& COV,const int YSIZE)
{
	matrix<double> cov = COV;
	int a_x,a_y,b_x,b_y;
	for(size_t i = 0;i < cov.nrow();++i)
	{
		if(i  < 11)
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
		(i < 11)?cov[i][j] = cov[j][i] = 0.64 * exp(-L/7):cov[i][j] = cov[j][i] = 0.16 * exp(-L/10);

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
    //DebugConfig(DA_config);

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
    int    nlag = 1; // number of lag	
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
    int Class = 23;
    size_t M = Class,M1= XSIZE * YSIZE; // size of state vector;
    vector<double> x_init(M, 1.0);
    string restart_f = DA_config["restart_file"];
    /*
     * generate the first ensemble of state vector x 
     *
     *  DONE by PROCESS-0, then broadcast to MPI_COMM_WORLD
     */
    MatrixXf x_b(M, K);
    matrix<double> scal_x_a(M1,K),cov(Class, Class);
    vector<double> rand_t(M),x_b_p_tre(M1 * K);
    vector<double> rand_s(M);
    double x_b_p[M1 * K];
    int P;
    debug("haha");
    if(mpi_rank == 0) {
    	    //get_covariance(cov,YSIZE);
	    for(size_t lag = 0; lag < nlag; ++lag) 
            {    
			MatrixXf x_b_temp(Class,K);
			x_b_temp = MakeNewmatrix(lag,nlag,cov,K);
	    		debug(x_b_temp);
			for(size_t i = 0; i < Class; ++i) 
	    		for(size_t k = 0; k < K; ++k) 
	    			x_b(i,k) = x_b_temp(i,k);
	    }
	
    		
		for(size_t i = 0; i < K; ++i){
				for(size_t x = 0;x < XSIZE;++x){
				for(size_t y = 0;y < YSIZE;++y){
					int index = XSIZE * y + x,index2 = regions[y][x];       
					scal_x_a[index][i] = x_b(index2,i);
		    		}
				}
	    
	    
	    	for(size_t j=0;j<M1;++j)    
	    		x_b_p[j + M1 * i] = scal_x_a[j][i];

	    	}
    	MPI_Bcast(x_b_p, M1 * K, MPI_DOUBLE, 0, MPI_COMM_WORLD);
   }
   else{
    
    	MPI_Bcast(x_b_p, M1 * K, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	for(size_t j = 0;j< M1 * K;++j)
    	       scal_x_a[j % M1][j / M1] = x_b_p[j];
   }

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
        double   tau_start= DA_start.tau();
        double   tau_end  = DA_end  .tau();
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
        char buf_s[128],cmd_str[128];
        double dx = atof(DA_config["res_lon"].c_str());
        double dy = atof(DA_config["res_lat"].c_str());
		sprintf(cmd_str, "/data1/xubx/GEOS-Chem/ools/generate_prior_co2_dir/apri_changeflux/%s0000.bpch", DA_start.str("YYYYMMDD").c_str());
		printf("reading stations output :%s\n", cmd_str);
		bpch station(cmd_str);
	
                for(size_t j = 0; j < station.size(); ++j){
			int ii = 0;
                        vector<double> rands(150);
			vector_randn_boost(rands, 150, 0, 0.01, -0.1, 0.1);
			if(kmp_match(station[j].category, "IJ-AVG-$") >= 0)
			    for(int i = 0;i<station[j].dim[0];i=i+10)
			    	for(int m = 0;m<station[j].dim[1];m=m+10){
					buf_o.geosI = i + 1;
					buf_o.geosJ = m + 1;
					buf_o.alt   = 1;
					buf_o.tau = station[j].tau0;
					buf_o.mdm = rands[ii] * 1e-6;
					buf_o.value = station[j].data[i][m][0] + rands[ii] * 1e-6;
					obs.push_back(buf_o);
					++ii;
				}
                	}
        /* build up index to speed up searching */
        OBS_INDEX obs_index = build_index(obs);
	/*
         * GENERATE y_b for each ensemble member
         * to prepare the following files:
         *      *aprior_scaling_factor
         *      *input.geos
         */
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
                write_scaling_factor(cmd_str, scal_x_a.get_col(i), XSIZE, YSIZE, datetime(DA_start.str("YYYY-MM-DD 00:00:00")).tau(), tau_end, DA_step / 60.0);
        

                /*
                 * prepare restart file
                 */
              sprintf(cmd_str, "cp %s %s/ensemble-%.3lu/", restart_f.c_str(), DA_dir, i );
                cmd(cmd_str);

                /* 
                 * write input.geos file
                 */
               sprintf(cmd_str, "%s/ensemble-%.3lu/input.geos", DA_dir, i);
                 write_input_geos(cmd_str, obs, DA_start, DA_end);

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
	 	
		char mk_str[128];
            	sprintf(mk_str, "statevector.%s.nc", DA_i.str("YYYYMMDDhhmmss").c_str());
    	    	bool flag = save_apri_vector(mk_str,x_b,K);	
	    	
            for(size_t i = 0; i < K; ++i) {
                char cmd_str[128];//string buffer of command 
                /*
                 * read modeled observation 
                 */
			sprintf(cmd_str, "%s/ensemble-%.3lu/%s0000.bpch", DA_dir, i, DA_start.str("YYYYMMDD").c_str());
                    printf("reading stations output :%s\n", cmd_str);
                    bpch station(cmd_str);
                    for(size_t j = 0; j < station.size(); ++j){
                        if(kmp_match(station[j].category, "IJ-AVG-$") >= 0){
			    	for(int i = 0;i<station[j].dim[0];i=i+10)
			    		for(int m = 0;m<station[j].dim[1];m=m+10){
					buf_mod.geosI = i + 1;
					buf_mod.geosJ = m + 1;
					buf_mod.tau = station[j].tau0;
					buf_mod.value = station[j].data[i][m][0];
					mod.push_back(buf_mod);
					}
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
  	    debug(x_b); 
            MatrixXf X_b(M,K);
	    X_b = x_b;
	    //debug(X_b);
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
                        R_test(i,i) = pow(obs[i].mdm ,2);

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
                          	index_x =int(int(obs[n].t.tau() - tau_start) / 24);
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
			    /*int a_x,a_y,b_x,b_y;
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
			    		for(size_t x = 0;x < YSIZE;++x){
						if(regions_index[i][x] != 0){
							b_x =  x; 
							b_y = regions_index[i][x];
						}
			    		}	
					L = sqrt(pow(a_x - b_x,2) +  pow(a_y - b_y,2));
				 	(i < 11)?Coef(i,0) = exp(-L/7):Coef(i,0) = exp(-L/10);
			   }

			   Coef(index_r) = 1;*/
		            double L;
			    MatrixXf Coef(M,1);
			    for(size_t i = 0; i < M; ++i) 
			    {
				    Coef(i,0) = 0;
				 	L = sqrt(pow(R_c[index_r].first - R_c[i % 23].first,2) +  pow(R_c[index_r].second - R_c[i % 23].second,2));
				 	((i % 23) < 11)?Coef(i,0) = exp(-L/7):Coef(i,0) = exp(-L/10);
			   }
			   Coef(index_r) = 1;
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
					X_b(i,k) = X_b(i,k) - Coef(i,0) * alpha(index_r,index_r) * K_gain(i,n) *(Y_b(n,k));
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
				L = sqrt(pow(R_c[regions[obs[n].geosJ][obs[n].geosI]].first - R_c[regions[obs[m].geosJ][obs[m].geosI]].first,2) +  pow(R_c[regions[obs[n].geosJ][obs[n].geosI]].second - R_c[regions[obs[m].geosJ][obs[m].geosI]].second,2));
			    	(regions[obs[m].geosJ][obs[m].geosI] < 11)?:length_L=7;length_L=10;
				MatrixXf a2(Y_b.row(n)),b2(Y_b.row(m)),Y_b_term(1,1);
			    	Y_b_term = a2 * b2.transpose();
				double fac =(Y_b_term(0,0) / HPHR[n][n]) / ( K - 1);
				y_b_bar(m,0) = y_b_bar(m,0) + exp(-L/length_L) * fac * (y_o(n,0) - y_b_bar(n,0));
				
				for(size_t k = 0;k < K;++k)
				Y_b(m,k) = Y_b(m,k) - exp(-L/length_L) * alpha(index_r, index_r) * fac *(Y_b(n,k));
			   }
			   //debug(y_b_bar);
			   //debug(Y_b);
			   for(size_t m = 0;m < n + 1;++m)
			   {
				L = sqrt(pow(R_c[regions[obs[n].geosJ][obs[n].geosI]].first - R_c[regions[obs[m].geosJ][obs[m].geosI]].first,2) +  pow(R_c[regions[obs[n].geosJ][obs[n].geosI]].second - R_c[regions[obs[m].geosJ][obs[m].geosI]].second,2));
			    	(regions[obs[m].geosJ][obs[m].geosI] < 11)?:length_L=7;length_L=10;
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
            	    vector<double> scal_posterior_x_a(XSIZE * YSIZE * F);
		    for(size_t f = 0;f < F;++f)
		    {
			for(size_t x = 0;x < XSIZE;++x){
		    		for(size_t y = 0;y < YSIZE;++y){
				int index = XSIZE * y + x,index2 = regions[y][x];       
	                   		scal_posterior_x_a[index + XSIZE * YSIZE * f] = x_b_bar(index2 + 23 * f,0);   
		    			for(size_t i = 0; i < K;++i)	
						scal_x_a[index + XSIZE * YSIZE * f][i] = x_b(index2 + 23 * f,i);
				}
            	    	}
		    }
	    	    for(size_t i = 0;i < K;++i)    
	    	    	for(size_t j = 0;j < M1;++j)    
	    			x_b_p[j + M1 * i] = scal_x_a[j][i];
    		    MPI_Bcast(x_b_p, M1 * K, MPI_DOUBLE, 0, MPI_COMM_WORLD);

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
		    write_scaling_factor(buf_str, scal_posterior_x_a, XSIZE, YSIZE, datetime(DA_start.str("YYYY-MM-DD 00:00:00")).tau(), tau_end, DA_step / 60.0);
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
		    write_scaling_factor(cmd_str, scal_posterior_x_a, XSIZE, YSIZE, datetime(DA_start.str("YYYY-MM-DD 00:00:00")).tau(), tau_end, DA_step / 60.0);
			
		    /*
		     * prepare restart file
		     */
		    sprintf(cmd_str, "cp %s %s/posterior/", restart_f.c_str(), DA_dir);
		    cmd(cmd_str);
		    /* 
		     * write input.geos file
		     */
		    sprintf(cmd_str, "%s/posterior/input.geos", DA_dir);
		    write_input_geos(cmd_str, obs, DA_start, DA_end);
		    /* 
		     * run geos model 
		     */
		    sprintf(cmd_str, "(cd %s/posterior; %s )",DA_dir, DA_config["geos"].c_str());
		    cmd(cmd_str, true);
		}	
   		else{
    
    		MPI_Bcast(x_b_p, M1 * K, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		for(size_t j = 0;j< M1 * K;++j)
    	       		scal_x_a[j % M1][j / M1] = x_b_p[j];
   		}


	    }// END OF DA CYCLE
		
	    MPI_Finalize();
	    return 0;
	}

			    


		
