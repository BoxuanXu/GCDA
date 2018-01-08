#include <iostream>
#include <algorithm>
#include <netcdfcpp.h>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include "bpch.h"
#include "xg_math_vector.h"
#include <deque>
#include "xg_debug.h"
#include <boost/lexical_cast.hpp>
#include <Eigen/Dense>
#include <random>
#include <boost/random.hpp>
#include <boost/random/random_device.hpp>

using namespace std;
using boost::lexical_cast;
using namespace Eigen;
using Eigen::MatrixXf;
boost::random::mt19937 gen;

template<typename Tm>
class statevector: public deque<matrix<Tm> >{ 
	     private:
		size_t _nlag,_ncol,_nrow; 
	     public:
		statevector(){}
		statevector(size_t nlag, size_t nrow, size_t nensemble,Tm v=0):_nlag(nlag),_ncol(nensemble),_nrow(nrow)
		{
			matrix<Tm> temp(nrow,nensemble,v);
			for(size_t lag=0;lag<nlag;++lag)
			{
			   (*this).push_back(temp);
			}
		}

		statevector(deque<matrix<Tm> > sv)
		{
			for(size_t lag = 0;lag < sv.size();++lag)
				(*this).push_back(sv[lag]);
					
		}
		Tm& operator()(size_t i, size_t j, size_t z){
			return (*this)[i][j][z]; 
		}
};
	void readpostvector(const char * file_r,int Class, int K)
	{
		NcFile s_f(file_r, NcFile::ReadOnly);
		if(!s_f.is_valid())
		{
			cout<<"couldn't open file"<<endl;
		//	return 0;
		}
		double temp[Class * K];
		matrix<double> r_mean(Class ,K);
	        s_f.get_var("post_state")->get(temp, 1, Class, K);	
		for(size_t i = 0;i < Class * K;++i)
			r_mean[i / K][ i % K] = temp[i];
		matrix<double> oneMean = rowmean(r_mean);	
		for(size_t i = 0;i < r_mean.nrow(); ++i)
            		for(size_t j = 0; j < r_mean.ncol(); ++j)
				r_mean[i][j] = r_mean[i][j] - oneMean[i][0]; 
		matrix<double> P_a = (r_mean * t(r_mean)) / (K - 1);
		debug(P_a);
	
		oneMean = rowmean(r_mean);	
		MatrixXf P_a_grid(Class ,1);
		for(size_t i = 0;i < r_mean.nrow(); ++i)
		{
            		for(size_t j = 0; j < r_mean.ncol(); ++j)
			P_a_grid(i,0) += (r_mean[i][j] - oneMean[i][0]) * (r_mean[i][j] - oneMean[i][0]);
			P_a_grid(i,0) = P_a_grid(i,0) / r_mean.ncol();
		}
			debug(P_a_grid);
		//debug(oneMean);
		//return oneMean;
	}
int main(int argc, char *argv[])
{	
	
    	       int  Class = lexical_cast<int>(argv[1]);
    	       double  K = lexical_cast<double>(argv[2]);
    	       double  a = lexical_cast<double>(argv[3]);
    	       double  b = lexical_cast<double>(argv[4]);
    	       double  c = lexical_cast<double>(argv[5]);
    	       double  d = lexical_cast<double>(argv[6]);
    	       double  e = lexical_cast<double>(argv[7]);

    	       /*matrix<double> x_b(Class, 1);
               char cmd_str[128];//string buffer of command 
               sprintf(cmd_str, "%s", argv[1]);
	        readpostvector(cmd_str,Class, K);*/ 
				
	       vector<double> rands1(Class);
	       vector_randn_boost(rands1, Class, K, a, b, c,d,e);
	       
	       // rands1 = vector_randn(Class, K, a);  
            	double P_a_grid=0,mean=0;
		for(int j = 0; j < rands1.size(); ++j)
		{
			mean += rands1[j];
		}
		mean = mean / rands1.size();
		cout<<"mean:"<<mean<<endl;
		for(int j = 0; j < rands1.size(); ++j)
	       		P_a_grid  += (rands1[j] - mean) * (rands1[j] - mean);
	        P_a_grid = P_a_grid / rands1.size();
		cout<<"P_a:"<<P_a_grid<<endl;;
	       	
		return 0;		
}

