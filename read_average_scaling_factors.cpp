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
using namespace std;
using boost::lexical_cast;

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
	matrix<double> readpostvector(const char * file_r,int Class, int K)
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
		//debug(oneMean);
		return oneMean;
	}
int main(int argc, char *argv[])
{	
	
    	       int  Class = lexical_cast<int>(argv[2]);
    	       int  K = lexical_cast<int>(argv[3]);

    	       matrix<double> x_b(Class, 1);
               char cmd_str[128];//string buffer of command 
               sprintf(cmd_str, "%s", argv[1]);
	       x_b = readpostvector(cmd_str,Class, K); 
	       double sum;
	       debug(x_b);
		for(size_t i = 0;i < Class;++i)
	       {
		if(i <= 131)
	       		sum += x_b[i][0];
			
			
	       	if(i == (Class - 1))
		{
			cout<<sum / 131.0<<endl;	
			sum = 0;
		}
			
	       }
	       sum = 0;
	       cout<<endl;
	       cout<<endl;
	       cout<<endl;
	       for(size_t i = 0;i < Class;++i)
	       {
		if(i > 131)
	       		sum += x_b[i][0];
	       	if((i == 0 && i > 0) || i == (Class - 1))
		{
			cout<<sum / 99.0<<endl;	
			sum = 0;
		}
			
	       }
	       return 0;		
}

