#include <iostream>
#include <algorithm>
#include <netcdfcpp.h>
#include <Eigen/Dense>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include "bpch.h"
#include <Eigen/Dense>
#include "xg_math_vector.h"
using namespace Eigen;
using namespace std;

	void readvector(const char * file_r)
	{
		//statevector<Tm> r_tmp(5,247,80);
		NcFile s_f(file_r, NcFile::ReadOnly);
		if(!s_f.is_valid())
		{
			cout<<"couldn't open file"<<endl;
			//return 0;
		}
		double temp[5 * 247 * 80];
	        s_f.get_var("apri_state")->get(temp,5,80,247);	
	  	//debug(temp);	
		for(size_t i = 0;i < 5*80*247;++i)
			cout<<temp[i]<<" ";
			
	}	

int main(int argc, char *argv[])
{
	//statevector<double> test(5,229,10,0);
	//test.Initialization();
	//test.Propagate();
	//bool r = test.save_apri_vector(argv[1]);
	readvector(argv[1]);
	/*    ReadConfig(argv[1], DA_config);
    DebugConfig(DA_config);
	for(size_t i=0;i<5;++i)
		for(size_t j=0;j<229;++j)
			for(size_t z=0;z<80;++z)
		cout<<test(i,j,z)<<endl;
*/	return 0;		
}

