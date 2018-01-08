/* *************************
 * Author: xbx1992
 * Created Time:  Tue 4  July 2015
 * LastModified:  Tue 11 Dec 2014 04:00:27 PM CST
 * C File Name:  state.h
 * ************************/
#ifndef STATE_H
#  define STATE_H
#include <vector>
#include <deque>
#include "configclass.h"
#include <xg_math_matrix.h>
#include <netcdfcpp.h>
#include <Eigen/Dense>

using namespace Eigen;
using Eigen::MatrixXf;


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
		
		void Initialization(daconfig dc, map<int,vector<pair<int,int> > >& regions_index)
		{
				
			/*MatrixXf cov = get_covariance(dc,regions_index);	
			LLT<MatrixXf> llt(cov);
			MatrixXf L = llt.matrixL();	
			matrix<double> C(L.rows(),L.cols());
			for(size_t i = 0;i < C.nrow();++i)
				for(size_t j = 0;j < C.ncol();++j)
					C[i][j] = L(i,j);*/
			vector<double> oneMean(_nrow,1.0);
			for(size_t lag=0;lag< _nlag;++lag)
			for(size_t k = 0;k < _ncol;++k)
			{
				vector<Tm> rands1(11),rands2(11);
				//vector<Tm> rands1(dc.Class);
				//vector_randn_boost(rands, cov.nrow(), 0, 1, -1, 1);
				vector_randn_boost(rands1, 11, 0, 0.25, 0,0);
				vector_randn_boost(rands2, 11, 0, 0.01, 0, 0);
				//matrix<Tm> ans = dot(C,rands1);
			        (*this)[lag][0][k] = 1;	
				for(size_t j = 0;j < 11;++j)
					(*this)[lag][j + 1][k] = rands1[j] + 1; 
				for(size_t j = 0;j < 11;++j)
					(*this)[lag][j + 12][k] = rands2[j] + 1;
			}
		}
		MatrixXf get_covariance(daconfig dc,map<int,vector<pair<int,int> > >& regions_index)
		{
			/*define regions_index*/
		 	int index[] = {0,1,14,29,37,46,58,67,81,97,107,115,131,230};

			MatrixXf cov(dc.Class,dc.Class);
			int a_x,a_y,b_x,b_y;
			for(size_t ct = 1;ct < 13; ++ct)
			{
				if(ct != 1 && ct!=7 && ct != 8 )
				{
					int cov_i;
				       (ct==12)?cov_i = 0.16:cov_i = 0.64;	
					for(size_t i = index[ct];i < index[ct+1];++i)	
					{
			    			int index_v = regions_index[i].size() / 2;
			    			a_x =  regions_index[i][index_v].first; 
			    			a_y =  regions_index[i][index_v].second;
						for(size_t j = index[ct];j < i;++j)
						{
			    				index_v = regions_index[j].size() / 2;
							b_x =  regions_index[j][index_v].first;; 
							b_y = regions_index[j][index_v].second;
							int L = sqrt(pow(a_x - b_x,2) +  pow(a_y - b_y,2));
							cov(i,j) = cov(j,i) = cov_i * exp(-L/10);

						}
					}	
				}
				for(size_t i = 0;i < dc.Class;++i)
				{
				if(i  < 131 && i!=0)
					cov(i,i) = 0.64;
				else
					cov(i,i) = 0.16;
				}
			}
			for(size_t i = 0;i < 230;++i)	
			{
				for(size_t j = 0;j < 230;++j)
			{
				cout<<cov(i,j);
			}	
			cout<<endl;
			}
			return cov;
		}	
		bool Propagate(daconfig dc)
		{
			bool rt = false;
	    		MatrixXf X_b = (*this).transfer(); 
	    		
			MatrixXf x_b_bar = X_b.rowwise().mean(); 
	    		debug(x_b_bar);	
			for(size_t i = 0; i < X_b.cols(); ++i)
            		{
				
				for(size_t j = 0; j < x_b_bar.rows();++j){
					X_b(j,i) = X_b(j,i) - x_b_bar(j,0); 
				} 
	    		}  
			
			matrix<Tm> temp(_nrow, _ncol);
			
					
			for(size_t k = 0;k < _ncol;++k)
			{
				vector<Tm> rands1(11),rands2(11);
				vector_randn_boost(rands1, 11, 0, 0.25, 0, 0);
				vector_randn_boost(rands2, 11, 0, 0.01, 0, 0);
				(*this)[0][0][k] = 1; 

				for(size_t j = 0;j < 11;++j)
					(*this)[0][j + 1][k] = rands1[j] / 2 + X_b(j + 1,0) / 2 + (x_b_bar(j + 1,0) + 1) / 2; 
				//for(size_t j = 0;j < 11;++j)
				//	(*this)[0][j + 1][k] = X_b(j + 1,0) + x_b_bar(j + 1,0); 
				for(size_t j = 0;j < 11;++j)
					(*this)[0][j + 12][k] = rands2[j] / 2 + X_b(j + 1,0) / 2 + (x_b_bar(j + 12,0) + 1) /2;
			}
			rt = true;	
			return rt;		
		}
            	vector<Tm> get_col(size_t c){
			vector<Tm> ans;
			for(size_t lag=0;lag< _nlag;++lag)
				for(size_t i = 0;i < _nrow;++i)
				ans.push_back((*this)[lag][i][c]);
			return ans;
				
		}
		
		Tm *expand()
		{
			Tm *ans = (Tm*)malloc(_ncol * _nlag * _nrow * sizeof(Tm));
			for(size_t k = 0;k < _ncol;++k)
				for(size_t lag=0;lag< _nlag;++lag)
					for(size_t i = 0;i < _nrow;++i)
					ans[i + lag * _nrow + k * _nlag * _nrow] = ((*this)[lag][i][k]);
			//ans = NULL;
			return ans;	
		}
		MatrixXf transfer()
		{
			MatrixXf temp(_nrow * _nlag , _ncol);
			for(size_t lag=0;lag< _nlag;++lag)
				for(size_t i = 0;i < _ncol;++i)
					for(size_t j = 0;j < _nrow;++j)
					temp(j + _nrow * lag, i) = ((*this)[lag][j][i]);
			return temp; 
		}	
	       	
		int size()
		{
			return this->size();
		}
		Tm& operator()(size_t i, size_t j, size_t z){
			return (*this)[i][j][z]; 
		}
	
		statevector<Tm> operator=(deque<matrix<Tm> > s)
		{
			assert(s.size() == this->size());
			for(size_t lag = 0;lag < s.size();++lag)
				(*this)[lag] = s[lag];
			return *this;
		}	
		bool save_apri_vector(const char * file_s)
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



			NcDim* nlag = s_f.add_dim("nlags",_nlag); 
			NcDim* nmembers = s_f.add_dim("nmembers",_ncol); 
			NcDim* nstate = s_f.add_dim("nstate",_nrow);
			//NcVar* meanstate = s_f.add_var("meanstate",ncDouble,nstate,nmembers);
			NcVar* state = s_f.add_var("apri_state",ncDouble,nlag,nstate,nmembers);
			s_f.get_var("apri_state")->add_att("dims","3");

			double temp[_nlag * _nrow * _ncol]; 
			for(size_t lag =0;lag<_nlag;++lag)
				for(size_t i =0;i<_nrow;++i)
					for(size_t j =0;j<_ncol;++j)
					temp[lag * _nrow * _ncol+i*_ncol+j] = (*this)[lag][i][j];
			s_f.get_var("apri_state")->put(temp,_nlag,_nrow,_ncol);	
		 	//free(temp);
			return 1;		
		}

		bool save_post_vector(const char * file_s)
		{

			NcFile s_f(file_s, NcFile::Write);
			if(!s_f.is_valid())
			{
				cout<<"couldn't open file"<<endl;
				return 0;
			}
			NcDim* nlag = s_f.get_dim("nlags"); 
			NcDim* nmembers = s_f.get_dim("nmembers"); 
			NcDim* nstate = s_f.get_dim("nstate");

			NcVar* state = s_f.add_var("post_state",ncDouble,nlag,nstate,nmembers);
			s_f.get_var("post_state")->add_att("dims","3");

			double temp2[_nlag * _nrow * _ncol]; 
			for(size_t lag =0;lag<_nlag;++lag)
				for(size_t i =0;i<_nrow;++i)
					for(size_t j =0;j<_ncol;++j)
					temp2[lag * _nrow * _ncol+i*_ncol+j] = (*this)[lag][i][j];
			s_f.get_var("post_state")->put(temp2,_nlag,_nrow,_ncol);	
		 	//free(temp);
		 	return 1;		
		}
		
		statevector<double> state2grid(REGIONS_MAP& regions, daconfig dc)
		{
			statevector<double> scal_x_a(dc.nlag,dc.F * dc.XSIZE * dc.YSIZE, dc.K);
	    		for(size_t lag = 0; lag < dc.nlag; ++lag) 
            		{    
    				for(size_t i = 0; i < dc.K; ++i)
				{	
					for(size_t x = 0;x < dc.XSIZE;++x){
					for(size_t y = 0;y < dc.YSIZE;++y){
						int index = dc.XSIZE * y + x,index2 = regions[y][x];       
						scal_x_a(lag,index,i) = (*this)[lag][index2][i];
		    			}
					}
				}
	    		}
		       	return scal_x_a;	
		}

	};

//template<typename Tm>
	statevector<double> reverse_transfer(const MatrixXf& x_b, int _nlag, int _nrow)
	{
		statevector<double> temp(_nlag, _nrow, x_b.cols());
		for(size_t i = 0;i < x_b.cols();++i)
			for(size_t j = 0;j < (_nrow * _nlag);++j)
				temp(j / _nrow,j % _nrow,i) = x_b(j,i);
		return temp; 
	}	
	
	void reverse_expand(statevector<double>& temp,double *x_b_p, int _nlag, int _nrow, int _ncol)
	{
		
		//statevector<double> temp(_nlag, _nrow, _ncol); 
		for(size_t j = 0;j<  _nlag * _nrow * _ncol;++j)
		{
    	       		temp((j % (_nrow * _nlag)) / _nrow, (j % (_nrow * _nlag)) % _nrow, j / (_nrow * _nlag)) = x_b_p[j];
		}
		debug("lala");
		//return temp;
	}
template<typename Tm>
	ostream& operator<<(ostream& o, const statevector<Tm>& st)
	{
		o << "statevector("<<st.size() * st[0].nrow() * st[0].ncol() <<"){\n";
		for(size_t i = 0;i < st.size();++i)
		{
			o << "statevector("<< i <<"){\n";
			for(size_t z = 0;z < st[0].nrow();++z)
			{
			for(size_t j = 0;j < st[0].ncol();++j)
			{
				o << st[i][z][j];
				
			}
				o << "\n";
			}
			o << "}";
		}
		return o<<"}";
	}
//template<typename Tm>
	statevector<double> readvector(const char * file_r,daconfig dc)
	{
		statevector<double> r_tmp(dc.nlag,dc.Class,dc.K);
		NcFile s_f(file_r, NcFile::ReadOnly);
		if(!s_f.is_valid())
		{
			cout<<"couldn't open file"<<endl;
		//	return 0;
		}
		double temp[dc.nlag * dc.Class * dc.K];
	        s_f.get_var("apri_state")->get(temp,dc.nlag,dc.Class,dc.K);	
		for(size_t i =0;i<dc.nlag * dc.Class * dc.K;++i)
			r_tmp(i / (dc.K * dc.Class),(i % (dc.K * dc.Class)) / (dc.K),(i % (dc.K * dc.Class)) % (dc.K)) = temp[i];
		return r_tmp;
			
	}	

	vector<double> poststate2grid(REGIONS_MAP& regions, MatrixXf& x_b_bar, daconfig dc)
	{
			
            	    vector<double> scal_posterior_x_a(dc.XSIZE * dc.YSIZE);

		    for(size_t x = 0;x < dc.XSIZE;++x){
		    		for(size_t y = 0;y < dc.YSIZE;++y){
				int index = dc.XSIZE * y + x,index2 = regions[y][x];       
	                   		scal_posterior_x_a[index] = x_b_bar(index2,0);   
				}
            	    	}
		    debug(scal_posterior_x_a);
		    return scal_posterior_x_a;
	}
	/*template<typename mT>
	MatrixXf m2M(matrix<mT>& m)
	{
		MatrixXf M(m.nrow(),m.ncol());
		for(size_t i = 0;i < m.nrow();++i)
			for(size_t j = 0;j < m.ncol();++j)
				M(i,j) = m[i][j];
		return M;
	}*/	
	matrix<double> M2m(MatrixXf& M)
	{
		matrix<double> m(M.rows(),M.cols());
		for(size_t i = 0;i < m.nrow();++i)
			for(size_t j = 0;j < m.ncol();++j)
				m[i][j] = M(i,j);
	
		return m;
	}

#endif /* ifndef STATE_H */
