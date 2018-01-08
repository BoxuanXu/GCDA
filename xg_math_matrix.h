/*
 * $Id$
 */
#ifndef XG_MATH_MATRIX_H
#  define XG_MATH_MATRIX_H

#include <vector>
#include <fstream>
#include <cassert>

#include "xg_math_vector.h"
using namespace std;

template<typename mT>
	class col_ref: public vector< mT*>{
		public:
			using vector< mT* >::at;
			using vector< mT* >::push_back;
			col_ref(){}
			col_ref(vector<mT>& col){
			for(size_t i = 0;i < col.size();++i)
				this->push_back(&col[i]);		
		     	}
			col_ref operator=(const mT& v){
				assert(v.size() = size());
				for(size_t i =0 ;i < size();++i)
					*at(i) = v[i];
				return *this;
			}
			col_ref operator=(const vector<mT>& v){
				assert(v.size() = size());
				for(size_t i = 0;i < v.size();++i)
				{
					*at(i) = v[i];
					cout<<"come here"<<endl;
				}
				return *this;
			}
			col_ref operator=(col_ref& v){
				assert(v.size() == size());
				for(size_t i = 0;i < size(); ++i)
					*at(i) = v[i];
				return *this;
			}
			size_t size(){
				return this -> vector<mT*>::size();
			}
			mT& operator[](size_t i ){
				return *at(i);
			} 
			mT& operator()(size_t i){
				return *at(i);
			}
			col_ref& operator+=(const vector<mT>& col){
				assert(this->size() = col.size());
				for(size_t i = 0; i < this->size();++i)
					*at(i) += col[i];
				return *this;
			}
			col_ref& operator-=(const vector<mT>& col){
                assert(this->size() == col.size());
                for(size_t i = 0; i < this->size(); ++i)
                    *at(i) -= col[i];
                return *this;
            }
            col_ref& operator*=(const vector<mT>& col){
                assert(this->size() == col.size());
                for(size_t i = 0; i < this->size(); ++i)
                    *at(i) *= col[i];
                return *this;
            }
            col_ref& operator/=(const vector<mT>& col){
                assert(this->size() == col.size());
                for(int i = 0; i < this->size(); ++i)
                    *at(i) /= col[i];
                return *this;
            }
		col_ref& operator+=(mT v){
			for(int i = 0;i<this->size();++i)
				*(at(i)) += v;
			return *this;
		}
		col_ref& operator-=(mT v){
                for(int i = 0; i < this->size(); ++i)
                    *(at(i)) -= v;
                return *this;
            }
            col_ref& operator*=(mT v){
                for(int i = 0; i < this->size(); ++i)
                    *(at(i)) *= v;
                return *this;
            }
            col_ref& operator/=(mT v){
                for(int i = 0; i < this->size(); ++i)
                    *(at(i)) /= v;
                return *this;
            }
};

template<typename T>
	col_ref<T> select_ref(vector<T>& v, const vector<size_t>& idx){
		col_ref<T> ans;
		for(size_t i=0;i<idx.size();++i){
			ans.push_back(&(v[idx[i]]));
		}
		return ans;
	}
template<typename mT>
    ostream& operator<<(ostream& o, const col_ref<mT>& v){
        o << "col("<< v.size() <<"X1){\n";
        for(typename vector< mT* >::const_iterator _v = v.begin(); _v != v.end(); ++_v){
            o << *_v << "\n";
        }
        return o<<"}";
    }
template<typename mT>
    class matrix: public vector< vector< mT> >{
        public:
            matrix(){}
            matrix(size_t _nrow, size_t _ncol, mT v = 0 , int flag = 1):
                vector< vector<mT> > (_nrow)
        {
	   if(flag == 1)
	   {
#if __cplusplus >= 201103L
            for(auto& _row: (*this))
                _row.resize(_ncol, v);
#else
            for(typename vector< vector< mT> >::iterator _row =  (*this).begin(); _row !=  (*this).end(); ++_row)
                _row->resize(_ncol, v);
#endif
	   }
           else
	  {
		vector<mT> col_test(13,v);
		for(size_t i = 0;i < 23 ;++i)
			(*this)[i] = col_test;
	  }	
        }
            matrix(const vector<mT> vec, size_t _ncol = 1):
                vector< vector<mT> > (vec.size() / _ncol +(vec.size() % _ncol != 0) ){
                    for(size_t i = 0; i < size(); ++i){
                        at(i).resize(_ncol, 0);
                        for(size_t j = 0; j < at(i).size();++j){
                            at(i)[j] = vec[i * _ncol + j];
                        }
                    }
                }
            matrix(const vector<vector<mT> >& _V):
                vector<vector<mT> >(_V){
                    *this = t(*this);
                }
            using vector< vector< mT> >::size;
            using vector< vector< mT> >::at;
            size_t nrow()const{
                return size();
            }
            size_t ncol()const{
                if(size() > 0)
                    return at(0).size();
                else 
                    return 0;
            }
            vector<mT> col(size_t i){
                vector<mT> ans;
#if __cplusplus >= 201103L
                for(auto _row :(*this)){
                    ans.push_back(_row[i]);
                }
#else
                for(typename vector< vector< mT> >::const_iterator _row =  this->begin(); _row !=  this->end(); ++_row){
                    ans.push_back((*_row)[i]);
                }
#endif
                return ans;
            }
            vector<mT> row(size_t i){
                return at(i);
            }
	
	    vector<mT> get_row(size_t r){
	        return at(r);
	    }
	    vector<mT> get_col(size_t c){
		vector<mT> ans;
		size_t _nrow=this->nrow();
		for(size_t j = 0;j < _nrow; ++j){
			ans.push_back(at(j)[c]);
		}
	        return ans;
	    }
	    matrix<mT> get_rows(const vector<size_t> select){
		matrix<mT> ans(select.size(),ncol());
		for(size_t i = 0;i < select.size();++i)
		{
			//ans.row(i) = at(select[i]);
			for(size_t j = 0;j < ncol();++j)
				ans[i][j] = at(select[i])[j];
		}
		return ans;
 		//return_matrix.setrow()
 	    }
 	    template<typename aa>
	    matrix& save_array(const aa savearr){
		for(size_t i=0;i<sizeof(savearr);++i)
		{
			at(i/at(0).size())[i%at(0).size()] = savearr[i];
		}	
		return *this;
	    }
	    template<typename aa>
  	    matrix& load_array(aa loadarr){
	    	for(size_t i=0;i<size()*at(0).size();++i)
		{
			loadarr[i] = at(i/at(0).size())[i%at(0).size()];
		}
	    }
            void resize(size_t _nrow, size_t _ncol){
                vector< vector<mT> > :: resize(_nrow);
#if __cplusplus >= 201103L
                for(auto& _row: (*this))
                    _row.resize(_ncol);
#else
                for(typename vector< vector< mT> >::iterator _row =  (*this).begin(); _row !=  (*this).end(); ++_row )
                    _row->resize(_ncol);
#endif
            }
            void setcol(size_t j, const vector<mT>& _col){
                assert(_col.size() == ncol());
                for(size_t i = 0; i < nrow(); ++i)
                    at(i)[j] = _col[j];
            }
            void setrow(size_t i, const vector<mT>& _row){
                assert(_row.size() == nrow());
                at(i) = _row;
            }
            matrix& swapcol(size_t i, size_t j){
                //cout << "swap col " << i << ", " << j << '\n';
                mT _cell;
#if __cplusplus >= 201103L
                for(auto& _row :(*this)){
                    _cell = _row[i];
                    _row[i] = _row[j];
                    _row[j] = _cell;
                }
#else
                for(typename vector< vector< mT> >::iterator _row =  this->begin(); _row !=  this->end(); ++_row ){
                    _cell = (*_row)[i];
                    (*_row)[i] = (*_row)[j];
                    (*_row)[j] = _cell;
                }
#endif
                return *this;
            }
            matrix& swaprow(size_t i, size_t j){
                //cout << "swap row " << i << ", " << j << '\n';
                vector<mT> _row = at(i);
                at(i) = at(j);
                at(j) = _row;
                return *this;
            }
            mT operator()(size_t i, size_t j)const{
                return at(i)[j];
            }
	    mT& operator()(size_t i, size_t j){
		return at(i)[j];
	    }
	    void load_array(const mT* a){
		size_t nc = ncol(),nr = nrow();
		for(size_t i = 0;i < nr;++i){
			for(size_t j = 0; j < nc;++j){
				at(i)[j] = a[i * nc + j];
			}
		}
	    }
	    void save_array(mT* a){
	     	size_t nc =ncol(),nr = nrow();
		for(size_t i = 0; i < nr; ++i){
			for(size_t j = 0; j < nc ; ++j)
			{
				a[i * nc +j ] =  at(i)[j];
			}
		}
	    }
	    void saveA(const char* f){
	    	fstream file(f,ios::out);
		size_t nc= ncol(),nr = nrow();
		file << nr << " " << nc << endl;
		for(size_t i =0 ;i < nr; ++i){
			for(size_t j = 0; j < nc; ++j){
					file << at(i)[j] << " ";
			}
			file<< endl;
		} 
		file.close();
	    }
	    void loadA(const char* f){
		fstream file(f,ios::in);
		size_t nc ,nr;
		mT buf_v;
		file >> nr >> nc;
		resize(nr,nc);
		for(size_t i=0;i<nr;++i){
			for(size_t j = 0;j < nc; ++j){
				file >> at(i)[j];
			}
		}
		file.close();
	   }
	   void saveB(const char* f){
                fstream binary_file(f,ios::out|ios::binary); 
                size_t nc, nr;
                mT buf_v;
                char *bufs;
                nr = nrow(); binary_file.write(reinterpret_cast<char *>(&nr),sizeof(size_t));
                nc = ncol(); binary_file.write(reinterpret_cast<char *>(&nc),sizeof(size_t));
                for(size_t i = 0; i < nr; ++i){
                    for(size_t j = 0; j < nc; ++j){
                        buf_v = at(i)[j];
                        binary_file.write(reinterpret_cast<char *>(&buf_v),sizeof(mT));
                    }
                }
                
	}
	 void loadB(const char* f){
                fstream binary_file(f,ios::binary|ios::in);
                size_t nc, nr;
                mT buf_v;
                binary_file.read(reinterpret_cast<char *>(&nr),sizeof(size_t));
                binary_file.read(reinterpret_cast<char *>(&nc),sizeof(size_t));
                cout << "nr = " << nr << ", nc = " << nc << endl;
                resize(nr, nc);
                for(size_t i = 0; i < nr; ++i){
                    for(size_t j = 0; j < nc; ++j){
                        binary_file.read(reinterpret_cast<char *>(&buf_v),sizeof(mT));
                        at(i)[j] = buf_v;
                    }
                }
                binary_file.close();
            }

	   
    };
template<typename mT>
    matrix<mT> operator+(const matrix<mT>& a,const matrix<mT>& b){
        assert(a.nrow() == b.nrow());
        assert(a.ncol() == b.ncol());
        matrix<mT> ans(a);
        for(size_t i = 0; i < ans.nrow(); ++i){
            for(size_t j = 0; j < ans.ncol(); ++j){
                ans[i][j] = a[i][j] + b[i][j];
            }
        }
        return ans;
    }
template<typename mT>
    matrix<mT> operator-(const matrix<mT>& a,const matrix<mT>& b){
        assert(a.nrow() == b.nrow());
        assert(a.ncol() == b.ncol());
        matrix<mT> ans(a);
        for(size_t i = 0; i < ans.nrow(); ++i){
            for(size_t j = 0; j < ans.ncol(); ++j){
                ans[i][j] = a[i][j] - b[i][j];
            }
        }
        return ans;
    }

template<typename mT>
   matrix<mT> operator+(const matrix<mT>& a,const vector<mT>& b){
	matrix<mT> ans(a.nrow(),a.ncol());
   	for(size_t i=0;i< a.nrow();++i)
		for(size_t j=0;j < a.ncol();j++)
            	    ans[i][j] = a[i][j]+ b[j];
	return ans;
   }

template<typename mT>
    matrix<mT> operator*(const matrix<mT>& a,const matrix<mT>& b){
        assert(a.ncol() == b.nrow());
        matrix<mT> ans(a.nrow(), b.ncol());
        for(size_t i = 0; i < ans.nrow(); ++i){
            for(size_t j = 0; j < ans.ncol(); ++j){
                ans[i][j] = 0;
                for(size_t k = 0; k < a.ncol(); ++k)
                    ans[i][j] += a(i,k) * b(k,j);
            }
        }
        return ans;
    }
template<typename mT>
   matrix<mT> operator*(const matrix<mT>& a,const mT& b){
	matrix<mT> ans(a);
	for(size_t i=0;i<ans.nrow();++i){
		for(size_t j=0;j < ans.ncol(); ++j){
			ans[i][j] = ans[i][j] * b;
		
		}
	}
	return ans;
}
template<typename T1,typename T2>
   matrix<T1> operator*(const matrix<T1>& a,const T2& b){
	matrix<T1> ans(a);
	for(size_t i=0;i<ans.nrow();++i){
		for(size_t j =0;j < ans.ncol();++j){
			ans[i][j] = ans[i][j] * (T1)b;
		}
	}
}
template<typename mT>
    matrix<mT> dot(const matrix<mT>& a,const vector<mT>& b){
	matrix<mT> ans(a.nrow(),1);
	assert(a.ncol() == b.size());
	for(size_t i = 0;i < ans.nrow();++i)
	{
		ans[i][0] = 0;
		for(size_t j = 0;j < a.ncol();++j)
		   {
			ans[i][0] += a[i][j] * b[j];
		   }

	}	
	return ans;
    }

template<typename mT>
    matrix<mT> operator/(const matrix<mT>& a,const matrix<mT>& b){
        assert(a.nrow() == b.nrow());
        assert(a.ncol() == b.ncol());
        matrix<mT> ans(a);
        for(size_t i = 0; i < ans.nrow(); ++i){
            for(size_t j = 0; j < ans.ncol(); ++j){
                ans[i][j] = a[i][j] / b[i][j];
            }
        }
        return ans;
    }

template<typename T1,typename T2>
    matrix<T1> operator/(const matrix<T1>& a,const T2& b){
        //assert(a.nrow() == b.nrow());
        //assert(a.ncol() == b.ncol());
        matrix<T1> ans(a);
        for(size_t i = 0; i < ans.nrow(); ++i){
            for(size_t j = 0; j < ans.ncol(); ++j){
                ans[i][j] = a[i][j] / b;
            }
        }
        return ans;
    }

/*
 * 矩阵转置
 */
template<typename mT>
    matrix<mT> t(const matrix<mT> & m){
        matrix<mT> ans(m.ncol(), m.nrow());
        for(size_t i = 0; i < ans.nrow(); ++i){
            for(size_t j = 0; j < ans.ncol(); ++j){
                ans[i][j] = m(j,i);
            }
        }
        return ans;
    }
 /*
 * 矩阵 求逆
 */
template<typename mT>
    matrix<mT> operator~(const matrix<mT> & m){
        matrix<mT> ans(m);
        //cout << " ans = " << ans << "\n";
        vector<pair<pair< size_t, size_t> , bool> > swaplist;
        // start swap
        for(size_t _t = 0; _t < ans.nrow(); ++_t){
           //cout << "_t = " << _t << '\n';
            // select main element
            size_t maxi = _t, maxj = _t;
            mT maxv = ::abs(ans(_t,_t));
            for(size_t i = _t; i < ans.nrow(); ++i){
                for(size_t j = _t; j < ans.ncol(); ++j){
                    if( ::abs(ans(i, j)) > maxv){
                        maxi = i; maxj = j; maxv = ::abs(ans(i,j));
                    }
                }
            }
            if(maxi != _t) {
                ans.swaprow(_t, maxi);
                swaplist.push_back(make_pair(make_pair(_t,maxi),true));
            }
            if(maxj != _t){
                ans.swapcol(_t, maxj);
                swaplist.push_back(make_pair(make_pair(_t,maxj),false));
            }
            //cout << "before ans = " << ans << "\n";
            ans[_t][_t] = 1 / ans(_t, _t);
            for(size_t j = 0; j < ans.nrow(); ++j)if(j != _t)
                ans[_t][j] *= ans(_t, _t);
	    for(size_t i = 0;i<ans.nrow();++i)if(i!= _t)
               for(size_t j = 0; j < ans.ncol(); ++j)if(j != _t)
                    ans[i][j] -= ans(i, _t) * ans(_t, j);
            for(size_t i = 0; i < ans.nrow(); ++i)if(i != _t)
                ans[i][_t] = - ans(i, _t) * ans(_t, _t);
            //cout << "after: ans = " << ans << "\n";
        }
        reverse(swaplist.begin(), swaplist.end());
#if __cplusplus >= 201103L
        for(auto _l :swaplist ) {
            if(_l.second) ans.swapcol(_l.first.first, _l.first.second);
            else ans.swaprow(_l.first.first, _l.first.second);
        }
#else
        for(vector<pair<pair< size_t, size_t> , bool> >::iterator _l =  swaplist.begin(); _l !=  swaplist.end(); ++_l){
            if(_l->second) ans.swapcol(_l->first.first, _l->first.second);
            else ans.swaprow(_l->first.first, _l->first.second);
        }
#endif
        //cout << "swap back: ans = " << ans << "\n";
        return ans;
    }
/*
* 方阵 求值
 */
template<typename mT>
    mT det(matrix<mT>  m){
        if(m.nrow() != m.ncol()){
            fprintf(stderr, "the matrix(%luX%lu) cannot calculate det()!", m.nrow(), m.ncol());
            exit(1);
        }
        mT ans(1);
        //cout << " m = " << m << "\n";
        // start swap
        for(size_t _t = 0; _t < m.nrow(); ++_t){
            //cout << "_t = " << _t << '\n';
            // select main element
            size_t maxi = _t, maxj = _t;
            mT maxv = ::abs(m(_t,_t));
            for(size_t i = _t; i < m.nrow(); ++i){
                for(size_t j = _t; j < m.ncol(); ++j){
                    if( ::abs(m(i, j)) > maxv){
                        maxi = i; maxj = j; maxv = ::abs(m(i,j));
                    }
                }
            }
            if(maxi != _t) {
                m.swaprow(_t, maxi);
            }
            if(maxj != _t){
                m.swapcol(_t, maxj);
            }
            ans *= m(_t, _t);
            m[_t][_t] = 1 / m(_t, _t);
            for(size_t j = 0; j < m.nrow(); ++j)if(j != _t)
                m[_t][j] *= m(_t, _t);
            for(size_t i = 0; i < m.nrow(); ++i)if(i != _t)
                for(size_t j = 0; j < m.ncol(); ++j)if(j != _t)
                    m[i][j] -= m(i, _t) * m(_t, j);
            for(size_t i = 0; i < m.nrow(); ++i)if(i != _t)
                m[i][_t] = - m(i, _t) * m(_t, _t);
        }
        return ans;
    }
/*
 * 取出对角线元素作为数组
 */
template<typename mT>
    vector<mT> diag(const matrix<mT>& m)
    {
        int r = m.nrow(), c = m.ncol();
        int i = 0;
        vector<mT> ans;
        while(i < r && i < c){
            ans.push_back(m(i,i));
            ++i;
        }
        return ans;
    }
template<typename mT>
    matrix<mT> sqrt(const matrix<mT>& m)
    {
        matrix<mT> ans(m);
        for(size_t i = 0; i < m.nrow(); ++i)
            ans[i] = sqrt(m[i]);
        return ans;
    }
template<typename eT>
    matrix<eT> matrix_one(size_t k)
  {
      matrix<eT> ans(k,k,0);
      for(size_t ii=0; ii<k; ++ii)
      {
         ans[ii][ii] = eT(1);
      }

      return ans;
  }
template<typename mT>
	matrix<mT> rowmean(matrix<mT>& m)
	{
		vector<mT> vv;
		mT cc;
		for(size_t i=0;i<m.nrow();++i)
		{
			cc = mean(m.get_row(i));
			vv.push_back(cc);
			
		}
		matrix<mT> ans(vv);
		return ans;
	}

template<typename mT>
	matrix<mT> colmean(const matrix<mT>& m){
	matrix<mT> ans(1,m.ncol());
	for(size_t i = 0; i < ans.ncol(); ++i)
	{
		ans(0,i) = mean(m.get_col(i));
		
	}
	return ans;
}
/**Cholesky Decomposition**/
template<typename mT>
    matrix<mT> chol(const matrix<mT>& A)
    {
	matrix<mT> L(A.nrow(), A.ncol(), .0);
        assert(A.nrow() == A.ncol());
        size_t n = A.nrow();
        for(size_t i = 0; i < n; ++i)
            for(size_t j = 0; j < i + 1; ++j){
                mT s = 0;
                for(size_t k = 0; k < j; ++k)
                    s += L(i,k) * L(j,k);
                      L(i,j) = (i == j) ?
                        	sqrt(A(i,i) - s) :
                                 ((A(i,j) - s) / L(j,j));
             }
             return L;	

}

template<typename mT>
    ostream& operator<<(ostream& o, const matrix<mT>& v){
        o << "matrix("<< v.nrow() <<"X"<< v.ncol() <<"){\n";
#if __cplusplus >= 201103L
        for(auto _v:v){
            for(auto __v: _v) o << __v << ",";
            o<<'\n';
        }
#else
        for(typename vector< vector< mT> >::const_iterator _v = v.begin(); _v != v.end(); ++_v){
            for(typename vector<mT>::const_iterator __v = _v->begin(); __v != _v->end(); ++__v)
                o << *__v << ",";
            o<<'\n';
        }
#endif
        return o<<"}";
    }

#endif /* ifndef XG_MATH_MATRIX_H */

