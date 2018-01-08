/*
 * $Id$
 */
#ifndef XG_MATH_STAT_H
#  define XG_MATH_STAT_H

#include <boost/math/distributions/fisher_f.hpp>
#include <boost/math/distributions/students_t.hpp>
#include "xg_math_vector.h"
#include "xg_math_matrix.h"
#include "xg_str.h"
#include "xg_debug.h"
#include <vector> 
using boost::math::students_t;
using namespace xg;

template<typename T>
T cor(const vector<T>& x, const vector<T>& y, string method="pearson"){
    if(string("pearson").find(tolower(method)) != string::npos){
        //fprintf(stderr, "using pearson correlation coefficient!\n");
        return ( mean(x * y) - mean(x)*mean(y) ) / sqrt( mean(x * x) - mean(x)*mean(x) )/ sqrt( mean(y * y) - mean(y)*mean(y) );
    }
    else if (string("spearman").find(method) != string::npos){
        //fprintf(stderr, "using spearman correlation coefficient!\n");
        return cor( get_rank(x), get_rank(y) );
    }
    else{
        fprintf(stderr, "ERROR `method` parameter in cor()!\n");
        exit(1);
    }
}

template<typename T>
T cor(const vector<T>& x, const vector<T>& y, T & p_value, string method="pearson"){
    if(string("pearson").find(tolower(method)) != string::npos){
        size_t N = x.size();
        T r = cor(x,y);
	(abs(r) > 100000)||(isnan(abs(r)))?r=0:r=r;//xubx_test
        T t_value = r * sqrt( (N - 2) / (1 - r * r) );
	if(boost::math::isnan(t_value)){
            debug(x);
            debug(y);
            debug(cor(x,y));
            debug(1 - r * r);
            debug(t_value);
	  
        }
        p_value = 2 * cdf(students_t(N-2), -fabs(t_value));
        return r;
    }
    else if (string("spearman").find(method) != string::npos){
        //fprintf(stderr, "using spearman correlation coefficient!\n");
        return cor( get_rank(x), get_rank(y), p_value );
    }
    else{
        fprintf(stderr, "ERROR `method` parameter in cor()!\n");
        exit(1);
    }
}

template<typename T>
class pcorr{
    public:
        vector<vector<T> > data;
        matrix<T> r, pr, t_value, p_value;
        pcorr(){};
        pcorr(const vector<vector<T> > & d, string method="pearson"){
            operator()(d, method);
        }
        pcorr& operator()(const vector<vector<T> > & d, string method="pearson"){
            data = d;
            //for_i(i,0,d.size())
            //cout << d[i] << endl;
            r.resize(d.size(), d.size());
            t_value.resize(d.size(), d.size());
            p_value.resize(d.size(), d.size());
            for(size_t i=0; i < d.size(); ++i)
                for(size_t j=i; j < d.size(); ++j){
                    r[i][j] = r[j][i] = cor(d[i],d[j], method);
                }
            if(fabs(det(r) - 0 ) < 1e-9){
                pr = matrix<T>(d.size(), d.size());
                return (*this);
            }
            pr = ~r;
            //debug(~r);
            //debug(sqrt(sum(diag(pr) * diag(pr))));
            pr = (pr) / sqrt((matrix<T>(diag(pr)) * t(matrix<T>(diag(pr))))) * -1;
            for(size_t i=0; i < pr.nrow(); ++i)
                pr[i][i] = 0;
            for(size_t i=0; i < pr.nrow(); ++i)
                for(size_t j=i + 1; j <  pr.ncol(); ++j){
                    t_value[i][j] = t_value[j][i] = pr[i][j] * sqrt(data[0].size() - data.size()) / sqrt(1 - pr[i][j] * pr[i][j]);
                    if(boost::math::isnan(t_value[i][j])){
                        //debug(data[0]);
                        //debug(data[1]);
                        //debug(data[2]);
                        //debug(r);
                        //debug(pr);
                        //debug(t_value);
                        p_value[i][j] = p_value[j][i] = 1;
                    }
                    else{
                        p_value[i][j] = p_value[j][i] = 
                            //2 * ( 1 - gsl_cdf_tdist_P(fabs(t_value[i][j]) , data[0].size() - data.size()));
                            2 * cdf(students_t(data[0].size() - data.size()), -fabs(t_value[i][j]));
                    }
                }
            return (*this);
        }
};

template<typename T>
ostream& operator<<(ostream& o, const pcorr<T> & r){
    o << "partial correlation report:\n";
    for(size_t i=0; i < r.data.size(); ++i)
        o << "DATA_" << i << " : " << r.data[i] << endl;
    o << r.data.size() << " sample and " << r.data[0].size() << " observation \n";
    o << "cor matrix: \n" << r.r << endl;
    o << "pcor matrix: \n" << r.pr << endl;
    o << "t_value matrix: \n" << r.t_value<< endl;
    o << "p_value matrix: \n" << r.p_value<< endl;
    return o;
}

template<typename T>
T cov(const vector<T>& x, const vector<T>& y){
    return mean(x*y) - mean(x) * mean(y);
}

template<typename T>
T var(const vector<T>& x, long long df = -1){
    if(df == -1)
        return cov(x, x);
    else{
        vector<T> dx = x - mean(x);
        return sum(dx * dx) / df;
    }
}

template<typename T>
T sd(const vector<T>& x, size_t df = -1){
    return sqrt(var(x,df));
}

template<typename T>
T computeR2(const vector<T>& _y,const vector<T>& fited){
    vector<T> dy = _y - mean(_y), _resi = _y - fited ;
    T TSS = sum(dy * dy);
    T ESS = sum(_resi * _resi);
   
    return max(1 - ESS / TSS, 0.0);// if ESS > TSS, let R2 = 0
}

template<typename T>
T p_ftest(const T& fvalue, size_t df1, size_t df2){
    //fvalue = R2 / (K-1) * (N-K)/ (1 - R2),
    //debug(fvalue);
    //debug(df1);
    //debug(df2);
    if(std::isinf(fvalue))
        return 0;
    else if(fvalue > 0 &&  !std::isnan(fvalue))
        return 1 - cdf(boost::math::fisher_f_distribution<double>(df1,df2), fvalue);
    return 1;
}


/*
 * 回归模型基类
 */
template<typename T>
struct regress_base{
    size_t N/*观测数*/, K/*自变量数*/, df/*自由度*/;
    vector<vector<T> > x;
    vector<T> y, b, estimated, residuals;

    T R2, R2_adj;
    T  fvalue, pvalue,b2;
    T AIC;
    vector<T> tvalues, tvalues_p;

    regress_base(){}
    regress_base(const vector<vector<T> > & _x, const vector<T> & _y):
        N(_y.size()), K(_x.size()),df(N-K), x(_x), y(_y), b(K), tvalues_p(K)
    {}
    void estimate(){
        //debug("compute \'estimate\'");
        matrix<T> X = matrix<T>(x), Y=matrix<T>(y);
        rATA = ~(t(X)*X);
        b = (rATA * ( t(X) * Y )).col(0);
        estimated = (matrix<T>(x)*matrix<T>(b)).col(0);
        residuals = y - estimated;
        assessment();
        ttest();
        ftest();//回归系数使用t检验，而方程的整体性能使用F检验
        computeAIC();
        
    }
    void assessment(){
        //vector<T> dy = y - mean(y);
        //TSS = sum(dy * dy);
        ESS = sum(residuals * residuals);
        //RSS = TSS - ESS;
        //R2 = 1 - ESS / TSS;
        R2 = computeR2(y, estimated);//可决系数，拟合优度
        R2_adj = 1 - (1-R2)*(N-1)/df;
       
    }
    void ttest(){
        tvalues = b / sqrt( ESS / (N - K)) / sqrt(diag(rATA));
        //debug(ESS);
        //debug(b);
        //debug(N);
        //debug(K);
        //debug(rATA);
        //debug(tvalues);
        //debug(tvalues[0]);
        //debug(abs(tvalues[0]));
        //debug(fabs(tvalues[0]));
        //debug(N-K);
        //debug(tvalues[0]);
        //debug(cdf(students_t(N-K), -abs(tvalues[0])));
        for(size_t i = 0; i < K; ++i)
            tvalues_p[i] = (boost::math::isnan(tvalues[i]) || boost::math::isinf(tvalues[i])) ? 1 :2 * cdf(students_t(N-K), -abs(tvalues[i]));
        //debug(tvalues_p[0]);
    }
    void ftest(){
        fvalue = R2 / (K-1) * (N-K)/ (1 - R2);
        pvalue = p_ftest(fvalue, K-1, N-K);
        //if(std::isinf(fvalue))
            //pvalue = 0;
        //else if(fvalue > 0 &&  !std::isnan(fvalue))
            //pvalue = 1 - cdf(boost::math::fisher_f_distribution<double>(K-1,N-K), fvalue);
        //else
            //pvalue = 1;
        //assert(pvalue == p_ftest(fvalue, N, K));
        
    }
    void computeAIC(){
        T RSS = sum(residuals * residuals);
        //debug(residuals);
        //debug(RSS);
        int k = b.size() + 1, n = y.size();
        //debug(k);
        //debug(n);
        if(n <= k + 1)
            AIC = 1e100;
        else
            AIC= n * log( RSS / n) + 2 * k + 2 * k * (k + 1) / (n - k - 1);//防止过度拟合的等问题，AIC越小越好，其实这里应该是AICc，在样本比较小的情况下用。
        //cout<<"AIC have been calculated"<<endl;
        //debug(AIC);
    }
    template<typename tt>
        friend ostream& operator<<( ostream & o, const regress_base<tt>& reg);
    protected:
    T TSS, ESS, RSS;
    matrix<T> rATA;// A = { x00,.....xnn}, (A^TA)^{-1}
};

template<typename tt>
ostream& operator<<( ostream & o, const regress_base<tt>& reg){
    o << "regression model result:\n";
    o << "Coefficients:\n\tEstimate\tt value\tPr(>|t|)\n";
    for(size_t i = 0; i < reg.b.size(); ++i)
        o << "x" << i << "\t" << reg.b[i] << '\t' << reg.tvalues[i] << '\t' << reg.tvalues_p[i]<<endl;
    o << "estimated: " << reg.estimated << endl;
    o << "residuals: " << reg.residuals << endl;
    o << "R2= " << reg.R2 << ", adjusted R2  = " << reg.R2_adj << endl;
    o << "f statistic = " << reg.fvalue << " on df ("<< reg.K - 1 <<", " << reg.df <<") , pvalue= " <<reg.pvalue << endl;
    o << endl;
    return o;
}

/* 
 * 多元线性回归模型
 */
template<typename T>
struct regress: public regress_base<T>{
    using regress_base<T>::estimate;
    regress(){}
    regress(const vector<vector<T> >& _x, const vector<T>& _y):
        regress_base<T>(_x, _y){
            estimate();
            //debug("[regress] constructed");
        }
    template<typename tt>
        friend ostream& operator<<( ostream & o, const regress<tt>& reg);

};

template<typename tt>
ostream& operator<<( ostream & o, const regress<tt>& reg){
    o << (regress_base<tt>)reg;
    //o << "Multivariate regression model result:\n";
    //o << "Coefficients:\n\tEstimate\tt value\tPr(>|t|)\n";
    //for(size_t i = 0; i < reg.b.size(); ++i)
        //o << "x" << i << "\t" << reg.b[i] << '\t' << reg.tvalues[i] << '\t' << reg.tvalues_p[i]<<endl;
    //o << "estimated: " << reg.estimated << endl;
    //o << "residuals: " << reg.residuals << endl;
    //o << "R2= " << reg.R2<< endl;
    //o << endl;
    //cout<<"where is <<"<<endl;
    return o;
}

/*
 * 二元线性回归模型
 */
template<typename T>
struct regress2 : public regress<T> {
    //T b0, b1, R2;
    using regress_base<T>::x;
    using regress_base<T>::y;
    using regress_base<T>::N;
    using regress_base<T>::K;
    using regress_base<T>::estimated;
    using regress_base<T>::residuals;
    using regress_base<T>::estimate;
       
      
    T pvalue,b1;
      
        
    regress2(){}
    regress2(const vector<vector<T> >& _x, const vector<T>& _y){
                 regress<T> r2(_x,_y);
                 pvalue=r2.pvalue;
                 b1=r2.b[1];
                 //cout<<"end regress2!"<<endl;
              //regress<T>(vector<vector<T> >(vector<T>(_y.size(),1),_x), _y){}
         }
       /*regress<T>(vector<vector<T> >({vector<T>(_y.size(),1),_x}), _y){
            //debug(_x);
            //debug(_y);
            //assert(_x .size() == _y. size());
            //b1 = cov(_x,y)/ cov(_x,_x);
            //b0 = mean(y) - b1 * mean(_x);
            //estimated = _x * b1 + b0;
            //residuals = y - estimated;
        }*/
    //void ttest(){
        //T MSE = var(residuals, residuals.size() - 2);
    //}
    template<typename tt>
        friend ostream& operator<<( ostream & o, const regress2<tt>& reg);
    
};

template<typename T>
ostream& operator<<( ostream & o, const regress2<T>& reg){
    o << (regress<T>)reg;
    //o << "y = " << reg.b1 << "x + " << reg.b0 << "\n";
    return o;
}




/* 
 * 二元分段线性回归模型
 */

template<typename T>
struct regress2_seg2: public regress_base<T>{
    using regress_base<T>::x;
    using regress_base<T>::y;
    using regress_base<T>::N;
    using regress_base<T>::K;
    using regress_base<T>::estimate;
    using regress_base<T>::estimated;
    using regress_base<T>::tvalues_p;
    using regress_base<T>::R2;
    using regress_base<T>::AIC;

    size_t breakpoint;
    T seg1_R2, seg1_fvalue, seg1_pvalue;
    T seg2_R2, seg2_fvalue, seg2_pvalue;
    size_t N1, K1, N2, K2;
    double dAIC;

          


    regress2_seg2(){}
    regress2_seg2(const vector<T>& _x, const vector<T>& _y):
        regress_base<T>(vector<vector<T> >({ vector<T>(_y.size(),1), _x, _x}), _y), breakpoint(-1){
            //cout << "N = " << N << endl;
            T best_R2 = 0;
            for(size_t k = 3; k < _x.size() - 3; ++k){
                x[2] = (x[1] - x[1][k] ) * (x[1] > x[1][k]);
                //regress<T> mreg({vector<T>(N,1), _x, (_x - _x[k]) * (_x > _x[k])}, _y);
                estimate();
                if(R2 > best_R2){
                    best_R2 = R2;
                    breakpoint = k;
                }
                //cout << "k = " << k << ": \n" << mreg << endl;
            }
            if(breakpoint >= 0 && breakpoint < N){
                //debug(sub(y,0,breakpoint));
                //debug(subright(y,breakpoint));
                x[2] = (x[1] - x[1][breakpoint] ) * (x[1] > x[1][breakpoint]);
                //debug(tvalues_p[0]);
                estimate();
                //debug(tvalues_p[0]);
                seg1_R2 = computeR2(sub(y,0,breakpoint), sub(estimated, 0, breakpoint));
                seg2_R2 = computeR2(subright(y,breakpoint), subright(estimated, breakpoint));
                N1 = breakpoint + 1, N2 = N - N1 + 1;
                K2 = K, K1 = K2 - 1 ;
                //debug(N1);
                //debug(N2);
                //debug(K1);
                //debug(K2);
                seg1_fvalue = seg1_R2 / (K1 - 1)* (N1 - K1 )/ (1 - seg1_R2);
                seg2_fvalue = seg2_R2 / (K2 - 1)* (N2 - K2 )/ (1 - seg2_R2);
                seg1_pvalue = p_ftest(seg1_fvalue, N1, K1);
                seg2_pvalue = p_ftest(seg2_fvalue, N2, K2);
                dAIC = AIC - regress2<T>(_x, _y).AIC;
            }
        }
    T breakx()const{
        return x[1][breakpoint];
    }
};

template<typename tt>
ostream& operator<<( ostream & o, const regress2_seg2<tt>& reg){
    o << (regress_base<tt>)(reg);
    o << "Segmented regression model result:\n";
    o << " y = " << reg.b[0] << " + " << reg.b[1] << "x + " << reg.b[2] << "(x - " << reg.x[0][reg.breakpoint] <<") I(x>"<<  reg.x[0][reg.breakpoint] <<")\n";
    o << "breakpoint:x["<< reg.breakpoint <<"]= " << reg.breakx()<< endl;
    //o << "R2= " << reg.R2 << ", adjusted R2  = " << reg.R2_adj << endl;
    o << "For segment 1: R2 = " << reg.seg1_R2 << ", f_statistic= " << reg.seg1_fvalue << ", p-value=" << reg.seg1_pvalue << endl;
    o << "For segment 2: R2 = " << reg.seg2_R2 << ", f_statistic= " << reg.seg2_fvalue << ", p-value=" << reg.seg2_pvalue << endl;
    
    o << endl;
    return o;
}


#endif /* ifndef XG_MATH_STAT_H */

