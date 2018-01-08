/*
 * $Id$
 */
#ifndef XG_REGRESSION_H
#  define XG_REGRESSION_H


#include <iostream>
#include "xg_vector_math.h"
#include <boost/math/distributions/fisher_f.hpp>
#include <boost/math/distributions/students_t.hpp>
#include <cmath>

#ifdef DEBUG
//#include "xg_debug.h"
#endif

using namespace std;
using namespace boost::math;

template<typename T>
class xg_regress2{
    public:
        xg_regress2(){}
        xg_regress2(const vector<T>& x, const vector<T>& y):
            x(x),y(y){
                assert(x.size() == y.size());
                int n = x.size();
                T m_x = sum(x)/x.size(), m_y = sum(y)/y.size();
                vector<T> dx = x - m_x, dy = y - m_y;
                T m_dx2 = mean(dx * dx), m_dxdy = mean(dx * dy);

                a = m_dxdy/ m_dx2;
                b = m_y - a * m_x;

                estimated = x * a + b;
                Residuals = y - estimated;
                T d2 = sum(Residuals * Residuals);
                T s2 = d2 / (n - 2.0);

                T cov00 = s2 * (1.0 / n) * (1.0 + m_x * m_x / m_dx2);
                T cov11 = s2 * 1.0 / (n * m_dx2);
                //T cov01 = s2 * (-m_x) / (n * m_dx2);

                /*** start assessment   */
                b_std_err = sqrt(cov00);
                b_t_value = b / b_std_err;


                a_std_err=sqrt(cov11);
                a_t_value=a / a_std_err;
                //printf("n=%d, a_t_value = %f\n", n, a_t_value);
                if(!std::isnan(a_t_value) && n > 2 &&  isfinite(a_t_value))
                    a_p_value = 2 * cdf(students_t(n-2), -fabs(a_t_value));
                else
                    a_p_value = 1;

                T dl = n - 2;//degrees of liberty
                T sct = sum(dy*dy);
                R2 = 1 - d2 / sct;
                adj_R2 = 1 - T(n-1)/dl*(1-R2);
                F = R2 * dl / (1 - R2);
#ifdef DEBUG
                //out(x);
                //out(a);
                //out(b);
                //out(estimated);
                //out(y);
                //out(Residuals);
                //out(sum(dy*dy));
                //out(d2);
                //out(R2);
                //out(dl);
                //out(F);
#endif
                /* fisher_f_distribution exists in c++11 */
                if(std::isinf(F))
                    p_value = 0;
                else if(F > 0 &&  !std::isnan(F))
                    p_value= 1 - cdf(boost::math::fisher_f_distribution<double>(1,dl), F);
                else
                    p_value = 1;
                computeAIC();
            }
        vector<T> x,y;
        vector<T> estimated, Residuals;
        T AIC;
        T a, a_std_err, a_t_value, a_p_value;
        T b, b_std_err, b_t_value, b_p_value;
        T F, p_value;
        T R2, adj_R2;
    private:
        void computeAIC(){
            T RSS = sum(Residuals * Residuals);
            int k = 3, n = x.size();
            if(n <= k + 1)
                AIC = 1e100;
            else
                AIC= n * log( RSS / n) + 2 * k + 2 * k * (k + 1) / (n - k - 1);
        }
};


template<typename T>
ostream& operator<<(ostream& o, const xg_regress2<T>& r){
    o << "regress result\n";
    o << "x:" << r.x << endl;
    o << "y:" << r.y << endl;
    o << "y = " << r.a << "x + " << r.b << endl;
    o << "estimated:" << r.estimated << endl;
    o << "Residuals:" << r.Residuals << endl;
    o << "Coefficients\tEstimate\tStd. Error\tt value\tPr(>|t|)"<<endl;
    o << "Intercept\t"<< r.b <<"\t"<< r.b_std_err << "\t" << r.b_t_value <<"\t"<<r.b_p_value << endl;
    o << "x\t\t"<< r.a << "\t" << r.a_std_err <<"\t"<< r.a_t_value << "\t"<<r.a_p_value <<endl;
    o << "Multiple R-squared: "<< r.R2 <<",    Adjusted R-squared: "<< r.adj_R2 <<endl;
    o << "F-statistic:  "<< r.F <<" on 1 and "<< r.x.size() - 2 <<" DF,  p-value: "<<r.p_value<<endl;
    o << "RSS= " << sum(r.Residuals * r.Residuals) << ", AIC : " << r.AIC << endl;
    return o;
}

template<typename T>
class xg_seg2_regress2{
    public:
        xg_seg2_regress2(const vector<T>& x, const vector<T>& y):
            x(x), y(y){
                int n =x.size(), bp = 2;
                T m_x = mean(x), m_y = mean(y), r2;
                vector<T> dx = x - m_x, dy = y - m_y;
                //*******************//
                vector<T> x1, x2, y1, y2;

                R2 = -1e9;
                breakpoint = -1;
                for(bp = 3; bp < n - 2; ++bp)
                {
#ifdef DEBUG
                    //printf("tring bp=%d\n", bp);
#endif
                    x1 = sub(x,0,bp);
                    x2 = subright(x,bp);
                    y1 = sub(y,0,bp);
                    y2 = subright(y,bp);
                    xg_regress2<T> reg1(x1, y1);
                    xg_regress2<T> reg2(x2, y2);

                    vector<T> _estimated = reg1.estimated;
                    push_back(_estimated, reg2.estimated);
                    vector<T> _Residuals = _estimated - y;
                    T d2 = sum(_Residuals * _Residuals);
                    T sct = sum(dy*dy);
                    r2 = 1 - d2 / sct;
#ifdef DEBUG
                    //printf("_Residuals of SEG = ");cout << _Residuals << endl;
                    //printf("d2 of SEG = %lf\n", d2);
                    //printf("r2 of SEG = %lf\n", r2);
#endif
                    //printf("bp: %d, R2 : %f\n", bp, r2);

                    //cout << "reg1  R2:" << reg1.adj_R2 << " p : " << reg1.p_value<< endl;
                    //cout << "reg2  R2:" << reg2.adj_R2 << " p : " << reg2.p_value<< endl;
                    //if(  (reg1.p_value < 0.05 ||  reg2.p_value < 0.05) && r2 > R2)
                    if( r2 > R2)
                    {
                        Reg1 = reg1;
                        Reg2 = reg2;
                        R2 = r2;
                        breakpoint = bp;
                        estimated = reg1.estimated;
                        push_back(estimated, reg2.estimated);
                        Residuals = estimated - y;
                    }
                }
                computeAIC();
            };
        vector<T> x,y;
        int breakpoint;
        vector<T> estimated, Residuals;
        xg_regress2<T> Reg1, Reg2;
        T F;
        T R2, adj_R2;
        T AIC, dAIC;
    private:
        void computeAIC(){
            T RSS = sum(Residuals * Residuals);
            int k = 7, n = x.size();
            //printf(" n=%lf\n",  n);
            //printf(" RSS / n=%lf\n",  RSS / n);
            //printf("log( RSS / n)=%lf\n", log( RSS / n));
            AIC= n * log( RSS / n) + 2 * k + 2 * k * (k + 1) / (n - k - 1);
            dAIC = AIC - xg_regress2<T>(x, y).AIC;
        }
};

template<typename T>
ostream& operator<<(ostream& o, const xg_seg2_regress2<T>& r){
    if(r.breakpoint == -1)
    {
        o << "No significant breakpoint found! \n";
    }
    else{
        o << "break point is at x = "  << r.x[r.breakpoint] << endl;
        o << "x:" << r.x << endl;
        o << "y:" << r.y << endl;
        o << "FIRST SEGMENT:\n" << r.Reg1 << endl;
        o << "\n SECOND SEGMENT:\n" << r.Reg2 << endl;
        o << "estimated:" << r.estimated << endl;
        o << "Residuals:" << r.Residuals << endl;
        o << "R2 = " << r.R2 << endl;
        o << "RSS="<< sum(r.Residuals * r.Residuals) << ", AIC=" << r.AIC << ", dAIC= " << r.dAIC << endl;
    }
    return o;
}

#endif /* ifndef XG_REGRESSION_H */

