/*
 * $Id$
 */
#ifndef XG_MATH_BASE_H
#  define XG_MATH_BASE_H

#include <vector>

using namespace std;

template<typename T>
vector<T> rep(T v, size_t N){
    return vector<T>(N,v);
}

template<typename T>
vector<T> seq(const T& from, const T& to, const T& by = 1){
    vector<T> ans;
    T tmp = from;
    while(tmp <= to){
        ans.push_back(tmp);
        tmp+= by;
    }
    return ans;
}

template<typename T>
vector<T> seq(const T& from, const T& to, size_t length){
    return seq(from, to, (to - from) / (length - 1));
}

void get_next(const char* p, int plen ,int* next){
	int q=0,k=0;
	next[0] =0;
	for(q=1;q< plen;q++){
		while(k>0 && p[q] != p[k]){
			k=next[k-1];
		}
		if(p[q]==p[k])
			k++;
		next[q]=k;
	}
}

int kmp_match(const char* t,const char* p){
	int plen =strlen(p);
	int *next = (int*)malloc(sizeof(int)*plen);
	int q=0,i=0,tlen =strlen(t);
        if(NULL == t || tlen <= 0 || NULL == p || plen <= 0 )
		return -1;
	get_next(p,plen,next);
	for(i=0;i<tlen;i++){
		while(q>0 && p[q]!=t[i]){
			q=next[q-1];
		}
		if(p[q]==t[i])
			++q;
		if(q==plen){
			return i-plen+1;
		}
	}
	return -1;
 	free(next);
}
#endif /* ifndef XG_MATH_H */

