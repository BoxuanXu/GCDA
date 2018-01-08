/*
 * $Id$
 */
#ifndef XG_MAP_MATH_H
#  define XG_MAP_MATH_H

#include <map>

using namespace std;

template<typename T1, typename T2>
map<T1, T2> operator+(map<T1, T2> a, const map<T1, T2> &b){
    for(map<T1, T2>::iterator it = b.begin(); it != b.end(); ++it)
        a[it->first] = a[it -> first] + it -> second;
    return a;
}

#endif /* ifndef XG_MAP_MATH_H */

