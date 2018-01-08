/* *************************
 * Author: xg1990
 * Created Time:  
 * LastModified:  Mon 02 Dec 2013 03:27:42 PM CST
 * C File Name: 
 * ************************/
#include <iostream>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include "bpch.h"
using namespace std;



int main(int argc, char *argv[])
{
    bpch v1(argv[1]);
    //v1.readF(argv[1]);
    v1.showheader();
    cout<<v1.size()<<endl;
    
    for(size_t i = 0; i < v1.size(); ++i){
        v1[i].showheader();
       // v1[i].showdata();
    }

    return 0;
}


