#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <vector>
using namespace std; 
double gaussrand()
{
	    static double V1, V2, S;
	    static int phase = 0;
	    double X;
	    srand(4385); 
		        if ( phase == 0 ) {
				        do {
						            double U1 = (double)rand() / RAND_MAX;
							                double U2 = (double)rand() / RAND_MAX;
									             
									            V1 = 2 * U1 - 1;
										                V2 = 2 * U2 - 1;
												            S = V1 * V1 + V2 * V2;
													            } while(S >= 1 || S == 0);
					         
					        X = V1 * sqrt(-2 * log(S) / S);
					    } else
						            X = V2 * sqrt(-2 * log(S) / S);
					             
					        phase = 1 - phase;
						 
						    return X;
}

int main()
{
	for(size_t i= 0;i < 20;++i)
	cout<<gaussrand()<<endl;
	return 0;
}
