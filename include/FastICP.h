// FastICP.h 
//
//Reference:
//Besl, P., & McKay, N. (1992). A method for registration of 3-D shapes. 
// IEEE Transactions on pattern analysis and machine intelligence, 239-256.


/*

Copyright (c) 2017, Xu CHEN, Fujian , China.(2017/06)
You can contact the author at <cx495086@outlook.com>

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the <organization> nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL ETH-ASL BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/


#ifndef FASTICP_H_
#define FASTICP_H_


#include "omp.h"
#include <tchar.h>
#include <ctime> 
#include <stdlib.h>
#include <string>
#include <fstream>
#include <iostream>
#include <algorithm> //
//#include <CString>
#include "nabo/nabo.h"
#include "nabo/Eigen/Eigenvalues"
#include <windows.h>
#include <comdef.h>
//using namespace std;
//using namespace Nabo;
using namespace Eigen;

	typedef struct{
		Matrix3d RT;      //3*3 rotation matrix
		Vector3d TR;      //3*1 translate matrix
		double* data_out;//target print out data
		double res;      //error
		int  iter;   //iter 
	} ICPOUT;

class FastICP
{
public:
	FastICP();


private:
	
	void rotationToQuaternion(double* RT,double* Q);
	void quaternionToRotation(double* Q,double* R);
	void extrapolateInTransformSpace(int t,double* dq,double* qs,double* err);
	double extrapolate(double* v,double* d,double vmax);
	void datafit(double *x,double *y,int n,double *a,int m,double *dt);
	
	
public:
	__int64 CountLines(char *filename)  ;
	int ReadTXT(char* fileName, double* matData);
	double mulknn(double* target,double* model,__int64 trow,__int64 mrow,double* dataclose,int k=1);
	void icp_alg(double* model,__int64 mrow,double* target,__int64 trow,ICPOUT* result,double res,int maxiter,int miniter);

	ICPOUT* DoICP(char* modelpath,char* datapath);

};



#endif //FASTICP_H_