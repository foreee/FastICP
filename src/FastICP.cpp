// FastICP.cpp : 

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



#include "FastICP.h"

using namespace std;
using namespace Nabo;
using namespace Eigen;

#define ITERMAX 200
#define ITERMIN 5
#define pi 3.1415926535

FastICP::FastICP()
{
	Eigen::initParallel();
	//Eigen::setNbThreads(6);

}
void FastICP::rotationToQuaternion(double* RT,double* Q)
{
	double* R=(double*) calloc(3*3,sizeof(double));
	Map<MatrixXd> RTBUF(RT,3,3);
	Map<MatrixXd> RBUF(R,3,3);
	RBUF=RTBUF.transpose();

	double t=R[0]+R[4]+R[8];
	double r,s;
	
	if(t>=0)
	{
		r=sqrt(1+t);s=0.5/r;
		Q[0]=0.5*r;Q[1]=(R[5]-R[7])*s;Q[2]=(R[6]-R[2])*s;Q[3]=(R[1]-R[3])*s;
	}
	else
	{
		double maxv=max(R[0],max(R[4],R[8]));
		if(maxv==R[0])
		{
			r =sqrt(1+R[0]-R[4]-R[8]);s=0.5/r;
			Q[0]=(R[5]-R[7])*s;Q[1]=0.5*r;Q[2]=(R[1]+R[3])*s;Q[3]=(R[6]+R[2])*s;
		}
		else if(maxv==R[4])
		{
		    r =sqrt(1+R[4]-R[0]-R[8]);s=0.5/r;
			Q[0]=(R[6]-R[2])*s;Q[1]=(R[1]+R[3])*s;Q[2]=0.5*r;Q[3]=(R[5]+R[7])*s;
		}
		else
		{
		    r =sqrt(1+R[8]-R[0]-R[4]);s=0.5/r;
			Q[0]=(R[1]-R[3])*s;Q[1]=(R[6]+R[2])*s;Q[2]=(R[5]+R[7])*s;Q[3]=0.5*r;
		}
	}
	free(R);
}
void FastICP::quaternionToRotation(double* Q,double* R)
{
	R[0]=Q[0]*Q[0]+Q[1]*Q[1]-Q[2]*Q[2]-Q[3]*Q[3];
	R[1]=2*(Q[1]*Q[2]-Q[0]*Q[3]);
	R[2]=2*(Q[1]*Q[3]+Q[0]*Q[2]);
	R[3]=2*(Q[1]*Q[2]+Q[0]*Q[3]);
	R[4]=Q[0]*Q[0]+Q[2]*Q[2]-Q[1]*Q[1]-Q[3]*Q[3];
	R[5]=2*(Q[2]*Q[3]-Q[0]*Q[1]);
	R[6]=2*(Q[1]*Q[3]-Q[0]*Q[2]);
	R[7]=2*(Q[3]*Q[2]+Q[0]*Q[1]);
	R[8]=Q[0]*Q[0]+Q[3]*Q[3]-Q[1]*Q[1]-Q[2]*Q[2];
}
void FastICP::extrapolateInTransformSpace(int t,double* dq,double* qs,double* err)
{
//Extrapolation in quaternion space. Details are found in:
//Besl, P., & McKay, N. (1992). A method for registration of 3-D shapes. 
// IEEE Transactions on pattern analysis and machine intelligence, 239-256.
	double dtheta1,dtheta2,n1,n2,n3;
	double angletheshold=10,scalefactor=25;
	Map<MatrixXd> QS(qs,7,Dynamic);
	Map<MatrixXd> DQ(dq,7,Dynamic);
	n1=DQ.col(t-3).norm();
	n2=DQ.col(t-2).norm();
	n3=DQ.col(t-1).norm();
	dtheta1=(180/pi)*acos(DQ.col(t-2).dot(DQ.col(t-3))/n1/n2);
	dtheta2=(180/pi)*acos(DQ.col(t-1).dot(DQ.col(t-2))/n2/n3);
	if(dtheta1<angletheshold && dtheta2<angletheshold )
	{
		double d[3]={err[t-1],err[t-2],err[t-3]};
		double v[3]={0,-n3,-n2-n3};
		double vmax=scalefactor*n3;
		double dv=extrapolate(v,d,vmax);
		if(dv!=0)
		{
			QS.col(t-1)=QS.col(t-1)+DQ.col(t-1)*dv/n3;
			QS.block<4,1>(0,t-1)=QS.block<4,1>(0,t-1)/QS.block<4,1>(0,t-1).norm();
		}
	}
}
double FastICP::extrapolate(double* v,double* d,double vmax)
//Extrapolation in quaternion space. Details are found in:
//Besl, P., & McKay, N. (1992). A method for registration of 3-D shapes. 
// IEEE Transactions on pattern analysis and machine intelligence, 239-256.
{
	double dv=0,v1,v2;
	double* err=(double*) calloc(3,sizeof(double));
	double sco1[3],sco2[3],list[4]={0,0,0,vmax};
	datafit(v,d,3,sco1,2,err); //linear fit
    datafit(v,d,3,sco2,3,err); //parabolic fit
	v1=-sco1[0]/sco1[1];
	v2=-sco2[1]/sco2[2]/2;
	list[1]=v1;list[2]=v2;
	sort(list,list+4);
	if(list[0]==0)
	{
		dv=list[1];return dv;
	}
	else if(list[1]==0)
	{
	    dv=list[2];return dv;
	}
	else if(list[2]==0)
		dv=list[3];
	free(err);
	return dv;
}
void FastICP::datafit(double *x,double *y,int n,double *a,int m,double *dt)
//  int n,m;   m=num(a);
//  double x[],y[],a[],dt[];
//y=a[0]+a[1]*(x-avg(x))+a[2]*(x-avg(x))^2+.....   
  { int i,j,k;
    double z,p,c,g,q,d1,d2,s[20],t[20],b[20];
    for (i=0; i<=m-1; i++) a[i]=0.0;
    if (m>n) m=n;
    if (m>20) m=20;
    z=0.0;
    for (i=0; i<=n-1; i++) z=z+x[i]/(1.0*n);
    b[0]=1.0; d1=1.0*n; p=0.0; c=0.0;
    for (i=0; i<=n-1; i++)
      { p=p+(x[i]-z); c=c+y[i];}
    c=c/d1; p=p/d1;
    a[0]=c*b[0];
    if (m>1)
      { t[1]=1.0; t[0]=-p;
        d2=0.0; c=0.0; g=0.0;
        for (i=0; i<=n-1; i++)
          { q=x[i]-z-p; d2=d2+q*q;
            c=c+y[i]*q;
            g=g+(x[i]-z)*q*q;
          }
        c=c/d2; p=g/d2; q=d2/d1;
        d1=d2;
        a[1]=c*t[1]; a[0]=c*t[0]+a[0];
      }
    for (j=2; j<=m-1; j++)
      { s[j]=t[j-1];
        s[j-1]=-p*t[j-1]+t[j-2];
        if (j>=3)
          for (k=j-2; k>=1; k--)
            s[k]=-p*t[k]+t[k-1]-q*b[k];
        s[0]=-p*t[0]-q*b[0];
        d2=0.0; c=0.0; g=0.0;
        for (i=0; i<=n-1; i++)
          { q=s[j];
            for (k=j-1; k>=0; k--)
              q=q*(x[i]-z)+s[k];
            d2=d2+q*q; c=c+y[i]*q;
            g=g+(x[i]-z)*q*q;
          }
        c=c/d2; p=g/d2; q=d2/d1;
        d1=d2;
        a[j]=c*s[j]; t[j]=s[j];
        for (k=j-1; k>=0; k--)
          { a[k]=c*s[k]+a[k];
            b[k]=t[k]; t[k]=s[k];
          }
      }
    dt[0]=0.0; dt[1]=0.0; dt[2]=0.0;
    for (i=0; i<=n-1; i++)
      { q=a[m-1];
        for (k=m-2; k>=0; k--)
          q=a[k]+q*(x[i]-z);
        p=q-y[i];
        if (fabs(p)>dt[2]) dt[2]=fabs(p);
        dt[0]=dt[0]+p*p;
        dt[1]=dt[1]+fabs(p);
      }
	//change scoeff
	double avg=(x[0]+x[1]+x[2])/3;
	if(m==3)  
	{
	    a[0]=a[0]-a[1]*avg+a[2]*avg*avg;
		a[1]=a[1]-2*a[2]*avg;
	}
  	if(m==2)
	{
		a[0]=a[0]-a[1]*avg;
	}
    return;
  }
__int64 FastICP::CountLines(char *filename)  
{  
    ifstream ReadFile;  
    __int64 n=0;  
    string tmp;  
    ReadFile.open(filename,ios::in);//ios::in means file can read only  
    if(ReadFile.fail())//file openned fail, return 0 
    {  
        return -1;  
    }  
    else//file exist  
    {  
        //#pragma omp parallel 
		while(getline(ReadFile,tmp,'\n'))  
        {  
            n++;  
        }  
        ReadFile.close();  
        return n;  
    }  
}
int FastICP::ReadTXT(char* fileName, double* matData)  
{  
	FILE *file;  
    __int64 LINES;    
    if(fopen_s(&file,fileName,"r"))
      {
           printf("Can not open the file!\n");
           return -1;
      } 
    else//file exist
    {  
        LINES=CountLines(fileName);  
                 
        int i;  
		for(i=0;i<LINES*3;i++) 
        {  
			fscanf_s(file,"%lf",&matData[i]);
			
        }  
		//cout<<matData[LINES*3-1]<<endl;
		fclose(file); //close file  
	return 0;
    }  
	
}
double FastICP::mulknn(double* target,double* model,__int64 trow,__int64 mrow,double* dataclose,int k)
{
	Map<MatrixXd> Mt(model,3,mrow);
	Map<MatrixXd> Tt(target,3,trow);
	MatrixXd M=Mt;  
	MatrixXd T=Tt;
	NNSearchD* nns = NNSearchD::createKDTreeLinearHeap(M);
	MatrixXi indices(k, T.cols());
	MatrixXd dists2(k, T.cols());
	nns->knn(T, indices, dists2, k, 0.2, NNSearchF::SORT_RESULTS);
	Map<MatrixXd> DATACLOSE(dataclose,3,trow);
	for(__int64 i=0;i<trow;++i)
	    DATACLOSE.col(i)=Mt.col(indices(0,i));

	delete nns;
	return dists2.row(0).cwiseSqrt().mean();//return the average dis of nearest points
}
void FastICP::icp_alg(double* model,__int64 mrow,double* target,__int64 trow,ICPOUT* result,double res,int maxiter,int miniter)
{
	//model means the static pointcloud
	//target is the moved pointcloud

	Matrix3d RT=Matrix3d::Identity();
	Vector3d TR=Vector3d::Zero();
	
	double* dataclose=(double*) calloc(trow*3,sizeof(double));
	Map<MatrixXd> DATACLOSE(dataclose,3,trow);
	//Map<MatrixXd> MODEL(model,3,mrow);
	Map<MatrixXd> TARGET(target,3,trow);
	Map<MatrixXd> DATAOUT(result->data_out,3,trow);
	//MatrixXd Dataclose=DATACLOSE.transpose();
	MatrixXd TARGETbuff=TARGET; //TARGET backup
	MatrixXd Target=TARGET.transpose();
	
	//for extrapolateInTransformSpace
	double* qs=(double*) calloc(maxiter*7,sizeof(double));
	double* dq=(double*) calloc(maxiter*7,sizeof(double));
	double* err=(double*) calloc(maxiter,sizeof(double));
	Map<MatrixXd> QS(qs,7,maxiter);
	Map<MatrixXd> DQ(dq,7,maxiter);

	for(int t=1;t<=maxiter;++t)
	{
		result->iter=t;
		clock_t start=clock();
		mulknn(target,model,trow,mrow,dataclose); 

		clock_t stop=clock();
		cout<<"once libnabo during time: "<<((double)(stop - start)) / CLOCKS_PER_SEC<<"s."<<"iter :"<<result->iter<<endl;
	    
		Vector3d t_mean= Target.colwise().mean().transpose();
		Vector3d d_mean= DATACLOSE.rowwise().mean();

		// caculate Covariance matrix
		//mat_anti=[(DATACLOSE-d_mean)*(Target-t_mean)]'
		Matrix3d mat_anti=((DATACLOSE-d_mean*VectorXd::Ones(trow).transpose())*(Target-VectorXd::Ones(trow)*t_mean.transpose())).transpose();
		//cout<<mat_anti<<endl<<endl;
		//svd
		JacobiSVD<MatrixXd> svd(mat_anti, ComputeFullU | ComputeFullV);

		//caculate RT=v*u' in each iteration
		RT=svd.matrixV()*svd.matrixU().transpose();
		if(RT.determinant()<0)
		{
			Matrix3d V=svd.matrixV();
			V.col(2)=-V.col(2);
			RT=V*svd.matrixU().transpose();
		}
	
		//caculate TR
		TR=d_mean-RT*t_mean;

		//caculate new target point cloud
		DATAOUT=Target.transpose();
 		Target=(RT*DATAOUT+TR*VectorXd::Ones(trow).transpose()).transpose();
		//update result
		result->RT=RT*result->RT;
		result->TR=TR+RT*result->TR;
		//caculate err
		err[t-1]=(Target-DATAOUT.transpose()).rowwise().norm().mean();
		result->res=err[t-1];
		
		//extrapolateInTransformSpace
		//double* resultRT=result->RT.data(); //save index by column
		//result->RT.transpose().data() is the same with result->RT.data(),so need to add a buff

		Matrix3d rtbuff;
		rtbuff=result->RT.transpose();
		rotationToQuaternion(rtbuff.data(),&qs[7*(t-1)]);
		QS.block<3,1>(4,t-1)=result->TR;
		if(t>1)
		{
			DQ.col(t-1)=QS.col(t-1)-QS.col(t-2);
			if(t>3)
			{
				extrapolateInTransformSpace(t,dq,qs,err);   //With extrapolation, we might be able to converge faster
			    // Update transformation and data

				quaternionToRotation(&qs[7*(t-1)],rtbuff.data());//updata result->RT
				result->RT=rtbuff.transpose();
				result->TR=QS.block<3,1>(4,t-1);//updata result->TR

				
				DATAOUT=result->RT*TARGETbuff+result->TR*VectorXd::Ones(trow).transpose();//updata result->dataout
				Target=DATAOUT.transpose();//update Target
			}
		}
		TARGET=Target.transpose();

		/*cout<<result->RT<<endl<<endl<<result->TR<<endl<<endl;
		  cout<<"qs: "<<qs[7*(t-1)]<<" "<<qs[7*(t-1)+1]<<" "<<qs[7*(t-1)+2]<<" "<<qs[7*(t-1)+3]<<" "<<qs[7*(t-1)+4]<<" "<<qs[7*(t-1)+5]<<" "<<qs[7*(t-1)+6]<<" "<<endl;
		*/
		stop=clock();
		cout<<"once ICP during time: "<<((double)(stop - start)) / CLOCKS_PER_SEC<<"s."<<"iter :"<<result->iter<<endl
		    <<"result->res: "<<endl<<result->res<<endl;
		if(t>2)
		{	
			if(abs(err[t-1]-err[t-2])<res && t>=miniter)
			{
				cout<<"Done"<<endl;
				break;
			}
		}

	}

	free(dataclose);free(err);free(qs);free(dq);
}
FastICP::ICPOUT* FastICP::DoICP(char* modelpath,char* datapath)
{
	clock_t start, stop;
	double res=0.0001;
	__int64 maxiter=ITERMAX, miniter=ITERMIN;

	//Read the points number of model and data
	cout << endl << "reading model points num ..." << endl;
	__int64 num_t = CountLines(modelpath);
	cout << " num_m: "<<num_t << endl;
	cout << endl << "reading data points num..." << endl;
	__int64 num_d = CountLines(datapath);
	cout << " num_d: "<<num_d << endl;

	//Read pointcloud of model and data
	double* T=(double*) calloc(3*num_t,sizeof(double));
	double* Tbuf=(double*) calloc(3*num_t,sizeof(double));
	double* D=(double*) calloc(3*num_d,sizeof(double));

	cout << endl << "reading model points ..." << endl;
	ReadTXT(modelpath, T);
	cout << endl << "reading data points ..." << endl;
	ReadTXT(datapath, D);

	ICPOUT* result=new ICPOUT;//C++ struct inition 
	result->RT=Matrix3d::Identity();
	result->TR=Vector3d::Zero();
	result->data_out=Tbuf;
	result->iter=0;
	result->res=0.0;

	start = clock();
	cout << endl << "Running ICP (point-to-point, no outliers)..." << endl;
	icp_alg(D,num_d,T,num_t,result, res,maxiter,miniter);// T has changed after ICP		
	stop = clock();

	// show the time cost
	printf("Run time = %g seconds\n",((double)(stop - start)) / CLOCKS_PER_SEC);

	cout << " model points number: "<<num_t << endl;
	cout << " data points number:  "<<num_d << endl;

	free(T);free(D);free(Tbuf);
	return result;
}

