
//DEMO FastICP



#include "FastICP.h"

using namespace std;
using namespace Nabo;
using namespace Eigen;

#define ITERMAX 200
#define ITERMIN 5
#define pi 3.1415926535

int _tmain(int argc, _TCHAR* argv[])
{
	
	char* modelpath="model.txt";
	char* datapath="data.txt";

	FastICP DoingICP;
	FastICP::ICPOUT* result=DoingICP.DoICP(modelpath,datapath);

	if (result != NULL)
		cout<<endl<< "ICP Done Successful, the result shown as below:"<<endl<<endl;

	cout<<"result->RT: "<<endl<<result->RT<<endl<<endl;
	cout<<"result->TR: "<<endl<<result->TR<<endl<<endl;
	cout<<"result->res: "<<endl<<result->res<<endl<<endl;
	

	system("pause");
	return 0;
}
