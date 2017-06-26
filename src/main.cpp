
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
	ICPOUT* result=DoingICP.DoICP(modelpath,datapath);

	if (result != NULL)
		cout<<endl<< "ICP Done Successful, the result shown as below:"<<endl<<endl;

	cout<<"result->RT: "<<endl<<result->RT<<endl<<endl;
	cout<<"result->TR: "<<endl<<result->TR<<endl<<endl;
	cout<<"result->res: "<<endl<<result->res<<endl<<endl;
	
	ofstream outputTXT;
	outputTXT.close();
	outputTXT.clear();
	outputTXT.open( "RotateT.txt", fstream::in | fstream::trunc);
	outputTXT   <<result->RT.row(0)<<" "<<result->TR(0)<<endl
				<<result->RT.row(1)<<" "<<result->TR(1)<<endl
				<<result->RT.row(2)<<" "<<result->TR(2)<<endl
				<<0<<" "<<0<<" "<<0<<" "<<1<<endl;
	outputTXT.close();


	system("pause");free(result);
	return 0;
}
