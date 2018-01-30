#include<iostream>
#include<fstream>
#include<cstdlib>
using namespace std;
void generateNodes(int );
float getRandom()
{
	return (float) rand()/ (float) RAND_MAX;
}
/*void printPair(float x,float y)
{
	cout<<x<<","<<y<<endl;
}*/

int main()
{
    //srand(10);
	ofstream output;
        output.open ("100.txt");
	int node;
	cerr<<"enter the no. of nodes"<<endl;
	cin>>node;
	float no_x,no_y;
	for(int i=0;i<node;i++)
	{
		no_x=getRandom();
		no_y=getRandom();
		output<<no_x<<" "<<no_y<<endl;

	}
	//generateNodes(node);

	return 0;
}
/*void generateNodes(int n)
{
	float no_x,no_y;
	for(int i=0;i<n;i++)
	{
		no_x=getRandom();
		no_y=getRandom();
		//printPair(no_x,no_y);
	}
}*/
