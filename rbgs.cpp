#include <iostream>
#include <fstream>
#include <vector>
#include "matrix.cpp"
#include "Timer.h"
#include <cmath>
#define PI 3.14159265358979323846

using namespace std;
int main(int argc,char *argv[])
{
	if (argc==4)
	{
	int nx,ny,c;
	nx = atoi(argv[1]);
	ny = atoi(argv[2]);
	c  = atoi(argv[3]);
	cout<<"Number of grids in x: "<<nx<<endl;
	cout<<"Number of grids in y: "<<ny<<endl;
	cout<<"Total number of iterations "<<c<<endl;
	double time =0;
	siwir::Timer timer;
	timer.reset();
	matrix u(nx,ny);
	u.setrhs();
	u.setBoundary();
	u.RBGS(c);
	time = timer.elapsed();
	cout << "Calculation for PDE took " << time<< " seconds\n";
	u.residual();
	u.write();
	u.release();
	}
}
