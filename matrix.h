#ifndef MATRIX_H
#define MATRIX_H
#include <iostream>
#include <vector>
#include <assert.h>
#include <cmath>
#include <fstream>
#define PI 3.14159265358979323846

using namespace std;
class matrix
{
	private:
	double *v=NULL; 
	double *rhs=NULL;
	int nx,ny;
	double hx,hy;
	double temp_hx;
	double temp_hy;
	double hx_2,hy_2;
	double temp_pi;
	double k_2;
	double *res=NULL;
	
	
	public:
	
	matrix(int ix, int iy);
	
	double get(int i,int j);
	
	void set(int i,int j,double val);
	
	void setrhs();
	
	void setBoundary();
	
	void RBGS(int);
	
	void write();
	
	void residual();

	void release();
	
};

#endif
