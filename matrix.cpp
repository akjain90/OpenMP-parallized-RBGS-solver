#include <assert.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include "matrix.h"
#include<omp.h>
#define PI 3.14159265358979323846

using namespace std;
        
	matrix::matrix(int ix, int iy)
	{
		#pragma omp parallel 
		{
		#pragma omp single
		{
		cout<<"Solving the PDE with "<<omp_get_num_threads()<<" threads \n";
		nx = ix;
		ny = iy;
		hx = 2.0/(nx);
		hy = 1.0/(ny);
		// Dynamic memory allocation to ensure golden touch policy
		v = new double [(nx+1)*(ny+1)];     
		res = new double [(nx*ny)];
		rhs = new double [(nx+1)*(ny+1)];
		temp_hx = 2*PI*hx;
	    temp_hy = 2*PI*hy;
		k_2 = 4*PI*PI;
		hx_2 = hx*hx;
		hy_2 = hy*hy;
		temp_pi = k_2*hx_2*hy_2;
		}
		}
	};
	
	//calculation of f(x,y) (RHS)//
	void matrix::setrhs()
	{
		int i,j;
		double val=0;		
		double temp_j = 0;
		#pragma omp parallel for private(j,i,val,temp_j) schedule (static)
		for (j=0;j<=ny;j++)
			{temp_j = temp_hy*j;
				for (i =0;i<=nx;i++)
				{
					val= temp_pi*sin(temp_hx*i)*sinh(temp_j);
					rhs[i+((nx+1)*j)] = val;
				}
			}
	};
	
	void matrix::setBoundary()
	{   double temp;
		double pi2 = 2*PI;
		#pragma omp parallel for private(temp) schedule (static)
		for( int i=0; i<=nx; i++)
		{ set(i,0,0);       // bottom boundary
			temp = sin(2*PI*hx*i)*sinh(pi2);
		  set(i, ny, temp);  // top boundary
		}

		#pragma omp parallel for schedule (static)
		for(int j=0; j<=ny; j++)
		{ set(0,j,0);       // left boundary
		  set(nx, j, 0);     // right boundary
		}
	};	

	double matrix::get(int i,int j)
	{
		return v[i+((nx+1)*j)];
	};

	void matrix::set(int i,int j,double val)
	{
		v[i+((nx+1)*j)] = val;
	};
	
	void matrix::RBGS(int k)
	{
		double temp=0; 
		double coeff= temp_pi+(2*hy_2+2*hx_2);
		double inv = 1/coeff;
		
		//RED//
		while(k>0)
		{
		#pragma omp parallel for private (temp) schedule (static)
			for (int j =1 ;j<ny;j+=2)
		{
			for (int i =1 ;i<nx;i+=2)
			{
				temp = inv*(rhs[i+((nx+1)*j)]+(hy_2*(get(i+1,j)+get(i-1,j)))+(hx_2*(get(i,j+1)+get(i,j-1))));
				set(i,j,temp);
			}
		}
		
		#pragma omp parallel for private (temp) schedule (static)
		for (int j =2 ;j<ny;j+=2)
		{
			for (int i =2 ;i<nx;i+=2)
			{
				temp = inv*(rhs[i+((nx+1)*j)]+(hy_2*(get(i+1,j)+get(i-1,j)))+(hx_2*(get(i,j+1)+get(i,j-1))));
				set(i,j,temp);
			}
		}
		
		//BLACK//
		#pragma omp parallel for private (temp) schedule (static)
		
		for (int j =1 ;j<ny;j+=2)
		{
			for (int i =2 ;i<nx;i+=2)
			{
				temp = inv*(rhs[i+((nx+1)*j)]+(hy_2*(get(i+1,j)+get(i-1,j)))+(hx_2*(get(i,j+1)+get(i,j-1))));
				set(i,j,temp);
			}
		}
		#pragma omp parallel for private (temp) schedule (static)
		for (int j =2 ;j<ny;j+=2)
		{
			for (int i =1 ;i<nx;i+=2)
			{
				temp = inv*(rhs[i+((nx+1)*j)]+(hy_2*(get(i+1,j)+get(i-1,j)))+(hx_2*(get(i,j+1)+get(i,j-1))));
				set(i,j,temp);
			}
		}
		--k;
			
		}
	};
	// Calculation of residual norm
	void matrix::residual()
	{
		double temp = 0;
		double temp_hx2 = (-1)/hx_2;
		double temp_hy2 = (-1)/hy_2;
		double temp_h = temp_hx2*temp_hy2;
		for(int j =1;j<ny;j++)
		{
			for (int i =1;i<nx;i++)
			{
					temp = ((temp_hx2*(get(i-1,j)-2*get(i,j)+get(i+1,j)))+(temp_hy2*(get(i,j-1)-2*get(i,j)+get(i,j+1)))+(k_2*get(i,j)));
					res[(i-1)+(nx*(j-1))] = temp -(temp_h*rhs[i+((nx+1)*j)]);
			}
			
		}
		
		double sum_r = 0;
		for (int j =0 ;j<nx*ny;j++)
		{

				sum_r = sum_r+res[j]*res[j];
				
		}
		sum_r = sum_r/(nx*ny);
		sum_r = pow(sum_r,0.5);
		cout<<"Residual norm :"<<sum_r<<endl;
	};
	
	void matrix::write()
	{
		cout<<"Writing the results to solution.txt :"<<endl;
		ofstream fout;
		fout.open("solution.txt");
		
		fout<<"# x y u(x,y)"<<endl;
		for(int j =0;j<=ny;j++)
		{
			for (int i =0;i<=nx;i++)
			{
				fout<<i*hx<<" "<<j*hy<<" "<<get(i,j)<<endl;
			}
			fout<<endl;
		}
	};

	void matrix::release()
	{
		delete [] v;
		delete[] res;
		delete[] rhs;
	};
