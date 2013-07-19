#include <iostream>
#include <cmath> 
#include <unistd.h>
#include <stdlib.h>
#include "reshotka.hpp"
using namespace std;

const double epsilon=0.001e-00;

 double u(double x, double y) {
	 return double(exp(x)*exp(y));}
 
 double g(double x, double y){
	 return double(-4.0*exp(x)*exp(y));}
	 
 double phi_0(double x, double y){
	 return double(exp(x)*exp(y)); }
	 
 double phi_1(double x, double y){
	 return double(2.0*exp(x)*exp(y));}

 int boundary(double &x, double &y) { 
	 if (x<epsilon){x=0.0; return 0;};
	 if (x>(1-epsilon)){x=1.0; return 0;};
	 if (y<epsilon){y=0.0; return 0;};
	 if (y> (1-epsilon)){y=1.0; return 0;};
	 return 1; // return true if point on boundary or mod y,x on 1/0
 }
  
int main (int argc, char* argv[])
{
	int N;
	for(int rez=getopt(argc,argv,"n:");rez!=-1;rez=getopt(argc,argv,"n:"))
	switch (rez){
		case 'n':
			N = atoi(optarg);
			if(!N)N=1000;
			break;
	};
	
	if(libInit(&argc,&argv)&LIB_ERROR)return 1;
	setFunctions(u,g,phi_0,phi_1);
	setBoundary(boundary);
	Flag flag;
	Rezult rez = libRunComputing(Point(0.5,0.5),flag,N);
	if(!rez.iStatus){
	cout << rez.dDisp << "\n";
	cout << rez.dRezSum << "\n";}
	libClose();
    return 0;
}
