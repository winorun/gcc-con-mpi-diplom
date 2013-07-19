#include <iostream>
#include <cmath> 
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
 
 //~ void printResult(){//времянка убрать принт нах
	//~ if(!rank){
		//~ if(flag&0x100){
			//~ printf("N - %i\n",N);
			//~ if((flag&0x1)&&(flag&0x10000000))printf("Presise solution:\t %f\n",u(porX,porY));
			//~ if(flag&0x1)printf("Numerical solution:\t %f\n",U);
			//~ if(flag&0x1)printf("Delta:\t\t\t %f\n",u(porX,porY)-U);
			//~ if(flag&0x1)printf("Disp:\t\t\t %f\n",Disp);
		//~ }else{
			//~ if(flag&0x200)printf("PS\t\t NS\t\t Delta\t\t Disp\n");
			//~ printf("%f\t %f\t %f\t %f\n",u(porX,porY),U,u(porX,porY)-U,Disp);
			//~ };
//~ }}
 
int main (int argc, char* argv[])
{
	
	if(libInit(&argc,&argv)&LIB_ERROR)return 1;
	setFunctions(u,g,phi_0,phi_1);
	setBoundary(boundary);
	libClose();
    return 0;
}
