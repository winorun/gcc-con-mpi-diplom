#include "reshotka.hpp"

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
	 
	POINT_ON_RESHOTKA porFist(&argc, &argv,0.5,0.5);
	
	porFist.setU(u);
	porFist.setG(g);
	porFist.setPhi_0(phi_0);
	porFist.setPhi_1(phi_1);
	porFist.setBoundary(boundary);
	porFist.voidMain();
//	porFist.printDebag();

    return 0;
}
