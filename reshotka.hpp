#ifndef RESHOTKA
#define RESHOTKA

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>  
#define DRAND double (rand())/ double(RAND_MAX)
#define DINAMIC_ALG 0x01000000
#define RESHOTKA_METOS 0x00100000
class POINT_ON_RESHOTKA{
	int rank,size;// mpi 
	int flag; // flag
	int porColDrave; // количество траекторий
	int N;//количество траекторий для одного процесса.
	int N_2; //количество точек для одного процесса
	double porX, porY;//=0.5e-00;
	double timeRun;

	static const double PI=3.141592e-00;
	
	double (*u) (double , double);
	double (*g) (double , double);
	double (*phi0) (double , double);
	double (*phi1) (double , double);
	int (*boundary) (double &, double &);
	double U,Disp;
	
	double voidSphere(int &N);
	double diam(double , double );
	void staticSphere();
	void dinamicSphere();
public:
	POINT_ON_RESHOTKA(int *, char***,double,double);
	~POINT_ON_RESHOTKA();
	
	void setU(double (*use) (double , double));
	void setG(double (*use) (double , double));
	void setPhi_0(double (*use) (double , double));
	void setPhi_1(double (*use) (double , double));
	void setBoundary(int (*use) (double &, double &));
	void printDebag();
	void printResult();
	int init(double,double);
	void mainRun();
	void setFlag(int fl){flag=flag | fl;};};
#endif
