#include "reshotka.hpp"
#include <unistd.h>
#include <mpi.h>
#include <stdlib.h>
#include <cmath>  

#define DRAND double (rand())/ double(RAND_MAX)

namespace{
	const double PI=3.141592e-00;
	
	Flag libGlobalFlag;
	int mpiRank,mpiSize; 
	FunWith2Par funU,funG,funPhi[2];
	int (*funBoundary) (double &, double &);
}

Point::Point(){
	dX=0.5;
	dY=0,5;
}

Point::Point(double x,double y){
	dX=x;
	dY=y;
}

//private funtion
 double diam(double x, double y){//return min R
	 if (x> fabs(1-x)) x= fabs(1-x);
	 if (y>fabs(1-y)) y=fabs(1-y); 
	 return (x>=y)?y:x;
}
	  
//public funtion

Flag libClose(){
	MPI_Finalize();
	libGlobalFlag=Flag(libGlobalFlag|LIB_MPI_CLOSE);
	return libGlobalFlag;}

Flag libInit(int *argc, char** argv[]){
	if(libGlobalFlag&LIB_INIT)return Flag(libGlobalFlag|LIB_ERROR);
	if(libGlobalFlag&LIB_ERROR)return libGlobalFlag;
	MPI_Init (argc, argv);
	MPI_Comm_rank (MPI_COMM_WORLD, &mpiRank);
	MPI_Comm_size (MPI_COMM_WORLD, &mpiSize);
	srand(mpiRank);
	libGlobalFlag=LIB_INIT;
	return libGlobalFlag;
}

Flag setFunctions(FunWith2Par U=0,FunWith2Par G=0,FunWith2Par Phi_0=0,FunWith2Par Phi_1=0){
if(libGlobalFlag&LIB_ERROR)return libGlobalFlag;
if(U){
	funU=U;
	libGlobalFlag=Flag(libGlobalFlag|U_EXISTS);}
if(G){
	funG=G;
	libGlobalFlag=Flag(libGlobalFlag|G_EXISTS);}
if(Phi_0&&Phi_1){
	funPhi[0]=Phi_0?Phi_0:funPhi[0];
	funPhi[1]=Phi_1?Phi_1:funPhi[1];
	libGlobalFlag=Flag(libGlobalFlag|PHI_FULL_EXISTS);}
return libGlobalFlag;
};

Flag setBoundary(int (*use) (double &, double &)){
	funBoundary=use;
	libGlobalFlag=Flag(libGlobalFlag|BOUNDARY_EXISTS);
	return libGlobalFlag;
}

//~ 
//~ double Point_ON_RESHOTKA::voidSphere(int &N){
	//~ double startTime = MPI_Wtime();
	//~ N_2 = 0;
//~ for(int i=N;i;i--){
	//~ double x=porX,y=porY;
	//~ double S1=0,S=0;
	//~ 
		//~ while (boundary(x,y)){
			//~ N_2++;
		   //~ double d=diam(x,y);
		   //~ 
		   //~ double alpha= DRAND;
		   //~ double omega1=cos(2*PI*alpha); 
		   //~ double omega2=sin(2*PI*alpha);
		   //~ 
		   //~ alpha= DRAND;
		   //~ double om1=cos(2*PI*alpha);
		   //~ double om2=sin(2*PI*alpha);
		   //~ 
		   //~ double alpha1=0;
		   //~ double alpha2=0;
		   //~ do{
			   //~ alpha1=DRAND;
			   //~ alpha2=DRAND;
			   //~ alpha2=4*alpha2/exp(1);
		   //~ }while(alpha2>(-4*alpha1*log(alpha1)));
		   //~ 
			//~ double nu=alpha1*d;
			//~ S+=(S1-((d*d-nu*nu-nu*nu*log(d/nu))/log(d/nu)))*d*d*g(x+nu*om1,y+nu*om2)/16;
			//~ S1-=d*d;
			//~ x=x+omega1*d;
			//~ y=y+omega2*d;
	   //~ };
	   //~ S+=(S1/4)*phi[1](x,y)+phi[0](x,y);
	   //~ U+=S/porColDrave;
	   //~ Disp+=(S*S)/porColDrave;
   //~ }
		//~ return MPI_Wtime()-startTime/N_2;
//~ }
//~ void Point_ON_RESHOTKA::mainRun(){
	//~ if(flag&DINAMIC_ALG)dinamicSphere();else staticSphere();
	//~ }
//~ void Point_ON_RESHOTKA::dinamicSphere(){
	//~ timeRun = MPI_Wtime();	
	//~ double minSpeedIndex=-1;
	//~ double speedIndex;
	//~ double kIndex=1;
	//~ double timeAndCol[2];
	//~ int modPril=size;
	//~ if(rank){
		//~ int N=porColDrave/(2*modPril);
		//~ porColDrave-=N*modPril;
		//~ while(N){
			//~ speedIndex=voidSphere(N);
			//~ MPI_Send(&speedIndex, 1,MPI_DOUBLE,0, 0,MPI_COMM_WORLD);
			//~ MPI_Recv(&N, 1,MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		//~ };
	//~ }
		//~ else {
			//~ do{
				//~ int N=0;
				//~ MPI_Status status;
				//~ MPI_Recv(&speedIndex, 1,MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD,&status);
				//~ if(minSpeedIndex==-1)minSpeedIndex=speedIndex;
				//~ if(minSpeedIndex>speedIndex)minSpeedIndex=speedIndex;
				//~ kIndex=minSpeedIndex/speedIndex;
				//~ if(porColDrave<50){modPril--;N=porColDrave;porColDrave=0;}
				//~ if(porColDrave>=50){
					//~ N=(porColDrave)/(2*modPril);// verni kIndex
					//~ if(!N)N=porColDrave;
					//~ porColDrave-=N;
					//~ }
				//~ MPI_Send(&N, 1,MPI_INT,status.MPI_SOURCE, 0,MPI_COMM_WORLD);
			//~ }while(modPril);
		//~ }
	//~ double inbuf[2],outbuf[2]={U,Disp};// <- говно. дублирование.
	//~ MPI_Reduce(outbuf, inbuf, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	//~ timeRun=MPI_Wtime()-timeRun;
	//~ if(!rank){
		//~ U=inbuf[0];
		//~ Disp=inbuf[1];
		//~ Disp=sqrt(fabs(Disp-U*U)/porColDrave);
	//~ }	
//~ }
//~ 

Rezult libRunComputing(Point startPoint,Flag & flag,int iN){
	Rezult rez={0,0,0};
	return rez;
	}

//~ void Point_ON_RESHOTKA::staticSphere(){
	//~ timeRun = MPI_Wtime();
	//~ if(!rank)N = porColDrave/size-size;else N= porColDrave/size+1;
	//~ voidSphere(N);
	//~ double inbuf[2],outbuf[2]={U,Disp}; // <- говно, дублирование.
	//~ MPI_Reduce(outbuf, inbuf, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	//~ timeRun=MPI_Wtime()-timeRun;
	//~ if(!rank){
		//~ U=inbuf[0];
		//~ Disp=inbuf[1];
		//~ Disp=sqrt(fabs(Disp-U*U)/porColDrave);
	//~ }	
//~ }
  
