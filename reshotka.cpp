#include "reshotka.hpp"
#include <unistd.h>

	void POINT_ON_RESHOTKA::setU(double (*use) (double , double))
	{u=use;flag=flag | 0x10;};
	void POINT_ON_RESHOTKA::setG(double (*use) (double , double))
	{g=use;flag=flag | 0x20;};
	void POINT_ON_RESHOTKA::setPhi_0(double (*use) (double , double))
	{phi0=use;flag=flag | 0x40;};
	void POINT_ON_RESHOTKA::setPhi_1(double (*use) (double , double))
	{phi1=use;flag=flag | 0x80;};
	void POINT_ON_RESHOTKA::setBoundary(int (*use) (double &, double &))
	{boundary=use;flag=flag | 0x08;};


POINT_ON_RESHOTKA::POINT_ON_RESHOTKA(int *argc, char** argv[],double x,double y){
MPI_Init (argc, argv);
MPI_Comm_rank (MPI_COMM_WORLD, &rank);
MPI_Comm_size (MPI_COMM_WORLD, &size);
srand(rank);
porX=x;
porY=y;
flag=0x1;
for(int rez=getopt(*argc,*argv,"n:flc");rez!=-1;rez=getopt(*argc,*argv,"n:flc"))
	switch (rez){
		case 'n': porColDrave = atoi(optarg);if(!porColDrave)porColDrave=1000;break;
		case 'f': flag=flag | 0x100;break;
		case 'l':flag=flag | 0x200;break;
		case 'c':flag=flag | 0x10000000;break;
		
	};

}

POINT_ON_RESHOTKA::~POINT_ON_RESHOTKA(){
MPI_Finalize();
}
int POINT_ON_RESHOTKA::init(double x,double y)
{
	porX=x;
	porY=y;
	flag=1;
	return 1;
	}
void POINT_ON_RESHOTKA::printDebag(){
		printf("Выполнение %i на  %i заняло %f секунд при %i среднем длене %i общий путь\n",rank, N,timeRun,N_2/N,N_2);
	}

	
void POINT_ON_RESHOTKA::printResult(){//времянка убрать принт нах
	if(!rank){
		if(flag&0x100){
			printf("N - %i\n",N);
			if((flag&0x1)&&(flag&0x10000000))printf("Presise solution:\t %f\n",u(porX,porY));
			if(flag&0x1)printf("Numerical solution:\t %f\n",U);
			if(flag&0x1)printf("Delta:\t\t\t %f\n",u(porX,porY)-U);
			if(flag&0x1)printf("Disp:\t\t\t %f\n",Disp);
		}else{
			if(flag&0x200)printf("PS\t\t NS\t\t Delta\t\t Disp\n");
			printf("%f\t %f\t %f\t %f\n",u(porX,porY),U,u(porX,porY)-U,Disp);
			};
}}

double POINT_ON_RESHOTKA::voidSphere(int &N){
	double startTime = MPI_Wtime();
	N_2 = 0;
for(int i=N;i;i--){
	double x=porX,y=porY;
	double S1=0,S=0;
	
		while (boundary(x,y)){
			N_2++;
		   double d=diam(x,y);
		   
		   double alpha= DRAND;
		   double omega1=cos(2*PI*alpha); 
		   double omega2=sin(2*PI*alpha);
		   
		   alpha= DRAND;
		   double om1=cos(2*PI*alpha);
		   double om2=sin(2*PI*alpha);
		   
		   double alpha1=0;
		   double alpha2=0;
		   do{
			   alpha1=DRAND;
			   alpha2=DRAND;
			   alpha2=4*alpha2/exp(1);
		   }while(alpha2>(-4*alpha1*log(alpha1)));
		   
			double nu=alpha1*d;
			S+=(S1-((d*d-nu*nu-nu*nu*log(d/nu))/log(d/nu)))*d*d*g(x+nu*om1,y+nu*om2)/16;
			S1-=d*d;
			x=x+omega1*d;
			y=y+omega2*d;
	   };
	   S+=(S1/4)*phi1(x,y)+phi0(x,y);
	   U+=S/porColDrave;
	   Disp+=(S*S)/porColDrave;
   }
		return MPI_Wtime()-startTime/N_2;
}
void POINT_ON_RESHOTKA::mainRun(){
	if(flag&DINAMIC_ALG)dinamicSphere();else staticSphere();
	}
void POINT_ON_RESHOTKA::dinamicSphere(){
	timeRun = MPI_Wtime();	
	double minSpeedIndex=-1;
	double speedIndex;
	double kIndex=1;
	double timeAndCol[2];
	int modPril=size;
	if(rank){
		int N=porColDrave/(2*modPril);
		porColDrave-=N*modPril;
		while(N){
			//~ printf("%i\t start %i \n",rank,N);
			speedIndex=voidSphere(N);
			MPI_Send(&speedIndex, 1,MPI_DOUBLE,0, 0,MPI_COMM_WORLD);
			MPI_Recv(&N, 1,MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			//~ printf("%i\t sent %i \n",rank,N);
		};
	}
		else {
			do{
				int N=0;
				MPI_Status status;
				MPI_Recv(&speedIndex, 1,MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD,&status);
				if(minSpeedIndex==-1)minSpeedIndex=speedIndex;
				if(minSpeedIndex>speedIndex)minSpeedIndex=speedIndex;
				kIndex=minSpeedIndex/speedIndex;
				if(porColDrave<50){modPril--;N=porColDrave;porColDrave=0;}
				if(porColDrave>=50){
					N=(porColDrave)/(2*modPril);// verni kIndex
					if(!N)N=porColDrave;
					porColDrave-=N;
					}
				MPI_Send(&N, 1,MPI_INT,status.MPI_SOURCE, 0,MPI_COMM_WORLD);
			}while(modPril);
		}
	double inbuf[2],outbuf[2]={U,Disp};// <- говно. дублирование.
	MPI_Reduce(outbuf, inbuf, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	timeRun=MPI_Wtime()-timeRun;
	if(!rank){
		U=inbuf[0];
		Disp=inbuf[1];
		Disp=sqrt(fabs(Disp-U*U)/porColDrave);
		printResult();
	}	
}

void POINT_ON_RESHOTKA::staticSphere(){
	timeRun = MPI_Wtime();
	if(!rank)N = porColDrave/size-size;else N= porColDrave/size+1;
	voidSphere(N);
	double inbuf[2],outbuf[2]={U,Disp}; // <- говно, дублирование.
	MPI_Reduce(outbuf, inbuf, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	timeRun=MPI_Wtime()-timeRun;
	if(!rank){
		U=inbuf[0];
		Disp=inbuf[1];
		Disp=sqrt(fabs(Disp-U*U)/porColDrave);
		printResult();
	}	
}

 double POINT_ON_RESHOTKA::diam(double x, double y){
	 if (x> fabs(1-x)) x= fabs(1-x);// fabs - абсолютное значение для аргумента с плавоющей точкой
	 if (y>fabs(1-y)) y=fabs(1-y); // пишем в глобальные переменные?(x)
	 return (x>=y)?y:x;// возвращаем минимальное R ?
	  }
  
