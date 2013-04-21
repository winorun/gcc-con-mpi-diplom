#include "reshotka.hpp"
#include <unistd.h>

POINT_ON_RESHOTKA::POINT_ON_RESHOTKA(int *argc, char** argv[],double x,double y){
MPI_Init (argc, argv);
MPI_Comm_rank (MPI_COMM_WORLD, &rank);
MPI_Comm_size (MPI_COMM_WORLD, &size);
srand(rank);
argcClass=argc;
argvClass=*argv;
porX=x;
porY=y;
flag=0x1;
for(int rez=getopt(*argc,*argv,"n:f");rez!=-1;rez=getopt(*argc,*argv,"n:f"))
	switch (rez){
		case 'n': porColDrave = atoi(optarg);if(!porColDrave)porColDrave=1000;break;
		case 'f': flag=flag | 0x100;break;
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
void POINT_ON_RESHOTKA::printDebag()
{

	if(!rank){
		printf("rank - %i \t N = %i\n",rank,N);
		printf("0x%x\n",flag);
//		printf("rank - %i \t size - %i\n",rank,size);
		printf("N=%i\n",porColDrave);
	}
}
	
void POINT_ON_RESHOTKA::printResult(int N){//времянка
	if(!rank){
		if(flag&0x100){
			printf("N - %i\n",N);
			if(flag&0x1)printf("Presise solution:\t %f\n",u(porX,porY));
			if(flag&0x1)printf("Numerical solution:\t %f\n",U);
			if(flag&0x1)printf("Delta:\t\t\t %f\n",u(porX,porY)-U);
			if(flag&0x1)printf("Disp:\t\t\t %f\n",Disp);
		}else{
			if(flag&0x200)printf("PS\t\t NS\t\t Delta\t\t Disp\n");
			printf("%f\t %f\t %f\t %f\n",u(porX,porY),U,u(porX,porY)-U,Disp);
			};
}}

void POINT_ON_RESHOTKA::voidMain()
{
	double S=0;
	N= porColDrave/(size-1);
	if(!rank){// менеджер распределения
		N = porColDrave%(size-1);
		}
for(int i=N;i;i--){
	double x=porX,y=porY;
	double S1=0;
		while (boundary(x,y)){
		   double d=diam(x,y);
		   
		   double alpha= DRAND;// шаг ?
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
//		if(!rank) printf("%i \t %f\n",i,S);
//		if(!rank) printf("%i \t %f\n",N-i+1,S/(N-i+1));
		
//	   U=S/size;
	   Disp+=(S*S)/porColDrave;
		//~ if(!rank) printf("%i \t %f\n",i,Disp);
		//~ if(!rank) printf("%i \t %f\n\n",N-i+1,Disp/(N-i+1));
	   //незабыть разделить на количество узлов
   }
		double inbuf[4],outbuf[4]={S/porColDrave,Disp,0};
		MPI_Reduce(outbuf, inbuf, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		int porColDraveAll;
		MPI_Reduce(&N, &porColDraveAll, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		if(!rank){
			U=inbuf[0];
			Disp=inbuf[1];
			//N=inbuf[2];
			Disp=sqrt(fabs(Disp-U*U)/porColDrave);
			printResult(porColDraveAll);
		}
	}

 double POINT_ON_RESHOTKA::diam(double x, double y){
	 if (x> fabs(1-x)) x= fabs(1-x);// fabs - абсолютное значение для аргумента с плавоющей точкой
	 if (y>fabs(1-y)) y=fabs(1-y); // пишем в глобальные переменные?(x)
	 return (x>=y)?y:x;// возвращаем минимальное R ?
	  }
  
