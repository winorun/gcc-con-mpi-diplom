#include "reshotka.hpp"

POINT_ON_RESHOTKA::POINT_ON_RESHOTKA(int *argc, char** argv[],double x,double y){
MPI_Init (argc, argv);
MPI_Comm_rank (MPI_COMM_WORLD, &rank);
MPI_Comm_size (MPI_COMM_WORLD, &size);
srand(rank);
argcClass=argc;
argvClass=*argv;
porX=x;
porY=y;
flag=1;
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
		printf("rank - %i \t size - %i\n",rank,size);
		if(*argcClass>=2){
			printf("%i\t \t arg1 = %s\n",*argcClass, argvClass[1]);
			int N = atoi(argvClass[1]);
			if(!N)printf("number is not correctly\n");else printf("N=%i\n",N);
		}
	}
}
	
void POINT_ON_RESHOTKA::printResult()
{
	if(!rank){
		printf("N - %i\n",size);
		if(flag&0x1)printf("Presise solution:\t %f\n",u(porX,porY));
		if(flag&0x1)printf("Numerical solution:\t %f\n",U);
		if(flag&0x1)printf("Delta:\t\t\t %f\n",u(porX,porY)-U);
		if(flag&0x1)printf("Disp:\t\t\t %f\n",Disp);
		}
}

void POINT_ON_RESHOTKA::voidMain()
{
	double x=porX,y=porY;
	double S=0,S1=0;
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
//		printf("%i) \t %f\n",rank,S);
//	   U=S/size;
	   Disp=S*S;
		double inbuf[3],outbuf[3]={S,Disp,0};
		MPI_Reduce(outbuf, inbuf, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		if(!rank){
			U=inbuf[0]/size;
			Disp=inbuf[1]/size;
			//Disp=sqrt(fabs(Disp-U*U)/size);
			printResult();
		}
	}

 double POINT_ON_RESHOTKA::diam(double x, double y){
	 if (x> fabs(1-x)) x= fabs(1-x);// fabs - абсолютное значение для аргумента с плавоющей точкой
	 if (y>fabs(1-y)) y=fabs(1-y); // пишем в глобальные переменные?(x)
	 return (x>=y)?y:x;// возвращаем минимальное R ?
	  }
  
