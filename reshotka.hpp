#ifndef RESHOTKA
#define RESHOTKA

// === Принятый способ именования ===		\\
// Имена функций и переменых -- с маленькой\\
// приставкой -- тип\\
// Типы -- с большой \\
// Макросы и константы -- заглавными	\\
// 					=== *** ===								\\

enum Flag{
	NOT_Flag = 0,
	DINAMIC_ALG = 0x40,
	SPHERE_METOS = 0x80,
	G_EXISTS=0x04,
	U_EXISTS=0x08,
	BOUNDARY_EXISTS=0x10,
	PHI_FULL_EXISTS=0x20,
	LIB_INIT=1,
	LIB_MPI_CLOSE=2,
	LIB_ERROR=0x100
};

struct Rezult{
	double dRezSum;
	double dDisp;
	int iStatus;
};

struct Point{
	double dX,dY;
	Point();
	Point(double,double);
};

Flag libInit(int *, char***);
Flag libClose();

typedef double (*FunWith2Par) (double , double);
Flag setFunctions(FunWith2Par U,FunWith2Par G,FunWith2Par Phi_0,FunWith2Par Phi_1);
Flag setBoundary(int (*use) (double &, double &));

Rezult libRunComputing(Point startPoint,Flag & flag,int iN);
Rezult libRunComputing(Point startPoint,Flag & flag,double dDispersion);
#endif
