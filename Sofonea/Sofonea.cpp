
//Scheme with central differences realisation. D2Q9 lattice.
//Scheme is based on divergent form of Boltzmann equation. 
//Modelling of the flow in rectangular cavern.
//Programmed by Gerasim Krivovichev

#include <cstdlib>
#include <iostream>
#include <math.h>
#include <conio.h>

#define BASIS 9

using namespace std;



double* setConstantsForEquilibriumFunction()
{
	double *W = new double[BASIS];
	W[0] = 4.0 / 9; W[1]= 1.0 / 9; W[2] = W[1];
	W[3] = W[1];    W[4] = W[1];   W[5] = 1.0 / 36;
	W[6] = W[5];    W[7] = W[5];   W[8] = W[5];

	return W;
}

int* setLatticeVelocityX()
{
	int *v_x = new int[BASIS];
	v_x[0] = 0;  v_x[1] = 1;  v_x[2] = 0;
	v_x[3] = -1; v_x[4] = 0;  v_x[5] = 1;
	v_x[6] = -1; v_x[7] = -1; v_x[8] = 1;

	return v_x;

}

int* setLatticeVelocityY()
{
	int* v_y = new int [BASIS];
	v_y[0] = 0; v_y[1] = 0;  v_y[2] = 1;
	v_y[3] = 0; v_y[4] = -1; v_y[5] = 1;
	v_y[6] = 1; v_y[7] = -1; v_y[8] = -1;
	return v_y;
}

double*** createDistributionFunction(int xNum, int yNum)
{
	double*** f = new double**[BASIS];
	for (int i = 0; i < BASIS; i++)
	{
		f[i] = new double*[xNum];
		for (int j = 0; j < xNum; j++)
			f[i][j] = new double[yNum];
	}
	return f;

}

double** createLayer(int xNum, int yNum)
{
	double** f = new double*[xNum];
	for (int i = 0; i < xNum; i++)
		f[i] = new double[yNum];
	return f;
}

void freeMemory(double*** arr, int xNum, int yNum)
{
	for (int i = 0; i < BASIS; i++) {
		for (int j = 0; j < xNum; j++) {
			delete [] arr[i][j];
		}
		delete [] arr[i];
	}
	delete [] arr;
}

void freeMemory(double** arr, int xNum, int yNum)
{
    for (int i = 0; i < xNum; i++) {
        delete [] arr[i];
    }
    delete [] arr;
}



int main(int argc, char *argv[])
{
	//что-то вынести в константы
	int xCount = 150, yCount = 150, timeCount = 800000;
	double timeLength = 400;
	double xStep, yStep, timeStep;
	double meanFreeTime, meanFreePath;
	double viscosity;
	double ReynoldsNumber = 400;
	double CourantNumber;
	double cavernLength = 1.0, cavernHeight = 1.0;
	double velocityUpBoundary = 0.1;
	double omegaCoeff;
	double s1;


	xStep = cavernLength /(xCount - 1);
	yStep = cavernHeight /(yCount - 1);
	timeStep = timeLength / (timeLength - 1);

	s1 = timeLength / xStep + 1;
	meanFreeTime = timeLength / (s1 - 1);

	viscosity = velocityUpBoundary * cavernLength / ReynoldsNumber;
	omegaCoeff = timeStep / (3 * viscosity);
	//*****************************************************
	//Lattice characteristics

	double *W = setConstantsForEquilibriumFunction();
	int *basisVx = setLatticeVelocityX();
	int *basisVy = setLatticeVelocityY();

	// massives for distribution function
	double*** fin = createDistributionFunction(xCount, yCount);
	double*** fin1 = createDistributionFunction(xCount, yCount);
	double*** F = createDistributionFunction(xCount, yCount);
	double*** feq = createDistributionFunction(xCount, yCount);
	double*** feq1 = createDistributionFunction(xCount, yCount);

	// macrocharacteristics
	double** density = createLayer(xCount, yCount);
	double** macroVelocityX = createLayer(xCount, yCount);
	double** macroVelocityY = createLayer(xCount, yCount);




	// очищение памяти
	delete[] W;
	delete[] basisVx;
	delete[] basisVy;
	freeMemory(fin, xCount, yCount);
	freeMemory(fin1, xCount, yCount);
	freeMemory(F, xCount, yCount);
	freeMemory(feq, xCount, yCount);
	freeMemory(feq1, xCount, yCount);
	freeMemory(density, xCount, xCount);
	freeMemory(macroVelocityX, xCount, yCount);
	freeMemory(macroVelocityY, xCount, yCount);


	system("PAUSE");
	return 0;
}

