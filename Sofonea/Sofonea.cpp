
//Scheme with central differences realisation. D2Q9 lattice.
//Scheme is based on divergent form of Boltzmann equation. 
//Modelling of the flow in rectangular cavern.
//Programmed by Gerasim Krivovichev

#include <cstdlib>
#include <iostream>
#include <math.h>
#include <conio.h>
#include <string.h>

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


double** setInitialDensity(int xNum, int yNum)
{
	double** den = createLayer(xNum, yNum);
	for (int i = 0; i < xNum; i++)
		for (int j = 0; j < yNum; j++)
			den[i][j] = 1;
	return den;
}

double** setInitialMacroVelocity(int xNum, int yNum)
{
	double** vel = createLayer(xNum, yNum);
	for (int i = 0; i < xNum; i++)
		for (int j = 0; j < yNum; j++)
			vel[i][j] = 0;
	return vel;
}

void computeEqDistribution(double ***feq, double *W, double **rho, double **Ux, double**Uy, int *v_x, int *v_y, int xNum, int yNum)
{
     double Sp;
       //variable for scalar production of lattice vector on macrospic velocity vector      
     for (int k = 0; k < BASIS; k++)
        for(int i = 0; i < xNum; i++)
              for(int j = 0; j < yNum; j++)
                {
                    Sp = v_x[k] * Ux[i][j] + v_y[k] * Uy[i][j];            
                    feq[k][i][j]=W[k] * rho[i][j] * (1 + 3 *Sp + 4.5 * Sp * Sp - 1.5 * (Ux[i][j] * Ux[i][j] + Uy[i][j] * Uy[i][j]));              
                }           
}

void Copy(double*** dest, double*** source, int xNum, int yNum)
{
	 for(int i = 0; i < BASIS; i++)
       for(int j = 0; j < xNum; j++)
		   for(int k = 0; k < yNum; k++)
             dest[i][j][k] = source[i][j][k];
	   
}

void computeDensity(double** den, double*** f, int xNum, int yNum)
{
	 for (int i = 1; i < xNum - 1; i++)
           for (int j = 1 ; j < yNum - 1; j++)
               den[i][j] = f[0][i][j] + f[1][i][j] + f[2][i][j]
		             + f[3][i][j] + f[4][i][j] + f[5][i][j]
					 + f[6][i][j] + f[7][i][j] + f[8][i][j];

}

void computeMacroVelocity(double** Ux, double** Uy, double** den, double*** f, int xNum, int yNum)
{
	 for (int i = 1; i < xNum - 1; i++)
           for (int j = 1 ; j < yNum - 1; j++)
		   {
			   Ux[i][j] = (f[1][i][j] + f[5][i][j] + f[8][i][j] 
			         - f[3][i][j] - f[6][i][j] - f[7][i][j]) / den[i][j];
               Uy[i][j] = (f[2][i][j] + f[5][i][j] + f[6][i][j] 
			         - f[4][i][j] - f[7][i][j] - f[8][i][j]) / den[i][j];
		   }

}

void setBoundaryCondMacroVelocity(double** Ux, double** Uy, double U0, int xNum, int yNum)
{
       
        for (int i = 0; i < xNum; i++)
       {
		   //Upper boundary
           Ux[i][yNum - 1] = U0;
           Uy[i][yNum - 1] = 0;
		   //lower boundary
		   Ux[i][0] = 0;
           Uy[i][0] = 0;
	   }
       
       for (int j = 0; j < yNum; j++)
       {
		   //left boundary
           Ux[0][j] = 0;
           Uy[0][j] = 0;
		   //right boundary
		   Ux[xNum - 1][j] = 0;
           Uy[xNum - 1][j] = 0;
	   }
           	   
}

void setBoundaryConditionsForRho(double ***fin,  double **rho, double **Ux, double**Uy, double U0, int xNum, int yNum)
{   
	  for (int i = 0; i < xNum; i++)
       {
		   // lower
           rho[i][0] =(fin[0][i][0] + fin[1][i][0] + fin[3][i][0] +
			   2.0 * (fin[4][i][0] + fin[7][i][0] + fin[8][i][0])) / (1.0 - Uy[i][0]);
		   // upper
           rho[i][yNum - 1] = (fin[0][i][yNum - 1] + fin[1][i][yNum - 1] + fin[3][i][yNum - 1]+ 
			   2.0 * (fin[4][i][yNum - 1] + fin[7][i][yNum - 1] + fin[8][i][yNum - 1])) / (1.0 + Uy[i][yNum - 1]);
	   }

       for (int j = 0; j < yNum; j++)
       {
		    //right
		    rho[xNum - 1][j] = (fin[0][xNum - 1][j] + fin[2][xNum - 1][j] + fin[4][xNum - 1][j]
			    + 2.0 * (fin[1][xNum - 1][j] + fin[5][xNum - 1][j] + fin[8][xNum - 1][j])) / (1.0+Ux[xNum - 1][j]);
            // left
            rho[0][j]=(fin[0][0][j] + fin[2][0][j] + fin[4][0][j]
			    + 2.0 * (fin[3][0][j] + fin[6][0][j] + fin[7][0][j])) / (1.0-Ux[0][j]);
	   
	   }
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
	double** density = setInitialDensity(xCount, yCount);
	double** macroVelocityX = setInitialMacroVelocity(xCount, yCount);
	double** macroVelocityY = setInitialMacroVelocity(xCount, yCount);

	computeEqDistribution(fin, W, density,
		macroVelocityX, macroVelocityY, basisVx, basisVy, xCount, yCount);
	//memcpy(F, fin, (sizeof(double) * BASIS * xCount * yCount));
	//memcpy(fin1, fin, (sizeof(double) * BASIS * xCount * yCount));
	Copy(F, fin, xCount, yCount);
	Copy(fin1, fin, xCount, yCount);


	// разгон
	for (int p = 0; p < 2; p++)
	{
		Copy(F, fin, xCount, yCount);
	    Copy(fin1, fin, xCount, yCount);
		computeDensity(density, fin, xCount, yCount);
		computeMacroVelocity(macroVelocityX, macroVelocityY, density, fin, xCount, yCount);
		setBoundaryCondMacroVelocity(macroVelocityX, macroVelocityY, velocityUpBoundary, xCount, yCount);



	}


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

