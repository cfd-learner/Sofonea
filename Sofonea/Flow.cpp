#include "Flow.h"

void normalizeVelocity(double** Ux, double** Uy, double U0, int xNum, int yNum)
{
	for (int i = 0; i< xNum; i++)
	{
		for (int j = 0; j < yNum;j++)
		{
			Ux[i][j] = Ux[i][j] / U0;
			Uy[i][j] = Uy[i][j] / U0;
		}
	}
}

void computeVortex(double** vortex, double** Ux, double** Uy, double hx, int xNum, int yNum)
{
	for (int i = 1; i < xNum - 1; i++)
		for (int j = 1; j < yNum - 1; j++)
			vortex[i][j] = (0.5 / hx) * (Uy[i + 1][j] - Uy[i - 1][j]
		- Ux[i][j + 1] + Ux[i][j - 1]);

}

double** setInitialStreamFunction(int xNum, int yNum)
{
	double** psi = createArray(xNum, yNum);
	for (int i = 0; i < xNum; i++)
		for (int j = 0; j < yNum; j++)
			psi[i][j] = 0;
	return psi;
}

void computeStreamFunction(double** psi, double** vortex, double hx, int countPoisson, int xNum, int yNum)
{
	for (int k = 0; k < countPoisson; k++)
		for (int i = 1; i < xNum - 1; i++)
			for (int j = 1; j < yNum - 1; j++)
				psi[i][j] = (hx * hx) * 0.25 * vortex[i][j] 
			+ 0.25*(psi[i + 1][j] + psi[i - 1][j] + psi[i][j + 1] + psi[i][j - 1]);
}

double findMinimalStreamFunction(double** psi, int xNum, int yNum)
{
	double psiMin;
	double* P = new double[xNum];
	double* Pmin = new double[yNum];
	for (int i = 0; i< xNum; i++)
	{
		Copy(P, psi[i], yNum);
		Pmin[i] = BubbleSort(P, yNum);
	}
	psiMin = BubbleSort(Pmin, xNum);

	return psiMin;
}

// если размерности сетки совпадают
double* findCoordinatesMinStreamFunctions(double** psi, double psiMin, double* x, int xNum, int yNum)
{
	double* coord = new double[2];
	double** X = createArray(xNum, yNum);
	double** Y = createArray(xNum, yNum);
	// аккуратно
	for (int i = 0; i < xNum; i++)
	{
		Copy(X[i], x, yNum);
		Copy(Y[i], x,  yNum);
	}
	for(int i = 0; i < xNum; i++)
		for(int j = 0; j < yNum; j++)
			if (psi[i][j] == psiMin)
			{
				coord[0] = X[i][j];
				coord[1] = Y[i][j];
				break;
			} 
	freeMemory(X, xNum);
	freeMemory(Y, xNum);
	return coord;

}