#include "Distribution.h"

double* setConstantsForEquilibriumFunction()
{
	double *W = new double[BASIS];
	W[0] = 4.0 / 9; W[1]= 1.0 / 9; W[2] = W[1];
	W[3] = W[1];    W[4] = W[1];   W[5] = 1.0 / 36;
	W[6] = W[5];    W[7] = W[5];   W[8] = W[5];

	return W;
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

//*********************************************************

void ApproximationForBGK_onTheUpperBoundary(double ***fin, double ***F,double ***feq,  double **rho, int *v_x, int *v_y, double **Ux, double **Uy, double omega, double tau, double hx, double hy, int j, int n)
{   

	for (int p = 0; p < BASIS; p++)
	{
		//Left(x)-left(y) approximation
		if ((p == 1) || (p == 5) || (p == 2))
		{
			for (int s = 1; s < n; s++) 
				fin[p][s][j]=(1 - (omega)) * F[p][s][j] 
			+ (omega) * feq[p][s][j] 
			- v_x[p] * tau / hx *(F[p][s][j] - F[p][s - 1][j]) 
				- v_y[p] * tau / hy *(F[p][s][j] - F[p][s][j - 1]);     
		}

		//Right(x)-left(y) approximation
		if ((p == 6) || (p == 3)) 
		{                            
			for (int s = 0; s < n - 1;s++) 
				fin[p][s][j] = (1-(omega)) * F[p][s][j] + (omega) * feq[p][s][j]
			- v_x[p] * tau / hx * (F[p][s + 1][j] - F[p][s][j])
				- v_y[p] * tau / hy * (F[p][s][j] - F[p][s][j - 1]);
		}
		if (p == 0)
		{
			for (int s = 0; s < n; s++)
				fin[p][s][j] = (1 - (omega)) * F[p][s][j] + (omega) * feq[p][s][j];
		}


	}

}


void ApproximationForBGK_onTheLeftBoundary(double ***fin, double ***F,double ***feq,  double **rho, int *v_x, int *v_y, double **Ux, double **Uy, double omega, double tau, double hx, double hy, int i, int m)
{   

	for (int p = 0; p < BASIS; p++)
	{
		//Right(x)-right(y) approximation
		if ((p == 7) || (p == 4))
		{
			for (int s = 0; s < m - 1; s++) 
				fin[p][i][s] = (1 - (omega)) * F[p][i][s] + (omega) * feq[p][i][s]
			- v_x[p] * tau / hx * (F[p][i + 1][s] - F[p][i][s])
				- v_y[p] * tau / hy * (F[p][i][s + 1] - F[p][i][s]);
		}


		//Right(x)-left(y) approximation
		if ((p == 2) || (p == 6)) 
		{ 
			for (int s = 1; s < m; s++) 
				fin[p][i][s] = (1 - (omega)) * F[p][i][s] + (omega) * feq[p][i][s]
			- v_x[p] * tau / hx * (F[p][i + 1][s] - F[p][i][s])
				- v_y[p] * tau / hy * (F[p][i][s] - F[p][i][s - 1]);
		}

		if ((p == 3) || (p == 0))
		{
			for (int s = 1; s < m; s++)
				fin[p][i][s]= (1 - (omega)) * F[p][i][s] + (omega) * feq[p][i][s]
			- v_x[p] * tau / hx * (F[p][i + 1][s] - F[p][i][s])
				- v_y[p] * tau / hy * (F[p][i][s] - F[p][i][s - 1]);
		}

	}
}


void ApproximationForBGK_onTheRightBoundary(double ***fin, double ***F,double ***feq,  double **rho, int *v_x, int *v_y, double **Ux, double **Uy, double omega, double tau, double hx, double hy, int i, int m)
{   
	for (int p = 0; p < BASIS; p++)
	{
		//Left(x)-right(y) approximation
		if ((p == 4) || (p == 8))
		{
			for (int s = 0; s < m - 1; s++)
				fin[p][i][s] = (1 - (omega)) * F[p][i][s] + (omega) * feq[p][i][s]
			- v_x[p] * tau / hx * (F[p][i][s] - F[p][i - 1][s])
				-v_y[p] * tau / hy * (F[p][i][s + 1] - F[p][i][s]);
		} 

		//Left(x)-left(y) approximation
		if ((p == 2) || (p == 5))
		{
			for (int s = 1; s < m; s++)
				fin[p][i][s] = (1 - (omega)) * F[p][i][s] + (omega) * feq[p][i][s]
			        -v_x[p] * tau / hx * (F[p][i][s] - F[p][i - 1][s])
			     	-v_y[p] * tau / hy * (F[p][i][s] - F[p][i][s - 1]);
		} 

		if (p == 1)
		{
			for (int s = 0; s < m; s++)
			{
				fin[p][i][s] = (1 - (omega)) * F[p][i][s] + (omega) * feq[p][i][s]
				   - v_x[p] * tau / hx * (F[p][i][s] - F[p][i - 1][s])
				   - v_y[p] * tau / hy * (F[p][i][s] - F[p][i][s - 1]);
			}
		} 			

		if (p == 0)
		{
			for(int s = 0; s < m; s++)
				fin[p][i][s] = (1 - (omega)) * F[p][i][s] + (omega) * feq[p][i][s];
		}


	}

}


void ApproximationForBGK_onTheLowerBoundary(double ***fin, double ***F,double ***feq, double **rho, int *v_x, int *v_y, double **Ux, double **Uy, double omega, double tau, double hx, double hy, int j, int n)
{   
	for (int p = 0; p < BASIS; p++)
	{
		//Left(x)-right(y) approximation
		if ((p == 1) || (p == 8))
		{
			for (int s = 1; s < n; s++)
				fin[p][s][j] = (1 - (omega)) * F[p][s][j] + (omega) * feq[p][s][j]
			        -v_x[p] * tau / hx * (F[p][s][j] - F[p][s - 1][j])
					-v_y[p] * tau / hy * (F[p][s][j + 1] - F[p][s][j]);
		} 

		if ((p == 4) || (p == 0))
		{
			for (int s = 1; s < n; s++)
				fin[p][s][j] = (1 - (omega)) * F[p][s][j] + (omega) * feq[p][s][j]
			          -v_x[p] * tau / hx * (F[p][s][j] - F[p][s - 1][j])
					  -v_y[p] * tau / hy * (F[p][s][j + 1] - F[p][s][j]);
		}

		//Right(x)-right(y) approximation
		if ((p == 3) || (p == 7))
		{
			for (int s = 0; s < n - 1; s++)
				fin[p][s][j] = (1 - (omega)) * F[p][s][j] + (omega) * feq[p][s][j]
			        -v_x[p] * tau / hx * (F[p][s + 1][j] - F[p][s][j])
			    	-v_y[p] * tau / hy * (F[p][s][j + 1] - F[p][s][j]);

		}


	}
}

//**********************************************************

void setZouHeBoundaryCondRight(double*** fin, double** rho, double** Ux, double** Uy, int xNum, int yNum)
{
	for (int i = 0; i < xNum; i++)
	{
		fin[4][i][yNum - 1] = fin[2][i][yNum - 1] - 2.0 / 3 * (rho[i][yNum - 1] * Uy[i][yNum-1]);
		fin[7][i][yNum - 1] = 0.5 * (fin[1][i][yNum - 1] - fin[3][i][yNum - 1])
			-0.5 * rho[i][yNum - 1] * Ux[i][yNum - 1]
		-1.0 / 6 * rho[i][yNum - 1] * Uy[i][yNum - 1] + fin[5][i][yNum - 1];
		fin[8][i][yNum - 1] = 0.5 * (fin[3][i][yNum - 1] - fin[1][i][yNum - 1]) 
			+ 0.5 * rho[i][yNum - 1] * Ux[i][yNum - 1]
		- 1.0 / 6 * rho[i][yNum - 1] * Uy[i][yNum -1] + fin[6][i][yNum - 1];
	}
}

void setZouHeBoundaryCondLeft(double*** fin, double** rho, double** Ux, double** Uy, int xNum, int yNum)
{
	for (int i = 0;i < xNum; i++)
	{
		fin[2][i][0] = fin[4][i][0] + 2.0 / 3 * (rho[i][0] * Uy[i][0]);
		fin[5][i][0] = fin[7][i][0] - 0.5 * (fin[1][i][0]-fin[3][i][0])
			+ 0.5 * rho[i][0] * Ux[i][0] + 1.0 / 6 * rho[i][0] * Uy[i][0];
		fin[6][i][0] = fin[8][i][0] + 0.5 * (fin[1][i][0] - fin[3][i][0])
			- 0.5 * rho[i][0] * Ux[i][0] + 1.0/ 6 * rho[i][0] * Uy[i][0];
	}
}

void setZouHeBoundaryCondLower(double*** fin, double** rho, double** Ux, double** Uy, int xNum, int yNum)
{
	for (int i = 0; i < xNum; i++)
	{
		fin[1][0][i] = fin[3][0][i] + 2.0 / 3 * rho[0][i] * Ux[0][i];
		fin[5][0][i] = fin[7][0][i] + 0.5 * (fin[4][0][i] - fin[2][0][i])
			+ 0.5 * rho[0][i] * Uy[0][i] + 1.0/6 *rho[0][i] * Ux[0][i];
		fin[8][0][i] = fin[6][0][i] + 0.5 * (fin[2][0][i] - fin[4][0][i])
			- 0.5 * rho[0][i] * Uy[0][i] + 1.0 / 6 *rho[0][i] * Ux[0][i];
	}
}

void setZouHeBoundaryCondUp(double*** fin, double** rho, double** Ux, double** Uy, int xNum, int yNum)
{
	for (int i = 0;i < xNum; i++)
	{
		fin[3][xNum - 1][i] = fin[1][xNum - 1][i]
		- 2.0/ 3 *rho[xNum - 1][i] * Ux[xNum - 1][i];
		fin[6][xNum - 1][i] = fin[8][xNum - 1][i] + 0.5 * (fin[4][xNum-1][i]
		-fin[2][xNum - 1][i]) + 0.5 * rho[xNum - 1][i] * Uy[xNum - 1][i]
		- 1.0 / 6 * rho[xNum - 1][i] * Ux[xNum - 1][i];
		fin[7][xNum - 1][i] = fin[5][xNum - 1][i] + 0.5 * (fin[2][xNum - 1][i] 
		- fin[4][xNum - 1][i]) - 0.5 * rho[xNum - 1][i] * Uy[xNum - 1][i]
		- 1.0 / 6 * rho[xNum - 1][i] * Ux[xNum - 1][i];
	}
}