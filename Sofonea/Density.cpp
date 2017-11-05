#include "Density.h"

double** setInitialDensity(int xNum, int yNum)
{
	double** den = createArray(xNum, yNum);
	for (int i = 0; i < xNum; i++)
		for (int j = 0; j < yNum; j++)
			den[i][j] = 1;
	return den;
}

void computeDensity(double** den, double*** f, int xNum, int yNum)
{
	 for (int i = 1; i < xNum - 1; i++)
           for (int j = 1 ; j < yNum - 1; j++)
               den[i][j] = f[0][i][j] + f[1][i][j] + f[2][i][j]
		             + f[3][i][j] + f[4][i][j] + f[5][i][j]
					 + f[6][i][j] + f[7][i][j] + f[8][i][j];

}

void setBoundaryCondDensity(double ***fin,  double **rho, double **Ux, double**Uy, double U0, int xNum, int yNum)
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
