#include "Macrovelocity.h"

double** setInitialMacroVelocity(int xNum, int yNum)
{
	double** vel = createLayer(xNum, yNum);
	for (int i = 0; i < xNum; i++)
		for (int j = 0; j < yNum; j++)
			vel[i][j] = 0;
	return vel;
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