#define BASIS 9
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