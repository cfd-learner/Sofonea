#include "Scheme.h"

void computeSingleStepScheme(double*** fin, double*** F, double*** feq, int xNum, int yNum, double omega, double tau, double hx, int* v_x, int* v_y)
{
     for (int p = 0;p < BASIS; p++)
           for (int i = 1;i < xNum - 1; i++)
               for (int k = 1;k < yNum - 1; k++)
			       fin[p][i][k] = F[p][i][k]-(omega) *(F[p][i][k] - feq[p][i][k])
				       - tau / hx * (F[p][i][k] - F[p][i - v_x[p]][k - v_y[p]]);
}

void computePredictor(double*** fin, double*** fin1, double*** F, double*** feq, 
					  int* v_x, int* v_y,
					  double tau, double hx, double hy, double omega,
					  int xNum, int yNum)
{
      for (int p = 0; p < 9; p++)
           for (int i = 1; i < xNum - 1; i++)
               for (int k = 1; k < yNum - 1; k++)
			         fin[p][i][k] = fin1[p][i][k] 
			        - 2.0 * tau / hx * ( F[p][i][k] - F[p][i - v_x[p]][k - v_y[p]])
				     -2.0 * (omega) * (F[p][i][k] - feq[p][i][k]);

}

void computeCorrector(double*** fin, double*** fin1,  double*** feq, 
					  int* v_x, int* v_y,
					  double tau, double hx, double hy, double omega,
					  int xNum, int yNum)
{
	  for (int p = 0; p < BASIS; p++)
           for (int i = 1;i < xNum - 1; i++)
               for (int k = 1; k< yNum - 1; k++)
			         fin[p][i][k] = 0.5 * (fin1[p][i][k] + fin[p][i][k])
					 - tau / hx * (fin[p][i + v_x[p]][k + v_y[p]] - fin[p][i][k])
					 -(omega) * (fin[p][i][k] - feq[p][i][k]);

}

void computeAverageLayer(double*** fin1, double*** fin, int xNum, int yNum)
{
for (int p = 0; p < BASIS; p++)
         for (int i = 0;i < xNum; i++)
		     for (int j = 0;j < yNum; j++)
                fin1[p][i][j] = 0.5 * (fin[p][i][j] + fin1[p][i][j]);
}