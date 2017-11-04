#define BASIS 9

void computeSingleStepScheme(double*** fin, double*** F, double*** feq, int xNum, int yNum, double omega, double tau, double hx, int* v_x, int* v_y)
{
     for (int p = 0;p < BASIS; p++)
           for (int i = 1;i < xNum - 1; i++)
               for (int k = 1;k < yNum - 1; k++)
			       fin[p][i][k] = F[p][i][k]-(omega) *(F[p][i][k] - feq[p][i][k])
				       - tau / hx * (F[p][i][k] - F[p][i - v_x[p]][k - v_y[p]]);
}