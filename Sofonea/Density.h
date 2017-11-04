double** setInitialDensity(int xNum, int yNum);
void computeDensity(double** den, double*** f, int xNum, int yNum);
void setBoundaryCondDensity(double ***fin,  double **rho,
							double **Ux, double**Uy, double U0, int xNum, int yNum);