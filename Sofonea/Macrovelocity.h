double** setInitialMacroVelocity(int xNum, int yNum);
void computeMacroVelocity(double** Ux, double** Uy, double** den, double*** f, int xNum, int yNum);
void setBoundaryCondMacroVelocity(double** Ux, double** Uy, double U0, int xNum, int yNum);