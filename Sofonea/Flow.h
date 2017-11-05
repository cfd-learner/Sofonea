#include "Array.h"

void normalizeVelocity(double** Ux, double** Uy, double U0, int xNum, int yNum);
void computeVortex(double** vortex, double** Ux, double** Uy, double hx, int xNum, int yNum);
double** setInitialStreamFunction(int xNum, int yNum);
void computeStreamFunction(double** psi, double** vortex, double hx, int countPoisson, int xNum, int yNum);
double findMinimalStreamFunction(double** psi, int xNum, int yNum);
double* findCoordinatesMinStreamFunctions(double** psi, double psiMin, double* x, int xNum, int yNum);
