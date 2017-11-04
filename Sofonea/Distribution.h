#include "Latticevelocity.h"

double* setConstantsForEquilibriumFunction();
void computeEqDistribution(double ***feq, double *W, double **rho, 
						   double **Ux, double**Uy, int *v_x, int *v_y, int xNum, int yNum);

void ApproximationForBGK_onTheUpperBoundary(double ***fin, double ***F,double ***feq, 
											double **rho, int *v_x, int *v_y, double **Ux, double **Uy, 
											double omega, double tau, double hx, double hy, int j, int n);

void ApproximationForBGK_onTheLeftBoundary(double ***fin, double ***F,double ***feq,
										   double **rho, int *v_x, int *v_y, double **Ux, double **Uy,
										   double omega, double tau, double hx, double hy, int i, int m);

void ApproximationForBGK_onTheRightBoundary(double ***fin, double ***F,double ***feq,
											double **rho, int *v_x, int *v_y, double **Ux, double **Uy,
											double omega, double tau, double hx, double hy, int i, int m);

void ApproximationForBGK_onTheLowerBoundary(double ***fin, double ***F,double ***feq,
											double **rho, int *v_x, int *v_y, double **Ux, double **Uy,
											double omega, double tau, double hx, double hy, int j, int n);

void setZouHeBoundaryCondRight(double*** fin, double** rho, double** Ux, double** Uy, int xNum, int yNum);
void setZouHeBoundaryCondLeft(double*** fin, double** rho, double** Ux, double** Uy, int xNum, int yNum);
void setZouHeBoundaryCondLower(double*** fin, double** rho, double** Ux, double** Uy, int xNum, int yNum);
void setZouHeBoundaryCondUp(double*** fin, double** rho, double** Ux, double** Uy, int xNum, int yNum);