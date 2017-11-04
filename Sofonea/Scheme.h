#include "Latticevelocity.h"

void computeSingleStepScheme(double*** fin, double*** F, double*** feq, 
							 int xNum, int yNum, double omega, double tau,
							 double hx, int* v_x, int* v_y);

void computePredictor(double*** fin, double*** fin1, double*** F, double*** feq, 
					  int* v_x, int* v_y,
					  double tau, double hx, double hy, double omega,
					  int xNum, int yNum);

void computeCorrector(double*** fin, double*** fin1,  double*** feq, 
					  int* v_x, int* v_y,
					  double tau, double hx, double hy, double omega,
					  int xNum, int yNum);

void computeAverageLayer(double*** fin1, double*** fin, int xNum, int yNum);