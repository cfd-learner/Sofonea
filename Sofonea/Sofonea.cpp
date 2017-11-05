
//Scheme with central differences realisation. D2Q9 lattice.
//Scheme is based on divergent form of Boltzmann equation. 
//Modelling of the flow in rectangular cavern.
//Programmed by Gerasim Krivovichev
//Modified by Maria Maschinskaya

#include <cstdlib>
#include <iostream>
#include <math.h>
#include <conio.h>
#include <string.h>

#include "Array.h"
#include "Density.h"
#include "Distribution.h"
#include "Macrovelocity.h"
#include "Scheme.h"
#include "Latticevelocity.h"

using namespace std;

int main(int argc, char *argv[])
{
	//что-то вынести в константы
	int xCount = 150, yCount = 150, timeCount = 8;
	double timeLength = 400;
	double xStep, yStep, timeStep;
	double meanFreeTime, meanFreePath;
	double viscosity;
	double ReynoldsNumber = 400;
	double CourantNumber;
	double cavernLength = 1.0, cavernHeight = 1.0;
	double velocityUpBoundary = 0.1;
	double omegaCoeff;
	double s1;


	xStep = cavernLength /(xCount - 1);
	yStep = cavernHeight /(yCount - 1);
	timeStep = timeLength / (timeLength - 1);

	s1 = timeLength / xStep + 1;
	meanFreeTime = timeLength / (s1 - 1);

	viscosity = velocityUpBoundary * cavernLength / ReynoldsNumber;
	omegaCoeff = timeStep / (3 * viscosity);
	//*****************************************************
	//Lattice characteristics

	double *W = setConstantsForEquilibriumFunction();
	int *basisVx = setLatticeVelocityX();
	int *basisVy = setLatticeVelocityY();

	// massives for distribution function
	double*** fin = createArray(BASIS, xCount, yCount);
	double*** fin1 = createArray(BASIS, xCount, yCount);
	double*** F = createArray(BASIS, xCount, yCount);
	double*** feq = createArray(BASIS, xCount, yCount);
	double*** feq1 = createArray(BASIS, xCount, yCount);

	// macrocharacteristics
	double** density = setInitialDensity(xCount, yCount);
	double** macroVelocityX = setInitialMacroVelocity(xCount, yCount);
	double** macroVelocityY = setInitialMacroVelocity(xCount, yCount);

	computeEqDistribution(fin, W, density,
		macroVelocityX, macroVelocityY, basisVx, basisVy, xCount, yCount);
	//memcpy(F, fin, (sizeof(double) * BASIS * xCount * yCount));
	//memcpy(fin1, fin, (sizeof(double) * BASIS * xCount * yCount));
	Copy(F, fin, BASIS, xCount, yCount);
	Copy(fin1, fin,  BASIS, xCount, yCount);


	// разгон
	for (int q = 0; q < 2; q++)
	{
		Copy(F, fin, BASIS, xCount, yCount);
		Copy(fin1, fin, BASIS, xCount, yCount);
		computeDensity(density, fin, xCount, yCount);
		computeMacroVelocity(macroVelocityX, macroVelocityY, density, fin, xCount, yCount);
		setBoundaryCondMacroVelocity(macroVelocityX, macroVelocityY, velocityUpBoundary, xCount, yCount);
		setBoundaryCondDensity(fin, density, macroVelocityX, macroVelocityY, velocityUpBoundary, xCount, yCount);
		computeEqDistribution(feq, W, density, macroVelocityX, macroVelocityY, basisVx, basisVy,  xCount, yCount);
		computeSingleStepScheme(fin, F, feq, xCount, yCount, omegaCoeff, timeStep, xStep, basisVx, basisVy);

		ApproximationForBGK_onTheUpperBoundary(fin, F, feq, density, basisVx, basisVy, macroVelocityX, macroVelocityY,
			omegaCoeff, timeStep, xStep, yStep, yCount - 1, xCount);  //
		setZouHeBoundaryCondRight(fin, density, macroVelocityX, macroVelocityY, xCount, yCount);
		ApproximationForBGK_onTheLowerBoundary(fin, F, feq, density, basisVx, basisVy, macroVelocityX, macroVelocityY,
			omegaCoeff, timeStep, xStep, yStep, 0, xCount);   //
		setZouHeBoundaryCondLeft(fin, density, macroVelocityX, macroVelocityY, xCount, yCount);
		ApproximationForBGK_onTheLeftBoundary(fin, F, feq, density, basisVx, basisVy, macroVelocityX, macroVelocityY,
			omegaCoeff, timeStep, xStep, yStep, 0, yCount); //
		setZouHeBoundaryCondLower(fin, density, macroVelocityX, macroVelocityY, xCount, yCount);
		ApproximationForBGK_onTheRightBoundary(fin, F,feq, density, basisVx, basisVy, macroVelocityX, macroVelocityY,
			omegaCoeff, timeStep, xStep, yStep, xCount - 1, yCount); //
		setZouHeBoundaryCondUp(fin, density, macroVelocityX, macroVelocityY, xCount, yCount);
	}

	// MAIN COMPUTATIONAL CYCLE
	for (int time = 0; time < timeCount; time++)
	{
		Copy(F, fin,  BASIS, xCount, yCount);
		computeDensity(density, fin, xCount, yCount);
		computeMacroVelocity(macroVelocityX, macroVelocityY, density, fin, xCount, yCount);
		setBoundaryCondMacroVelocity(macroVelocityX, macroVelocityY, velocityUpBoundary, xCount, yCount);
		setBoundaryCondDensity(fin, density, macroVelocityX, macroVelocityY, velocityUpBoundary, xCount, yCount);
		computeEqDistribution(feq, W, density, macroVelocityX, macroVelocityY, basisVx, basisVy,
			xCount, yCount);

		// PREDICTOR
		computePredictor(fin, fin1, F, feq, basisVx, basisVy,
			timeStep, xStep, yStep, omegaCoeff,
			xCount, yCount);

		ApproximationForBGK_onTheUpperBoundary(fin, F, feq, density, basisVx, basisVy, macroVelocityX, macroVelocityY,
			omegaCoeff, timeStep, xStep, yStep, yCount - 1, xCount);
		setZouHeBoundaryCondRight(fin, density, macroVelocityX, macroVelocityY, xCount, yCount);
		ApproximationForBGK_onTheLowerBoundary(fin, F, feq, density, basisVx, basisVy, macroVelocityX, macroVelocityY,
			omegaCoeff, timeStep, xStep, yStep, 0, xCount);
		setZouHeBoundaryCondLeft(fin, density, macroVelocityX, macroVelocityY, xCount, yCount);
		ApproximationForBGK_onTheLeftBoundary(fin, F, feq, density, basisVx, basisVy, macroVelocityX, macroVelocityY,
			omegaCoeff, timeStep, xStep, yStep, 0, yCount);
		setZouHeBoundaryCondLower(fin, density, macroVelocityX, macroVelocityY, xCount, yCount);
		ApproximationForBGK_onTheRightBoundary(fin, F, feq, density, basisVx, basisVy, macroVelocityX, macroVelocityY,
			omegaCoeff, timeStep, xStep, yStep, xCount - 1, yCount);
		setZouHeBoundaryCondUp(fin, density, macroVelocityX, macroVelocityY, xCount, yCount);


		// CORRECTOR
		computeDensity(density, fin, xCount, yCount);
		computeMacroVelocity(macroVelocityX, macroVelocityY, density, fin, xCount, yCount);
		setBoundaryCondMacroVelocity(macroVelocityX, macroVelocityY, velocityUpBoundary, xCount, yCount);
		setBoundaryCondDensity(fin, density, macroVelocityX, macroVelocityY, velocityUpBoundary, xCount, yCount);
		computeEqDistribution(feq, W, density, macroVelocityX, macroVelocityY, basisVx, basisVy, xCount, yCount);
		computeCorrector(fin, fin1, feq, basisVx, basisVy, timeStep, xStep, yStep, omegaCoeff, xCount, yCount);
		computeAverageLayer(fin1, fin, xCount, yCount);

		cout << "time = " << time << endl;
	}


	computeDensity(density, fin, xCount, yCount);
	computeMacroVelocity(macroVelocityX, macroVelocityY, density, fin, xCount, yCount);
	setBoundaryCondMacroVelocity(macroVelocityX, macroVelocityY, velocityUpBoundary, xCount, yCount);
	setBoundaryCondDensity(fin, density, macroVelocityX, macroVelocityY, velocityUpBoundary, xCount, yCount);

	// очищение памяти
	delete[] W;
	delete[] basisVx;
	delete[] basisVy;
	freeMemory(fin,  BASIS, xCount);
	freeMemory(fin1,  BASIS, xCount);
	freeMemory(F,  BASIS, xCount);
	freeMemory(feq,  BASIS, xCount);
	freeMemory(feq1,  BASIS, xCount);
	freeMemory(density, xCount);
	freeMemory(macroVelocityX, xCount);
	freeMemory(macroVelocityY, xCount);


	system("PAUSE");
	return 0;
}

