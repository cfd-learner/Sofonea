#include "Array.h"


double*** createArray(int basisNum, int xNum, int yNum)
{
	double*** f = new double**[basisNum];
	for (int i = 0; i < basisNum; i++)
	{
		f[i] = new double*[xNum];
		for (int j = 0; j < xNum; j++)
			f[i][j] = new double[yNum];
	}
	return f;

}

double** createArray(int xNum, int yNum)
{
	double** f = new double*[xNum];
	for (int i = 0; i < xNum; i++)
		f[i] = new double[yNum];
	return f;
}

void Copy(double*** dest, double*** source, int basisNum, int xNum, int yNum)
{
	 for(int i = 0; i < basisNum; i++)
       for(int j = 0; j < xNum; j++)
		   for(int k = 0; k < yNum; k++)
             dest[i][j][k] = source[i][j][k];
	   
}

void freeMemory(double*** arr, int basisNum, int xNum)
{
	for (int i = 0; i < basisNum; i++) {
		for (int j = 0; j < xNum; j++) {
			delete [] arr[i][j];
		}
		delete [] arr[i];
	}
	delete [] arr;
}

void freeMemory(double** arr, int xNum)
{
    for (int i = 0; i < xNum; i++) {
        delete [] arr[i];
    }
    delete [] arr;
}