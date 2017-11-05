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

void Copy(double** dest, double** source, int n, int m)
{
	 for(int i = 0; i < n; i++)
       for(int j = 0; j < m; j++)
             dest[i][j] = source[i][j];
	   
}

void Copy(double* dest, double* source, int n)
{
	 for(int i = 0; i < n; i++)
             dest[i] = source[i];
	   
}


// поиск минимального элемента
double BubbleSort(double *a, int n)
{
       double t;
       for (int i = n - 1; i >= 0; i--)
           for (j = 0; j < i; j++)
               if (a[j] > a[j+1])
               {
                   t = a[j];
                   a[j] = a[j + 1];
                   a[j + 1] = t;            
               }
       return a[0];
};

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

void freeMemory(double* arr)
{
    delete [] arr;
}

void freeMemory(int* arr)
{
    delete [] arr;
}